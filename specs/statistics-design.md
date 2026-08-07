# Statistics module design

Design note for `include/dd/exact/statistics/` and `src/dd/exact/statistics/`.
Written 2026-08-07, when the module was added.

## Purpose

Give `DwPackage` the same instrumentation surface MQT Core's `dd::Package` has, so that the two
QRAT backends' reports can be read side by side: unique-table hit ratios, memory-manager usage,
compute-table hit rates, garbage-collection counts, and peak/active memory, emitted as JSON.

Before this, exact-dd's entire metric surface was `vNodeCount()` / `mNodeCount()` plus two
function-local `vFreed`/`mFreed` counters in `garbageCollect()` that were computed and thrown away.
`ExactDDBench` printed wall time and node counts and nothing else. There was no way to answer, for
a Grover run at 20+ qubits, *where* the cost goes — unique-table thrash, compute-table miss rate,
GC churn, or chunk over-allocation.

The reference implementation is MQT Core's `include/mqt-core/dd/statistics/` (five headers, four
sources). Field names and JSON keys are deliberately identical to it. Three structural things could
not be ported literally; those are §3.

## 1. File map

| File | Contents |
|---|---|
| `statistics/StatisticsConfig.hpp` | `kStatisticsEnabled` — the compile-time gate |
| `statistics/Statistics.hpp` / `.cpp` | Abstract base: `reset()`, `json()`, `toString()`, `operator<<` |
| `statistics/TableStatistics.hpp` / `.cpp` | Generic hash-table counters; used for the two compute caches |
| `statistics/UniqueTableStatistics.hpp` / `.cpp` | Adds active-entry and gc counters |
| `statistics/MemoryManagerStatistics.hpp` / `.cpp` | Templated pooled-allocator counters |
| `statistics/PackageStatistics.hpp` / `.cpp` | The aggregate, plus `DwPackage::statistics()` |

Hierarchy:

```
Statistics                       (json / toString / operator<< / reset)
 +- TableStatistics              (entrySize, numBuckets, numEntries, peakNumEntries,
 |   |                            collisions, hits, lookups, inserts)
 |   +- UniqueTableStatistics    (+ numActiveEntries, peakNumActiveEntries, gcRuns)
 +- MemoryManagerStatistics<T>   (numAllocations, numAllocated, numUsed,
 |                                numAvailableForReuse, + two peaks)
 +- PackageStatistics            (the six records above, aggregated)
```

Ownership is by-value composition, as in MQT Core. `MemoryManager<T>` holds its own
`MemoryManagerStatistics<T>` and exposes it via `statistics()`. `DwPackage` holds four records
(`vUniqueStats_`, `mUniqueStats_`, `mvCacheStats_`, `mmCacheStats_`) and folds all six into a
`PackageStatistics` on demand. Nothing is heap-allocated or globally registered.

## 2. Why the module exists in this shape

MQT Core's design has two layers worth keeping and one worth discarding.

Keep: **the struct hierarchy and the JSON contract.** `json()` is the only thing subclasses
override; `toString()` is `json().dump(2)` and `operator<<` is `toString()`. That means one
serialization path, no format drift, and a report a benchmark harness can consume directly.

Keep: **the `"unused"` sentinel.** A table that was never probed (`lookups == 0`), or a manager
that never handed out an entry (`peakNumUsed == 0`), serializes as the JSON *string* `"unused"`
rather than an object of zeros. It keeps an idle 12-qubit report readable. The cost is that the
reported type is polymorphic — a consumer must handle string-or-object per section.

Discard: **the free-function aggregation.** See D1.

## 3. Design decisions

### D1. `statistics()` is a const member returning a snapshot, not MQT's free functions

MQT Core aggregates via header-only free functions — `getStatistics(package)`,
`printStatistics(package)` — that reach directly into `Package`'s members. That works only because
every member of `Package` is public.

Every member of `DwPackage` is private (`vUnique_`, `mUnique_`, `mvCache_`, `mmCache_`, `vMemory_`,
`mMemory_`). Making them public to enable a free function would trade real encapsulation for a
stylistic match; a `friend` declaration would be the same trade with extra indirection. Asking the
package for its own report is also how the rest of `DwPackage`'s introspection surface already
reads (`vNodeCount()`, `vectorToString()`, `printVectorDiagram()`), so `statistics()` is the
locally idiomatic choice.

Returning a **snapshot by value** is what lets the method be `const`. The map-derived fields (D2)
have to be refreshed at report time; refreshing them on the returned copy means nothing in the
package is mutated and no member needs to be `mutable`. The copy is six small structs of
`std::size_t`; it is a report-time cost only.

`DwPackage::statistics()` is defined in `statistics/PackageStatistics.cpp`, not `DwPackage.cpp`. A
member's definition may live in any translation unit, and putting it there keeps the full
`<nlohmann/json.hpp>` out of `DwPackage.cpp`'s include chain.

### D2. Map-derived fields are snapshotted at report time, never tracked

`numEntries`, `numBuckets` and `collisions` are properties of a map's *current state*, not of the
history of operations on it. `TableStatistics::snapshot(map)` re-reads all three from the map when
the report is built.

The alternative — MQT Core's approach of self-incrementing `numEntries` in `trackInsert()` — does
not survive contact with this codebase. `garbageCollect()` erases from `vUnique_`/`mUnique_`
directly, and `.clear()`s both compute caches, so a self-incremented counter would need a matching
decrement at every erase site and would silently drift the moment a new one was added. The map is
the single source of truth; `numEntries` is only a cache of it, refreshed both in `trackInsert()`
and again in `snapshot()`.

`peakNumEntries` is genuinely historical and cannot be recovered from the map, so it *is*
maintained incrementally.

`snapshot()` is O(numBuckets + numEntries) because of the collision walk. It runs once per report
and must not migrate anywhere else.

### D3. `collisions` is a current-state measure, not MQT Core's cumulative chain-step count

**This is the one number that is not comparable between the two backends despite sharing a name.**

MQT Core's tables are fixed-`NBUCKET` chained tables that never rehash, and its `collisions` counts
chain steps walked during unsuccessful comparisons in `searchTable()` — a cumulative event count.

`DwPackage`'s tables are `std::unordered_map`s. They expose no collision counter, and they rehash
as they grow, dissolving and re-forming every chain. So `collisions` is defined here as

```
sum over buckets of max(0, bucket_size(b) - 1)
```

i.e. the number of entries currently sharing a bucket with at least one other entry, computed in
`snapshot()`.

Consequently `colRatio()` (`collisions / lookups`) means "current chain overhead per lookup
performed so far", not "chain steps per lookup".

**Why not an insert-time probe?** The obvious cumulative alternative is to call `map.bucket(key)`
before each `emplace` and record whether the bucket was already occupied. That costs a second hash
of the `VKey`/`MKey` on every insert — and a `VKey` hash is two `Dw::hash()` calls, each over five
arbitrary-precision integers. `makeVEdge`/`makeMEdge` are the hottest path in the package; this is
the one place the module is not allowed to slow down. The state-based definition moves that cost to
report time, where it is free.

**Known blind spot:** `garbageCollect()` `.clear()`s both compute caches, so a report-time walk of
a just-cleared cache shows 0 collisions however badly it thrashed beforehand. If a cumulative count
is ever wanted for the compute tables specifically, the insert-time probe drops into the same
`if constexpr` gate with no other design change — the cost argument above applies to the unique
tables, not to the caches.

### D4. Gating: `if constexpr` inside inline methods, not `#ifdef` at call sites

CMake option `EXACT_DD_STATISTICS` (default ON) sets the compile definition `DD_EXACT_STATISTICS`,
mirroring how `EXACT_DD_WITH_GMP` sets `DD_EXACT_WITH_GMP`. `StatisticsConfig.hpp` turns that into
`inline constexpr bool kStatisticsEnabled`.

Every mutating method (`trackLookup`, `trackHit`, `trackInsert`, `snapshot`, `trackActiveEntry`,
`untrackActiveEntry`, `trackGcSweep`, and the five `MemoryManagerStatistics` trackers) is defined
inline in its header with its body wrapped in `if constexpr (kStatisticsEnabled)`. When the option
is OFF the bodies are discarded and the calls collapse to nothing.

Rejected alternatives:

- **`#ifdef` at every call site.** There are over twenty. It would make `makeVEdge`,
  `multiplyRec`, `incRef`/`decRef` and `garbageCollect` substantially harder to read for no gain.
- **Conditionally shaping the structs** (fields present only when enabled). Then `json()` and every
  derived getter would need conditional compilation too. Keeping the fields unconditional means one
  struct layout, one set of definitions, and code that compiles identically either way.

`trackInsert` and `snapshot` take the map **by reference** rather than a precomputed size, so that
the call site evaluates nothing when the gate is off — only a reference is formed.

**ODR hazard.** Because the tracking methods are inline, the definition must agree between the
library and every consumer, or those inline functions have differing definitions across translation
units — an ODR violation no compiler will diagnose. That is why `target_compile_definitions` uses
`PUBLIC`, not `PRIVATE`. It is also an additional reason `EXACT_DD_INSTALL` stays off: an
installed consumer would need the definition recorded in the export.

### D5. A disabled build reports a marker, not zeros

With `EXACT_DD_STATISTICS=OFF`, `PackageStatistics::json()` emits **only**

```json
{ "statistics_disabled": true, "note": "…reconfigure with -DEXACT_DD_STATISTICS=ON." }
```

and deliberately not the usual key tree filled with zeros. A consumer reading
`j["vector"]["unique_table"]["hits"]` then fails structurally instead of silently believing a zero.
`vNodeCount()`/`mNodeCount()` remain accurate in both configurations.

### D6. Templated statistics are defined in a `.cpp` with explicit instantiation

`MemoryManagerStatistics<T>`'s non-inline members live in `MemoryManagerStatistics.cpp`, which ends
with explicit instantiations for `DwVNode` and `DwMNode`. This is the same trick MQT Core uses, and
its purpose here is to keep `<nlohmann/json.hpp>` — ~25k lines — out of `MemoryManager.hpp`'s, and
therefore `DwPackage.hpp`'s, include chain. The stats headers include only
`<nlohmann/json_fwd.hpp>`.

The trade-off: a future `MemoryManager<SomeNewNode>` will not link until its instantiation is added
to that list. The note is at the instantiation site.

### D7. Garbage collection does not reset the compute-table counters

MQT Core's `ComputeTable::clear()` calls `stats.reset()`. Here it would be a no-op anyway —
`reset()` only zeroes `numEntries`, which is snapshotted from the map regardless (D2) — and keeping
`lookups`/`hits`/`inserts` cumulative over the package's whole lifetime is strictly more
informative than resetting them on every collection.

**Consequence when reading a report:** a compute table's `num_entries` reflects the cache *since
the last collection*, while its `lookups`/`hits`/`inserts` are lifetime totals. In the sample in §6
this shows up as `matrix_vector_mult` having `num_entries: 0` alongside `inserts: 99`.

## 4. Instrumentation points

All in `src/dd/exact/DwPackage.cpp` unless noted.

| Site | Calls |
|---|---|
| `DwPackage` ctor | seeds the four `entrySize` fields (`numBuckets` deliberately not seeded — the maps rehash) |
| `makeVEdge` | `trackLookup` before `find`, `trackHit` on the found branch, `trackInsert` after `emplace` |
| `makeMEdge` | same, against `mUniqueStats_` |
| `multiplyRec(mEdge, vEdge)` | same three, against `mvCacheStats_` |
| `multiplyRec(mEdge, mEdge)` | same three, against `mmCacheStats_` |
| `incRef` (both overloads) | `trackActiveEntry()` on the 0→1 transition |
| `decRef` (both overloads) | `untrackActiveEntry()` on the 1→0 transition |
| `garbageCollect` | `trackGcSweep(map)` after each sweep loop, *after* the early exit |
| `MemoryManager` ctor (hpp) | `trackInitialAllocation` |
| `MemoryManager::get` (hpp) | `trackReusedEntries` on the free-list branch, `trackUsedEntries` on the chunk branch |
| `MemoryManager::returnEntry` (hpp) | `trackReturnedEntry` |
| `MemoryManager::allocateNewChunk` (hpp) | `trackAllocation` |

Two placement details that matter:

- `trackGcSweep` sits **after** `garbageCollect`'s early exit (`!force && size < gcLimit`), so
  `gcRuns` counts sweeps that actually ran rather than every call. It also re-derives
  `numActiveEntries = numEntries`, valid because every survivor of a sweep has `ref > 0` by
  construction — this self-heals any drift in the incRef/decRef tracking.
- A node whose refcount saturates at `kMaxRefCount` is never decremented again, so it stays counted
  as active forever. MQT Core has the same asymmetry. `trackGcSweep`'s re-derivation bounds how far
  the number can drift between sweeps.

## 5. Deviations from MQT Core beyond D1–D7

Sections MQT Core reports and this does not, with the reason:

| Absent section | Why |
|---|---|
| `density_matrix` | exact-dd has no `dNode` analogue |
| `real_numbers` | `Dw` edge weights are neither hash-consed nor pooled — they are heavyweight arbitrary-precision values, not simple doubles — so there is no `RealNumber` unique table or memory manager to report |
| `compute_tables.vector_add` / `matrix_add` | `DwPackage::addRec` has no compute table at all. An "add cache hit ratio" cannot exist until an add cache does (and one was measured flat previously — see the perf notes) |
| `getDataStructureStatistics()` | static `sizeof`/`alignof` reporting; low value here, since the numbers understate reality (§7) |
| per-variable breakdown | MQT keeps one `UniqueTableStatistics` per qubit and reports `{total, 0, 1, …}`. `DwPackage` uses one map per node kind, so this would mean per-var counter vectors alongside a single map, leaving `numBuckets`/`loadFactor` global and the breakdown half-meaningful |

Not instrumented, deliberately: the per-call memo maps inside `kroneckerRec`, `outerProductRec` and
`innerProductRec`. They are function-local, not package state, and `innerProduct` is `const`.
Instrumenting them would mean either threading a stats reference through the recursion or promoting
them to members.

Two small behavioural deviations, both guarding against a degenerate case MQT does not hit:

- `MemoryManagerStatistics::trackReturnedEntry()` guards `numUsed > 0` before decrementing, so the
  counter cannot wrap if tracking was skipped for the matching handout.
- `getUsageRatio()` guards `numAllocated == 0` and returns `0.`; MQT's unguarded version would
  report `NaN` with statistics disabled.

## 6. Report schema

```
{ "vector":  { "unique_table":   "unused" | { …11 table keys…, num_active_entries,
                                              peak_num_active_entries, gc_runs },
               "memory_manager": "unused" | { …12 manager keys… } },
  "matrix":  { "unique_table": …, "memory_manager": … },
  "compute_tables": { "matrix_vector_mult": "unused" | { …11 table keys… },
                      "matrix_matrix_mult": "unused" | { … } },
  "active_memory_mib": <double>,
  "peak_memory_mib":   <double> }
```

The eleven table keys are `num_buckets, memory_MiB, num_entries, peak_num_entries, collisions,
hits, lookups, inserts, hit_ratio, col_ratio, load_factor`. The twelve manager keys are
`memory_allocated_MiB, memory_used_MiB, memory_used_MiB_peak, num_allocated, num_allocations,
num_available_for_reuse, num_available_for_reuse_peak, num_available_from_chunks,
num_available_total, num_used, num_used_peak, usage_ratio`. Both sets are MQT Core's, verbatim.
`matrix_vector_mult`/`matrix_matrix_mult` are MQT Core's names for the same two caches.

Sample, from `ExactDDBench`'s GHZ + outerProduct + forced-gc workload at 12 qubits (abridged):

```json
"vector": {
  "unique_table":   { "lookups": 110, "hits": 21, "inserts": 89, "collisions": 1,
                      "num_entries": 23, "peak_num_entries": 89,
                      "num_active_entries": 23, "peak_num_active_entries": 34, "gc_runs": 1 },
  "memory_manager": { "num_allocations": 1, "num_allocated": 2048, "num_used": 23,
                      "num_used_peak": 89, "num_available_for_reuse": 66 } },
"compute_tables": {
  "matrix_vector_mult": { "lookups": 120, "hits": 21, "inserts": 99, "num_entries": 0 },
  "matrix_matrix_mult": "unused" }
```

Reading it: 89 nodes were inserted, 23 survived the forced collection, and 89 − 23 = 66 went back
on the free list. The mv cache shows `num_entries: 0` because the collection cleared it, while its
lifetime counters survive (D7). `matrix_matrix_mult` is `"unused"` because the workload never
multiplies two matrices.

## 7. Reading caveats

**Every MiB figure is a lower bound**, in a way MQT Core's are not. Each node holds one `Dw` per
outgoing edge, and a `Dw`'s five arbitrary-precision integers own heap storage that `sizeof` cannot
see. `memory_MiB`, `memory_used_MiB`, `active_memory_mib` and `peak_memory_mib` all understate
reality, potentially by a large factor at Grover scale where coefficients grow. **Do not read
`peak_memory_mib` as a memory budget.**

`getMemoryMiB()` keeps MQT's `numBuckets × entrySize` formula. For an `unordered_map` this
approximates node storage — `bucket_count() >= size()` at the default `max_load_factor` of 1.0, so
it is within roughly a factor of two — and ignores the bucket-pointer array entirely.

**`collisions` / `col_ratio` are not comparable to MQT Core's** — see D3. The identical key names
are a trap for anyone diffing the two backends' reports.

## 8. Measured overhead

Grover, `ExactDDBench <n>`, GMP backend, `EXACT_DD_STATISTICS` ON vs OFF:

| n | ON (ms) | OFF (ms) | reps |
|---|---|---|---|
| 14 | 755, 756, 762 | 754, 759, 777 | 3 |
| 16 | 2913, 2936 | 2917, 2925 | 2 |
| 18 | 15112, 15430 | 15257, 15426 | 2 |
| 20 | 105871, 106462 | 105794, 106521 | 2 |

Sorted, the runs interleave at every size — neither build is consistently on one side — and the
within-build spread (~2% at n=18, ~0.6% at n=20) swamps the between-build difference. At n=20 the
gap between the means is 9 ms out of 106 s, or 0.008%. The n=18/20 runs were interleaved
(ON, OFF, ON, OFF) so thermal drift over the ~10-minute sweep could not bias one build.

Node counts are byte-identical across builds (3409/54 at n=18, 2468/60 at n=20), confirming the
instrumentation is observation-only.

This is consistent with the recorded profile of this library, where malloc/free is ~45% and `Dw`
arbitrary-precision arithmetic dominates: an integer increment per table probe does not register.
**Do not turn `EXACT_DD_STATISTICS` off "for speed."** The only reason to turn it off is isolating
a `DwPackage`/`Dw` change from any instrumentation at all.

## 9. Tests

`test/test_statistics.cpp`, 29 cases. All but one live in a fixture that `GTEST_SKIP()`s when the
module is compiled out; the exception asserts the disabled-build marker (D5) and is skipped when it
is compiled in. So both configurations run a meaningful suite: 156/156 pass either way.

Coverage groups: unique tables (unused sentinel, entry-size seeding, `numEntries == vNodeCount()`,
hits-without-inserts on rebuild, load factor, the collision bound of D3), reference counting
(active-entry tracking, peak retention, gc re-derivation, gc-run counting past the early exit),
memory managers (seeding, `numUsed == numEntries`, second-chunk allocation, reuse after collection,
the available/used identities), compute tables (hits on repeat, unused sentinel, counter survival
across the gc clear per D7), and the report (section presence, JSON round-trip, stream operator,
memory ordering, `reset()` semantics, snapshot idempotence).

Two traps worth knowing when extending the suite:

- A `|+>` tensor product **creates no nodes at all**. `makeVEdge` sees identical children at every
  level and collapses the whole thing to a single terminal edge, so the unique table is never
  probed. Use a computational-basis state or GHZ when a test needs nodes.
- `ghz()` returns an incRef'd handle. Discarding it keeps the final root alive; what a forced
  collection reclaims is the gate DDs and superseded intermediate states.

## 10. Open items

- **Dangling comment reference.** `TableStatistics.hpp`'s `snapshot()` says "see the
  collision-semantics note on this class", but that note now lives here (D3) rather than in the
  header. Either repoint it at this document or restore a one-line version in the header.
- **Compute-table collision blind spot** (D3) — unmeasured by construction after a collection.
- **`addRec` has no compute table**, so there is nothing to report for addition. Adding one was
  measured flat previously; if that is ever revisited, the report gains two sections for free.
- **Exact-arithmetic metrics were considered and deliberately excluded** from this iteration: peak
  and mean bit-length of `Dw` coefficients, and counts of the expensive `Dw::inverse()`/`gcd()`
  normalization calls. These have no MQT analogue and would partly replace the missing
  `real_numbers` section — arguably the numbers that would best explain exact-backend cost. Scope
  for this pass was MQT parity only.
