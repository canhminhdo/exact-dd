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
  Note that §11 shows the caches genuinely are cleared (262 times at n=20), so this is a live gap,
  not a hypothetical one.
- **`addRec` has no compute table**, so there is nothing to report for addition. Adding one was
  measured flat previously; if that is ever revisited, the report gains two sections for free.
- **Exact-arithmetic metrics** — peak and mean bit-length of `Dw` coefficients, and counts of the
  expensive `Dw::inverse()`/`gcd()` normalization calls. Excluded from this pass, whose scope was
  MQT parity. §11 promotes these from "arguably useful" to **the clear next step**: they are the
  only instrumentation that would see the cost the measurements actually point at.

## 11. Findings

Measured 2026-08-07/08 via `ExactDDBench 14 16 18 20` (GMP backend, `EXACT_DD_STATISTICS=ON`).

### Where Grover's cost goes: not in DD bookkeeping

| n | time (ms) | vUnique peak | vUnique lookups | hit % | gc_runs | chunks | mv lookups | mv hit % | mv peak | peak MiB |
|---|---|---|---|---|---|---|---|---|---|---|
| 14 | 723 | 4107 | 99,259 | 34.7 | 15 | 2 | 134,112 | 48.2 | 4386 | 1.242 |
| 16 | 2,767 | 4110 | 264,039 | 35.7 | 41 | 2 | 352,163 | 48.1 | 4423 | 1.246 |
| 18 | 14,933 | 4111 | 679,722 | 36.6 | 105 | 2 | 896,426 | 48.0 | 4438 | 1.249 |
| 20 | 106,533 | 4114 | 1,712,596 | 37.6 | 262 | 2 | 2,234,961 | 47.9 | 4477 | 1.253 |

**Every structural metric is flat in n while runtime grows 147×.** Peak vector nodes 4107 → 4114
(+0.2%), peak mv cache 4386 → 4477 (+2%), peak memory 1.242 → 1.253 MiB (+0.9%), unique-table hit
ratio 34.7% → 37.6% (improving), mv hit ratio 48.2% → 47.9% (flat), chunk allocations 2 at every
size.

Operation *count* grows 17.3×; runtime grows 147×. The residual is cost per operation:

| n | 14 | 16 | 18 | 20 |
|---|---|---|---|---|
| µs per unique-table lookup | 7.3 | 10.5 | 22.0 | 62.2 |

**8.5× more expensive per probe at n=20 than n=14, on a table of identical size.** That is `Dw`
coefficient growth. This module confirms it by elimination, and the derivation needs no
instrumentation beyond what already exists.

Per-hypothesis verdict:

- **Chunk over-allocation — ruled out.** 2 allocations, 6144 entries, 4114 peak, across a 106 s run.
- **Unique-table thrash — ruled out as a scaling driver.** Load factor 0.38–0.57, a few hundred
  collisions on ~4000 entries, hit ratio constant in n. The 36% hit ratio is inherent to Grover
  (the state changes every gate), not pathology.
- **GC churn — real but cheap.** 262 sweeps of a ~4000-entry map ≈ 1M node visits, against 1.71M
  bignum-hashing lookups. Note this refutes an earlier guess that `gc_runs` might be 0: the table
  peaks at ~4110, just over `kInitialGcLimit` of 4096, at *every* size.
- **Compute-table miss rate — the one open lever, and a constant factor rather than a scaling one.**
  1.16M inserts / 262 clears ≈ 4400 per clear, matching the 4477 peak: the cache fills, is wiped,
  refills. At 48% flat in n it cannot explain 147×.

`matrix_matrix_mult` is `"unused"` at every size — Grover never multiplies two matrices.

### Two null results — do not retry either

**`vUnique_`/`mUnique_` `reserve(4096)` + `max_load_factor(0.75f)` in the constructor.** Cuts
collisions ~22% consistently (885 → 691, 547 → 455, 740 → 616, 425 → 334) and changes no timing:
+0.9%, +2.0%, +1.1%, −0.2%, all inside the ~1–3% noise band with no consistent direction. Reverted.
Two incidental findings: `mUnique_.reserve(4096)` is 68× oversized for a table peaking at 42–60
entries, and because `reserve(n)` is `rehash(ceil(n / max_load_factor()))` while the load factor is
set *afterwards*, the table rehashed anyway and ended at 8192 buckets — larger than the 6421 it
would have reached on its own.

**Replacing `find()`-then-`emplace()` with `try_emplace()` in `makeVEdge`/`makeMEdge`.** The pair
hashes the key twice on every miss; at n=20 that is 1,068,237 redundant `VKey` hashes out of 2.78M
(38.4%), roughly 10.7M `hash_value()` calls over arbitrary-precision integers. Removing all of them
changed nothing: 702 → 710, 2769 → 2811, 14853 → 14843, 104641 → 103715 ms, no consistent
direction. Reverted.

*Re-tested later under the paired design*, because the first pass used best-of-2 against a noise
floor that was underestimated by half — the same error that had wrongly condemned 1a/1b. The
re-test confirms the null: **5 of 12 paired wins** (median +0.5% / +1.4% / −0.3% at n=16/18/20),
i.e. indistinguishable from the 6/12 expected by chance, with a decision rule fixed before the data
was seen. The performance question is closed — do not re-measure it.

**`try_emplace` is nonetheless in the code**, in both `makeVEdge` and `makeMEdge`, kept on the same
grounds as the `Dw::hash` early-out: it does strictly less work (38.4% of this table's key hashes at
n=20 were recomputing a hash just computed), and being performance-neutral it costs nothing to
carry. It requires the placeholder-erase-on-throw guard, because inserting a `nullptr` placeholder
before allocating inverts the original exception ordering, and `p == nullptr` means *terminal* in
this representation — a leftover null entry would read back as a valid terminal edge rather than
failing cleanly. Behaviour-preservation is verified by the value-determined counters matching their
session-long recorded values exactly.

Note also that the first pass's derived bound was wrong. "Unique-table key hashing is at most ~2% of
runtime" assumed the experiment could resolve 2%; it could only exclude an effect above the noise,
~4.4 s of 104.6 s, giving `0.384 × H ≤ 4.4 s` → **H ≤ ~11%**. The correct reading is that hashing
is somewhere under ~11% and that removing 38% of it is not worth a measurable amount — not that
hashing is 2%. The 62 µs per lookup is dominated by the surrounding arithmetic either way. With the default `Inverse` normalization every one of those 1.07M node creations calls
`Dw::inverse()`, and every `Dw` result runs `canonicalize()`, whose `reduceRationalDenominator` does
a gcd across five bignums. That matches the recorded profile (malloc/free ~45%, gcd+division ~16%,
multiplication ~11%) far better than hashing does.

**Next step is a profile of `Dw::inverse`, `canonicalize`, and GMP temporary allocation — not
another micro-optimization of the table path.**

### Methodological note: `collisions` is not run-to-run stable

`VKeyHash` mixes in child node *pointers*, so ASLR changes the bucket distribution between runs. On
byte-identical code at n=14, `collisions` measured 847 / 857 / 879 / 887 (~5% spread) while
`num_buckets` 6421, `lookups` 99,259, `hits` 34,486, `inserts` 64,773 and `peak_num_entries` 4107
were identical every time. Treat `collisions` as a trend across several runs, never as an exact
check; everything value-determined *is* exact, and that set is what validates a refactor.

## 12. Profile: `Dw::inverse`, `canonicalize`, GMP allocation

Measured with `/usr/bin/sample` at 1 ms on a `Release -g` build (timing-parity checked against the
recorded baseline), all runs in one exclusive window with no concurrent load. Plus exact GMP
allocation counts from a scratch harness installing `mp_set_memory_functions`.

### The cost profile *inverts* between n=14 and n=20

| category | n=14 (5,821 samples) | n=20 (50,941 samples) |
|---|---|---|
| libgmp arithmetic | 37.8% | **94.2%** |
| malloc/free | **54.0%** | 4.9% |
| exact-dd's own code | 6.4% | 0.8% |

Top self-time symbols at n=20 are all bignum inner loops: `__gmpn_addmul_1` 20.5%,
`__gmpn_submul_1` 16.2%, `__gmpn_mul_1` 10.8%, `div2` 9.7%, `__gmpn_hgcd2` 9.7%, `__gmpn_sub_n`
5.4%, `__gmpn_add_n` 4.9%, `__gmpn_toom22_mul` 2.7% — multiplication, division and GCD on wide
operands.

**This corrects the previously recorded profile** (malloc/free ~45%, gcd+division ~16%,
multiplication ~11%). Those figures match the n=14 regime almost exactly and describe n=20 not at
all. The old profile was taken at a size where the answer is qualitatively different.

### The growing term is `Dw::inverse()`

Self-time attributed by call-tree ancestor. `Dw::inverse` is called across a TU boundary so it
survives `-O3` inlining; `canonicalize` does not, hence the separate build below.

| ancestor | n=14 | n=20 |
|---|---|---|
| under `Dw::inverse` | 18.5% | **45.4%** |
| under `multiplyRec`/`addRec` (excl. inverse) | 41.1% | 31.8% |
| under `makeVEdge`/`makeMEdge` (excl. inverse) | 36.4% | 22.5% |

`Dw::inverse()` goes from under a fifth of runtime to nearly half. It is called from
`makeVEdge`/`makeMEdge` under the default `NormalizationStrategy::Inverse`, once per node creation
whenever `eta != 1`.

From a `-fno-inline-functions` build at n=14 (**attribution only — 1.8× slower, timings not
comparable to anything**): `canonicalize` not reached via `inverse` accounts for 24.8%, `inverse`
15.4%, other `Dw::` frames 35.2%, outside `Dw::` 24.6%. So `canonicalize` is a real cost in its own
right, not merely a passenger inside `inverse`.

### GMP allocation: constant count per operation, growing width

| | n=14 | n=20 |
|---|---|---|
| allocations | 23,786,085 | 429,891,917 |
| **allocs per unique-table lookup** | **239.6** | **251.0** |
| bytes allocated (cumulative) | 445 MiB | 53.3 GiB |
| bytes per lookup | 4,702 | 32,649 |
| peak live | 2.77 MiB | 21.96 MiB |
| allocations of 8 bytes | 80.8% | 78.7% |
| allocations of 16 bytes | 9.5% | 11.5% |

Two things fall out. **Roughly 250 GMP temporaries per DD operation, essentially constant in n** —
a fixed algorithmic cost of how `Dw` arithmetic is written, not something that grows. And ~90% of
them are one or two limbs: GMP heap-allocates even a single 64-bit value, and cumulative traffic
(53 GiB) exceeds peak live (22 MiB) by ~2,400×, so it is pure churn.

Meanwhile *bytes* per operation grow 7× (4.7 KB → 32.6 KB). That is the whole story in one line:
the number of temporaries is fixed, their width grows, so malloc's share collapses while bignum
arithmetic's share explodes — and per-operation cost rises 8.5× exactly as §11 measured.

### Line-level: the cost is one gcd statement

Attributing libgmp self-time under `Dw::inverse` to its nearest source line at n=20:

| line | code | share of `inverse()` |
|---|---|---|
| **Dw.cpp:188–189** | the 4-way `gcd(gcd(gcd(a,b),gcd(c,d)), e)` | **64.6%** |
| Dw.cpp:191 | `a_ /= g` | 9.3% |
| Dw.cpp:195 | `e_ /= g` | 9.1% |
| Dw.cpp:247/257/267 | the two `mulTuple` convolutions | 10.2% combined |

**That single gcd statement is 28.8% of total n=20 runtime.** It lives in
`reduceRationalDenominator`, i.e. `canonicalize()` running on `inverse()`'s result — which
necessarily carries a large denominator, since `d = x² − 2y²` is quartic in the coefficients.

This also resolves an ambiguity flagged above: `__gmpn_addmul_1`/`submul_1` topping the n=20 symbol
list did **not** mean tuple multiplication dominated. Those primitives are the inner loops of GCD
and division too, and the line attribution shows that is where they were being called from.

### Counter probe on the gcd (temporary scaffolding, since removed)

| | n=14 | n=20 |
|---|---|---|
| phase-2 gcd calls | 111,011 | 1,712,596 (= exactly the unique-table lookup count) |
| result was 1 (nothing to reduce) | 18.3% | 16.4% |
| gcd ops: current tree → early-exit fold | 444,044 → 383,069 (1.16×) | 6,850,384 → 6,005,902 (1.14×) |
| **operand width, `e`** | **1,348 bits** | **16,400 bits** |
| operand width, max coefficient | 1,355 bits | 16,419 bits |
| phase-1 `while (e_ % 2 == 0)` iterations per call | 0.279 | 0.341 |

**Operand width grows 12× from n=14 to n=20 — that is the per-operation cost growth, located.**
The gcd count per DD operation is fixed at one; only the numbers get wider.

Two consequences for the candidate list:

- **Phase 1's division loop is not worth optimising.** It averages ~0.3 iterations per call, so
  replacing it with a shift can save essentially nothing. Ruled out for free.
- **An early-exit fold is not primarily an early-exit win.** The gcd is 1 only ~16% of the time, so
  the op count barely moves (1.14×). Any benefit would have to come from *operand sizes*: the
  current balanced tree performs two large×large gcds (`gcd(a,b)` and `gcd(c,d)`), whereas folding
  progressively from `e` performs one, after which the running gcd is the common factor and is
  typically small. That is a hypothesis about cost, not a measured fact.

### First measured win: skip `eta * etaInv` in normalization (2a), plus 2c

`makeVEdge`/`makeMEdge` took `eta` from the leftmost nonzero child weight and then multiplied
*every* child by `etaInv` — including the one `eta` came from, whose product is exactly one by
construction. That product ran a full Z[w] convolution and then canonicalized a fraction whose
denominator is `inverse()`'s quartic norm, i.e. one of the large gcds above, spent to rediscover the
value 1. Assigning `Dw::one()` instead (**2a**), together with replacing `reduceSqrt2Power`'s
subtraction-based parity guard with a low-bit comparison and its halvings with shifts (**2c**):

| n | baseline | 2a + 2c | change |
|---|---|---|---|
| 14 | 702 ms | 661 ms | −5.9% |
| 16 | 2,769 ms | 2,579 ms | −6.9% |
| 18 | 14,853 ms | 13,865 ms | −6.6% |
| 20 | 104,641 ms | 100,913 ms | −3.6% |

All nine value-determined counters byte-identical at every size, as required — 2a changes how a
weight is computed, not what it equals.

Four further changes landed on top, all provable work reductions inside `Dw`:

- **1a** — the norm tuple `z * conj(z)` is always `(x, y, 0, −y)`, so `normXY()` computes
  `x = a²+b²+c²+d²` and `y = ab+bc+cd−ad` directly: 8 products instead of `mulTuple`'s 16, with two
  components known structurally.
- **1b** — `mulByConjNormNumer()` specialises the multiply against `{x, −y, 0, y}`, whose zero
  component and repeated ±y collapse 16 products to 8.
- **2b** — `reduceRationalDenominator` phase 1 extracts powers of two with `lsb()` and one shift
  instead of a loop of full-width bignum divisions.
- `makeMEdge`'s normalization loop starts at the leftmost nonzero index (earlier children have zero
  weight by construction, and zero × etaInv is zero).

Measured on their own, on top of 2a+2c: **~1%** — median −1.3% / −1.0% / −0.6% at n=16/18/20, with
**11 of 12 paired wins (sign test p ≈ 0.003)**.

**Cumulative, original HEAD versus all of the above** (paired, alternating prebuilt binaries, BASE
running first in each pair so drift penalises FINAL — these are lower bounds):

| n | baseline median | final median | change | paired wins |
|---|---|---|---|---|
| 16 | 2,810 ms | 2,515 ms | **−10.5%** | 5/5 |
| 18 | 14,970 ms | 13,794 ms | **−7.9%** | 5/5 |
| 20 | 106,611 ms | 99,136 ms | **−7.0%** | 2/2 |

Twelve of twelve paired wins. Note the parts do not sum to the whole and should not be quoted as if
they did: 2a removes a gcd that 1a/1b would otherwise have fed, so the changes interact.

### Remaining redundant hashing, measured and closed

A counter probe (temporary, `-DDD_EXACT_HASH_PROBE`, since removed) answered what is left:

| | Grover n=18 | default suite |
|---|---|---|
| `Dw::hash` calls | 2,417,436 | 425,942 |
| — exactly `Dw::one()` | 63.1% | 62.5% |
| — exactly `Dw::zero()` | 4.3% | 5.4% |
| — **trivial total** | **67.4%** | **67.9%** |
| non-trivial calls, mean max-coefficient limbs | 26.0 (max 101) | 5.2 (max 19) |
| redundant key hashes from `find`-then-`emplace` | 35.8% of all key hashes | 36.0% |

**Two-thirds of every `Dw::hash()` call hashes a compile-time constant**, which confirms the
structural argument from `makeVEdge`'s postconditions: whichever branch runs, one child weight ends
as exactly `Dw::one()`, and when `children[0]` is zero both weights are constants.

**But counts are the wrong metric here, and the probe is what shows it.** Weighting each call by the
limbs it actually hashes (5 coefficients × mean limb count), the trivial calls are only **7.4% of
the limb-hashing work at Grover n=18** — 1.63M cheap single-limb calls against 788K calls averaging
26 limbs each. Since all key hashing is bounded at ≤11% of runtime, the ceiling for eliminating
*all* trivial hashing is **≤0.8%**, well under the 4.2% noise floor. The share rises to 28.9% of the
work on the small-circuit default suite, so it matters more for tiny workloads than at 20 qubits.

An `isZero()`/`isOne()` early-out **was implemented anyway** — it is provably less work and three
lines — and measured flat exactly as that ceiling predicts: **8/17 paired wins, sign test p = 0.69**
(median-of-ratios +0.9% / −0.2% / +0.7% at n=14/16/18). It is kept on the grounds that it does
strictly less work, not on measured performance. Do not expect it to show up in a benchmark.

Its correctness rests on `isOne()` testing all six fields (exact) and `isZero()` ignoring `k_`/`e_`,
which can only produce a legal hash collision between unequal "zero-ish" values — equal values still
hash equal, the only requirement — and `canonicalize()` forces `k_ = 0, e_ = 1` when all
coefficients vanish, so that case never arises. The hash-consing test is the behavioural check: an
equality-inconsistent hash would stop structurally identical DDs from sharing nodes.

This also corrected a badly-chosen criterion in the decision rule that gated the probe: it required
"non-trivial calls averaging >2 limbs", on the intuition that expensive non-trivial calls made the
trivial ones worth skipping. That is backwards — the *larger* the non-trivial calls, the *smaller*
the trivial share of total work, and the less an early-out can buy.

**The per-call memos are negligible.** `outerProductRec` and `innerProductRec` never appear in
Grover at all, and on the default suite (which includes the GHZ measure sweep, i.e. QRAT's actual
`measureOneQubit` path) they account for ~1,300 of 439,648 key hashes — 0.3%. The `try_emplace`
reentrancy caveat that applies to them is therefore moot.

### Measurement method: the single-run spread is 4.2–4.3%, not 1–3%

Attribution needed a paired design, because single-run spread here is 4.2–4.3%: five reps at each
size, the two prebuilt binaries run alternately so drift cannot favour either. **2a carries most of
it** (consistently faster at all four sizes on its own). **2c is nil at n=16 (3/5 paired wins,
+0.3% median) but ~2% at n=18 (5/5 paired wins, −2.3% median, sign-test p ≈ 0.03)** — consistent
with its mechanism, since the two bignum subtractions it removes cost O(operand width) and width
grows 12× from n=14 to n=20.

Also folded in: `makeMEdge`'s normalization loop now starts at the leftmost nonzero index rather
than 0, since every earlier child has zero weight by construction and zero × etaInv is zero.

### What this licenses, and what it does not

The two suspects named at the end of §11 are confirmed, with `Dw::inverse` the dominant and growing
one. Two concrete levers follow, both **unmeasured hypotheses** that must clear the same bar the
last two changes failed:

- **Avoid the field inverse where an exact division would do.** Under `Gcd` normalization the
  divisor divides every weight exactly in D[ω], so a field inverse is not required — yet
  `gcdNormalize` in `DwPackage.cpp` still calls `eta.inverse()`. Switching strategy as the code
  stands therefore does *not* dodge the cost.
- **Reduce temporaries.** ~250 per DD operation is a lot for the arithmetic being done; at n=14 that
  alone is over half the runtime.

What the profile does *not* say is that either lever will pay. It says where the time is, nothing
more.

### The weights were copied three times per insert (fixed), and the lazy hash is dead

Two questions were raised together: should `Dw::zero()`/`Dw::one()` (and the `DwEdge` equivalents)
be cached, and can the repeated `edge.w.hash()` in `VKeyHash`/`MKeyHash` be avoided? Investigating
both turned up a third thing, larger than either.

**First, a premise worth recording so it is not rediscovered:** the compute tables never hash a `Dw`.
`mvCache_`/`mmCache_` key on `std::pair<DwMNode *, DwVNode *>` via `PtrPairHash`, pointers only —
which is why `multiplyRec` re-applies `scale` to the cached result. *Every* `Dw` hash in the package
comes from a `makeVEdge`/`makeMEdge` unique-table probe.

**The fix that landed.** `makeVEdge`/`makeMEdge` copied every child weight three times per insert
(into `key`; again when `key` was copied into the map; again into the node) where the design needs
one. Copies 2 and 3 are removable with `std::move`, plus `eta` into the returned edge and
`children[0]` on the redundant-node-elimination path. Measured with a copy-counting probe on the
default suite:

| | `Dw` copies | `Dw` moves |
|---|---|---|
| moves reverted | 2,937,926 | 456,685 |
| moves in place | 2,493,849 | 900,762 |
| delta | **−444,077 (−15.1%)** | +444,077 |

Exactly conserved, so the moves demonstrably happen. Timing, paired design, baseline first:
**9/12 paired wins, median −1.21%** — n=16 4/5 (−1.21%), n=18 4/5 (−1.22%), n=20 1/2 (−0.10%),
sign-test p = 0.146. Real but small, and it vanishes at n=20 exactly as the profile predicts:
malloc/free is 54.0% of runtime at n=14 but 4.9% at n=20 (§12).

Note the subset matters. An earlier run of only the two unique-table moves measured **4/12 wins,
median +0.2%, p = 0.75** — a clean null. The win only appears with all of them.

**The lazy hash is dead — measured, not argued.** A probe added `mutable std::size_t probeHash_` to
`Dw`, memoizing on first `hash()`:

- **Warm hits: 0 of 83,549 non-trivial calls**, if warmed the obvious way. The reason is structural:
  the key is built by *copying* `children` before anything is hashed, so the cache lands on the key's
  copy — which libc++ never consults again, since it stores the hash in the node.
- Warmed instead on the weights that survive into the node, **20% of non-trivial arriving weights are
  already warm** (16,790/83,549). Non-trivial calls are 31% of all `hash()` calls (the rest exit via
  the zero/one early-out, measured 69% here vs 67.4% previously).

So the ceiling is 0.20 × 0.31 × (all key hashing, bounded at ~11% by the `try_emplace` experiment)
≈ **under 1%**, against a 4.2% noise floor — and it would need hand-written move operations plus a
`k_` type change to avoid growing every DD node. Not worth it. **Do not revisit.**

Two by-products worth keeping:

- **Soundness was confirmed, not just argued.** 83,549 warm hits, **0 mismatches**. `Dw` is
  immutable after construction: every write to `a_`–`e_`/`k_` is inside
  `canonicalize()`/`reduceSqrt2Power()`/`reduceRationalDenominator()`, and `canonicalize()` has
  exactly one call site, the 6-arg constructor.
- **The move hazard, if anyone does try this.** Boost's `gmp_int` move-assign is an `mpz_swap`, so
  the moved-from object keeps the destination's old limbs. An implicit move would *copy* a cache
  field while the limbs *swap*, leaving the source describing a value it no longer holds. Both move
  operations must be hand-written to swap the cache alongside.

**`eta == 1` on entry**, measured in the same probe: `makeVEdge` 43%, `makeMEdge` 67%. That is the
fraction of calls where children pass through the normalization block unrewritten.

### Deferred, with reasons

- **Caching `Dw::zero()`/`one()`/`DwEdge::zero()`** as `const &` to function-local statics. Helps
  only at *assignment* sites (`children[i].w = Dw::one()`, `c = mEdge::zero()`), where it turns 1–2
  malloc/free pairs into zero; at *construction* sites (`mEdge{p, Dw::one()}`, `return
  vEdge::zero()`) it is a wash, because the caller still copies. Source-compatible at all call sites
  in exact-dd and QRAT (checked). A small-n win against a 20-qubit target.
- **Pointer-keyed unique tables** (MQT Core's design: `unordered_set<DwVNode *>` with a dereferencing
  hash/eq) would delete copy #1 too and drop the entry from ~224 to 24 bytes. But it forces
  allocating the node *before* the lookup and returning it to the pool on the ~35% hit path, undoing
  part of the `try_emplace` work.
- **The gcd** at `reduceRationalDenominator` (64.6% of `inverse()` = **28.8%** of n=20 runtime,
  not the 62% quoted in earlier drafts) is still the only lever with a ceiling above ~6%.
  One detail: it is guarded by `if (e_ != 1)`, so it fires *only* on values carrying a rational
  denominator — i.e. everything downstream of `Inverse` normalization's `etaInv`. The untested
  restructure is a progressive fold seeded with `e_` (usually the smallest operand) that bails the
  moment the running gcd reaches 1, instead of the balanced tree `gcd(gcd(gcd(a,b),gcd(c,d)),e)`
  that always computes all four full-width gcds.

### A comment that promised a guard the code did not have

The `try_emplace` rework left a comment stating that a `bad_alloc` from `vMemory_.get()` erases the
null placeholder before propagating. The comment landed; the `try`/`catch` did not. Since
`p == nullptr` means *terminal* in this representation, a leftover null entry would be found by a
later probe and returned as a valid terminal edge — a silent wrong answer rather than a clean
failure. The guard is now actually present at both sites.

### Copy/move audit: the arithmetic operators were the whole story

A full copy/move audit of the library, prompted by the `std::move` result above. It split cleanly
into work that removes **bignum arithmetic** and work that only removes **allocations** — and only
the first kind measured.

**Tier 1 — `Dw::operator*`, `operator+`, `inverse()`.** `operator*` was making **12 needless
`Integer` copies per call**: eight to materialise the two `std::array<Integer, 4>` temporaries that
`mulTuple` took by reference (every caller has the coefficients as separate members), and four more
copying the dead result tuple into the constructor's by-value parameters. `mulTuple` now takes eight
`const Integer &`; the result is moved. `operator+` scaled one operand by a provably-zero shift
(`k` is the max of the two) and multiplied four coefficients by a denominator of 1; both are now
guarded, with a no-copy fast path. `inverse()` moves its dead numerator and denominator.

**Measured: 12/12 paired wins, median −5.36%, sign test p = 0.0005** — n=16 −7.84%, n=18 −5.04%,
**n=20 −2.72%**. The first result in this effort that survives at n=20, exactly as the profile
predicts for work removed from the 94.2% libgmp bucket rather than the 4.9% malloc one.

**Attribution caveat — this figure is NOT Tier 1 alone.** The baseline binary was built by
`git stash push -- src/dd/exact/Dw.cpp`, which reverts that file to HEAD, so the comparison also
picks up every earlier `Dw.cpp` change: the `inverse()` closed forms (`normXY`,
`mulByConjNormNumer`), the `reduceSqrt2Power` parity test, the `reduceRationalDenominator` shift,
and the `Dw::hash` zero/one early-out. Read −5.36% as *all `Dw.cpp` arithmetic work versus HEAD*.
Isolating Tier 1 needs a baseline carrying the earlier changes and only those; it has not been
measured.

**Attribution: it is essentially all `operator*`.** A counter probe (temporary, since removed) at
Grover n=16: `operator*` runs **1,368,959** times, 81% of which short-circuit on 0/1, leaving
**261,745** non-trivial multiplications — against `operator+`'s **15,380** calls in total. Even if
`operator+` were made free it could not account for a 5% result, so the `operator+` guards are kept
on strictly-less-work grounds, not as the cause.

The same probe corrected an assumption in the plan: `e != 1` is **not** rare. At Grover n=16 the
left operand has a non-unit denominator on **83%** of additions and only **16%** take the fast path
— under `Inverse` normalization, values downstream of `etaInv` dominate. "Everything that is not
downstream of `inverse()`" was the wrong intuition about which case is common.

**Tier 2–4 — allocation-only work: null, as predicted.** Gate matrices are now function-local
statics returned by `const &` (they were rebuilt per gate application, several also redoing a
negation or conjugation); the remaining dead locals are moved into `makeVEdge`/`makeMEdge`;
`wrapWithControl` takes its edge by value; `const` was dropped where it blocked a move.
**8/12 paired wins, median −0.35%, p = 0.39** — consistent at n=16 (5/5, −0.94%) and nothing beyond.
Kept on the same footing as `try_emplace` and the hash early-out.

**Cumulative, all `Dw.cpp` arithmetic work plus Tiers 2–4, versus HEAD's `Dw.cpp`: 12/12 paired
wins, median −5.86%, p = 0.0005** — −8.35% / −5.32% / −3.29% at n=16 / 18 / 20. Same baseline
caveat as above: this is the whole arithmetic effort, not the audit's items alone.

The one cleanly-attributed comparison here is **Tiers 2–4**, whose two binaries differed only in
those changes.

**The new tests were mutation-checked.** Tier 1 rewrites arithmetic rather than shuffling
ownership, so `test_dw.cpp` gained randomised distributivity, group-axiom and complex-evaluation
cross-checks that deliberately exercise both `operator+` paths (asserting each is reached at least
20 times, so neither can silently go untested). Forcing the fast path unconditionally fails 3 tests;
corrupting one `mulTuple` coefficient fails 8.

**Traps found and deliberately not "fixed":**

- `measureOneQubit`'s `return {projected, innerProduct(projected, projected)};` — braced-init is
  sequenced left to right, so moving `projected` inline would make `innerProduct` read a moved-from
  value. The probability is now hoisted into a local first, then both are moved.
- `resultUnit` in `multiplyRec`/`kroneckerRec`/`outerProductRec` is stored in the memo *and* read in
  the return expression, so it cannot be moved.
- `MemoryManager`'s defaulted move is correct: moving `vector<vector<T>>` transfers the outer
  buffer, leaving the inner vectors' data — which `chunkIt_`/`chunkEndIt_`/`available_` point into —
  untouched.

### A canonical-form regression the suite could not see

While summarising the above, `reduceRationalDenominator`'s phase 1 was found reduced from a loop to
a single shift:

```cpp
if (!boost::multiprecision::bit_test(e_, 0)) { e_ = e_ >> 1; k_ += 2; }   // wrong
```

That extracts **one** factor of 2 from `e_` instead of all of them, so `e_` can end up even and the
"e is odd" half of the canonical form is lost. The consequence is not a wrong arithmetic result but
a **split representation**: 1/4 canonicalizes to `k=2, e=2` when written with the 4 in `e`, and to
`k=4, e=1` when written with it in `k`. `operator==` and `hash()` compare the stored fields, not the
value, so those two compare unequal — hash-consing would treat one weight as two.

All 163 tests passed with this in place, because nothing constructed a value needing more than one
factor removed (the `2^1` case is indistinguishable). `Dw.RationalDenominatorIsAlwaysOdd` now
asserts the invariant directly for 2^0…2^6 and checks that both spellings of the same number
compare and hash equal; it fails on 2^2 and up against the broken version.

Phase 1 is restored to the `lsb()` form, which extracts every factor in one shift — correct, and
still cheaper than the original `while (e_ % 2 == 0)` loop of full-width divisions and modulos.

### Closing the two open items: the operator rewrite isolated, and the gcd lever

**Part A — the operator rewrite, measured against a correct baseline.** The earlier −5.36% used a
baseline built by reverting `Dw.cpp` to HEAD, so it contained every change in that file. Rebuilding
the baseline to differ *only* in the four operator-rewrite edits (`mulTuple`'s 8-argument form,
`operator*`'s moves, `operator+`'s guards and fast path, `inverse()`'s moves), with the restored
phase-1 `lsb()` extraction present in both arms:

**12/12 paired wins, median −1.16%, sign test p = 0.0005** — n=16 −2.84%, n=18 −1.09%,
n=20 −0.80%.

So the rewrite is real but small, and — importantly — it **shrinks with n**, because the work it
removes (an `mpz_init_set` per multiplication) is a per-operation constant while the work that
remains grows with coefficient width. The remaining ~4.2% of the −5.36%, and the part that survives
at n=20, belongs to the *earlier* `Dw.cpp` changes: the `inverse()` closed forms, the
`reduceSqrt2Power` parity test and the `Dw::hash` early-out. The narrative "the operator rewrite is
what survives at n=20" was wrong.

**Part B — the phase-2 gcd lever is closed, negatively.** A counter probe (temporary, since removed)
at Grover n=16 / n=18:

| | n=16 | n=18 |
|---|---|---|
| calls with `e != 1` | 223,346 | 553,509 |
| `g == e` (division was exact) | 14% | 14% |
| `g == 1` (nothing to reduce) | 22% | 22% |
| mean bits: `e` | 3,300 | 7,500 |
| mean bits: **`g`** | **2,100** | **4,700** |
| mean bits: `gcd(e, a)` — the fold's value after one step | 2,100 | 4,700 |

Two conclusions, both negative and both worth keeping:

- **`g == e` at 14% is below the 40% bar fixed in advance**, so the exact-division fast path (four
  `mpz_divisible_p` in front of the gcd) and its upstream variant in `makeVEdge` are **ruled out**.
- **The progressive fold's premise is false.** It assumed that after `gcd(e_, a_)` the running value
  would be small, making the remaining three steps cheap small × large gcds. Measured, that value is
  **~64% of `e`'s width** — every step stays large × large. Measuring it anyway confirmed the
  prediction: **5/12 paired wins, median +0.01%, p = 0.77**. The fold is kept (at most four gcds
  with a 22% early exit, versus the tree's unconditional four) but its comment now says plainly that
  it is not a speed-up.

**What this says about the remaining 28.8%.** The reduction is not wasted work to be avoided: `g > 1`
on 78% of calls and averages nearly two thirds of the denominator's width, so the fraction genuinely
needs reducing and the gcd genuinely finds a large factor. The cost is intrinsic to carrying `Q[ω]`
denominators through `Inverse` normalization, not an artifact of how the gcd is arranged. **Anything
further has to attack the denominators' existence, not the gcd** — and the two routes to that
(exact division at 14%, or a different normalization strategy) are respectively too rare and out of
scope by decision.

Two new canonical-form tests guard this area, since a wrong reducing factor does not corrupt values,
it splits one weight into two representations and silently stops hash-consing:
`Dw.RationalDenominatorIsAlwaysOdd` and `Dw.ReducingFactorUsesEveryCoefficient`. Both were
mutation-checked. Note that the obvious mutation (dropping an operand from the fold) is caught by
the *distributivity* and *Grover* tests instead, because it produces wrong values; the mutation that
produces the split-representation mode is skipping the reduction entirely, and that one makes the
algorithm suite hang — coefficients grow without bound — which is itself a measure of how
load-bearing this gcd is.

### `Inverse` normalization is a net loss at scale — 5.7× at n=20

The bignum analysis argued that Grover's amplitudes live entirely in D[ω] (H, X and multi-controlled
Z are exactly representable and D[ω] is closed under + and ×), so `e = 1` is all the mathematics
requires — and that the 16,400-bit denominators at n=20 are manufactured by `Inverse` normalization
dividing sibling weights through a **field** inverse, then removed again by the phase-2 gcd. That
predicted ~37% of runtime spent cancelling a fraction that need not exist. Measuring the three
strategies end to end shows the effect is far larger than 37%, because the fraction does not merely
cost time — **it makes the per-operation cost grow with n.**

`ExactDDSimulation` gained a defaulted strategy parameter and `ExactDDBench` a
`--strategy=none|inverse|gcd` selector. Grover, every run verified at `P(|1...1>) = 1`:

| n | `Inverse` ms | `None` ms | `Inverse` µs/lookup | `None` µs/lookup | none/inverse |
|---|---|---|---|---|---|
| 10 | 55.8 | 88.1 | 3.94 | 2.67 | 1.58× |
| 12 | 172.9 | 272.8 | 4.34 | 2.61 | 1.58× |
| 14 | 589.9 | 808.6 | 5.48 | 2.60 | 1.37× |
| 16 | 2,449 | 2,358 | 8.64 | 2.60 | **0.96×** |
| 18 | 14,096 | 6,462 | 19.49 | 2.63 | **0.46×** |
| **20** | **102,962** | **17,977** | **56.91** | **2.71** | **0.17×** |

**`None`'s per-operation cost is flat across the whole range (2.60–2.71 µs/lookup). `Inverse`'s
grows 14×, from 3.94 to 56.91.** That is the fraction cost scaling with coefficient width, and it is
what the profile had been measuring all along without the comparison to make it legible.

`Inverse` does buy sharing — it needs 2.3× fewer node operations at n=10, rising to 3.7× fewer at
n=20 — but that advantage grows far too slowly to offset a per-operation cost that grows 14×. The
crossover is at **n=16**, and by n=20 `None` is **5.7× faster**. Peak live nodes are ~4,100–4,300
under both (the GC limit binds either way), so the larger DD costs nothing in memory.

The mechanism was confirmed rather than inferred, with a temporary counter (since removed), summed
over n=10/12/14:

| | `Dw::inverse()` calls | phase-2 gcd entries |
|---|---|---|
| `Inverse` | 88,751 | 132,399 |
| **`None`** | **0** | **0** |
| `Gcd` | 3,397,416 | 4,627,835 |

`None` eliminates rational denominators completely, exactly as predicted.

**`Gcd` closes the exact-division project, negatively.** Its wall time is uninformative (it calls
`inverse()` 38× *more* than `Inverse`, since `gcdNormalize` and `Dw::gcd`→`roundedDivide` both go
through the field inverse), but its *size* numbers decide the question. At n=14 it needs 286,815
lookups against `None`'s 311,507 and `Inverse`'s 107,673 — i.e. **`Gcd`'s compaction is within 8% of
no normalization at all**, nowhere near `Inverse`'s. So a properly implemented `Gcd` (exact D[ω]
division, no field inverse) would land at roughly `None`'s operation count *and* `None`'s
per-operation cost — which is simply `None`, with a gcd computation added. **Not worth building.**

**Caveat, and the reason the default was not changed here.** This is one workload. Grover has very
few distinct amplitudes, so `Inverse`'s scalar-factor sharing has unusually little to find while
paying its full cost; a circuit with more amplitude diversity could behave differently. Before
changing the shipped default, repeat this on the GHZ measure sweep (QRAT's actual measurement path)
and on QFT/Shor. The strategy plumbing added here makes that a one-flag experiment.

A `Grover.AllNormalizationStrategiesAgreeOnTheFullSearch` test now guards all three strategies
through a full 5-qubit search, comparing all 32 amplitudes exactly against `Inverse`. It was
mutation-checked (suppressing the propagated `eta` makes it fail). Before it, the only strategy
coverage was `DwGateBuilders.NormalizationStrategyInvariance`, which exercises matrix gate builders
at 3 qubits and never reaches `makeVEdge`, `multiply` or `add`.
