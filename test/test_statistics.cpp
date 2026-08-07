#include "dd/exact/DwGateMatrixDefinitions.hpp"
#include "dd/exact/DwPackage.hpp"
#include "dd/exact/MemoryManager.hpp"
#include "dd/exact/statistics/StatisticsConfig.hpp"

#include <gtest/gtest.h>
#include <nlohmann/json.hpp>

#include <cstddef>
#include <sstream>
#include <vector>

using dd::exact::Dw;
using dd::exact::DwPackage;
using dd::exact::DwVNode;
using dd::exact::kStatisticsEnabled;
using dd::exact::MemoryManager;
using dd::exact::NormalizationStrategy;
namespace gates = dd::exact::gates;

namespace {

/// GHZ = (|0...0> + |1...1>)/sqrt(2), built through applyOperation so that
/// the unique tables, the matrix-vector cache, and the memory managers all
/// see traffic. Leaves the returned state incRef'd.
DwPackage::vEdge ghz(DwPackage &pkg, std::size_t n) {
    auto state = pkg.makeZeroState();
    pkg.incRef(state);
    const std::size_t control = n - 1;
    state = pkg.applyOperation(pkg.makeSingleQubitGateDD(control, gates::h()), state);
    for (std::size_t q = 0; q < control; ++q)
        state = pkg.applyOperation(pkg.makeControlledSingleQubitGateDD(control, q, gates::x()), state);
    return state;
}

/// A state with 2^levels distinct amplitudes, hence roughly that many
/// distinct vector nodes -- enough to exhaust MemoryManager's first chunk.
/// NormalizationStrategy::None keeps this fast: it skips the Dw::inverse()
/// that Inverse would run per node, which is by far the dominant cost.
DwPackage::vEdge manyDistinctNodes(DwPackage &pkg, std::size_t levels) {
    const std::size_t size = std::size_t{1} << levels;
    std::vector<Dw> amplitudes;
    amplitudes.reserve(size);
    for (std::size_t i = 0; i < size; ++i)
        amplitudes.emplace_back(static_cast<long>(i + 1));
    return pkg.makeStateFromVector(amplitudes);
}

/// Every test in the fixture inspects counters, so all of them are
/// meaningless -- and would fail -- in a build with the module compiled out.
class Statistics : public ::testing::Test {
protected:
    void SetUp() override {
        if constexpr (!kStatisticsEnabled) {
            GTEST_SKIP() << "built with EXACT_DD_STATISTICS=OFF";
        }
    }
};

} // namespace

// ---------------------------------------------------------------------
// Unique tables
// ---------------------------------------------------------------------

TEST_F(Statistics, FreshPackageUniqueTablesAreUnused) {
    const DwPackage pkg(3);
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.vUniqueTable.lookups, 0U);
    EXPECT_EQ(stats.mUniqueTable.lookups, 0U);
    EXPECT_EQ(stats.json()["vector"]["unique_table"], "unused");
    EXPECT_EQ(stats.json()["matrix"]["unique_table"], "unused");
}

TEST_F(Statistics, EntrySizeIsSeeded) {
    const DwPackage pkg(3);
    const auto stats = pkg.statistics();
    EXPECT_GT(stats.vUniqueTable.entrySize, 0U);
    EXPECT_GT(stats.mUniqueTable.entrySize, 0U);
    EXPECT_GT(stats.mvComputeTable.entrySize, 0U);
    EXPECT_GT(stats.mmComputeTable.entrySize, 0U);
}

TEST_F(Statistics, SnapshotNumEntriesMatchesNodeCount) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.vUniqueTable.numEntries, pkg.vNodeCount());
    EXPECT_EQ(stats.mUniqueTable.numEntries, pkg.mNodeCount());
    EXPECT_GT(stats.vUniqueTable.numEntries, 0U);
    pkg.decRef(state);
}

TEST_F(Statistics, InsertsEqualEntriesWhenNothingCollected) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    // Nothing has been collected yet, and makeVEdge/makeMEdge never discard
    // a freshly allocated node the way MQT Core's searchTable can, so every
    // insert is still present in the table.
    EXPECT_EQ(stats.vUniqueTable.inserts, stats.vUniqueTable.numEntries);
    EXPECT_EQ(stats.mUniqueTable.inserts, stats.mUniqueTable.numEntries);
    pkg.decRef(state);
}

TEST_F(Statistics, RepeatedIdenticalConstructionRegistersHits) {
    DwPackage pkg(3);
    // A computational-basis state, not a |+> product: the latter has
    // identical children at every level, so makeVEdge collapses it to a
    // single terminal edge and never reaches the unique table at all.
    const std::vector<bool> pattern{true, false, true};
    (void)pkg.makeBasisState(pattern);
    const auto first = pkg.statistics();
    ASSERT_GT(first.vUniqueTable.numEntries, 0U);

    (void)pkg.makeBasisState(pattern);
    const auto second = pkg.statistics();

    EXPECT_GT(second.vUniqueTable.hits, first.vUniqueTable.hits);
    EXPECT_EQ(second.vUniqueTable.inserts, first.vUniqueTable.inserts);
    EXPECT_EQ(second.vUniqueTable.numEntries, first.vUniqueTable.numEntries);
}

TEST_F(Statistics, HitRatioAndColRatioWithinRange) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    EXPECT_GE(stats.vUniqueTable.hitRatio(), 0.);
    EXPECT_LE(stats.vUniqueTable.hitRatio(), 1.);
    EXPECT_GE(stats.vUniqueTable.colRatio(), 0.);
    pkg.decRef(state);
}

TEST_F(Statistics, LoadFactorEqualsEntriesOverBuckets) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    ASSERT_GT(stats.vUniqueTable.numBuckets, 0U);
    EXPECT_GE(stats.vUniqueTable.numBuckets, stats.vUniqueTable.numEntries);
    EXPECT_DOUBLE_EQ(stats.vUniqueTable.loadFactor(), static_cast<double>(stats.vUniqueTable.numEntries) /
                                                          static_cast<double>(stats.vUniqueTable.numBuckets));
    pkg.decRef(state);
}

TEST_F(Statistics, CollisionsAreBoundedByEntries) {
    // Pins the current-state definition documented on TableStatistics: each
    // entry can contribute at most one to the per-bucket excess, and an
    // empty table has no collisions at all.
    const DwPackage empty(3);
    EXPECT_EQ(empty.statistics().vUniqueTable.collisions, 0U);

    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    EXPECT_LE(stats.vUniqueTable.collisions, stats.vUniqueTable.numEntries);
    pkg.decRef(state);
}

// ---------------------------------------------------------------------
// Reference counting and garbage collection
// ---------------------------------------------------------------------

TEST_F(Statistics, IncRefTracksActiveEntries) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    EXPECT_GT(stats.vUniqueTable.numActiveEntries, 0U);
    EXPECT_GE(stats.vUniqueTable.peakNumActiveEntries, stats.vUniqueTable.numActiveEntries);
    pkg.decRef(state);
}

TEST_F(Statistics, DecRefUntracksActiveEntries) {
    DwPackage pkg(4);
    auto state = ghz(pkg, 4);
    const auto peak = pkg.statistics().vUniqueTable.peakNumActiveEntries;
    ASSERT_GT(peak, 0U);

    pkg.decRef(state);
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.vUniqueTable.numActiveEntries, 0U);
    EXPECT_EQ(stats.vUniqueTable.peakNumActiveEntries, peak);
}

TEST_F(Statistics, GarbageCollectSyncsActiveEntriesToNumEntries) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    ASSERT_TRUE(pkg.garbageCollect(true));
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.vUniqueTable.numActiveEntries, stats.vUniqueTable.numEntries);
    EXPECT_EQ(stats.mUniqueTable.numActiveEntries, stats.mUniqueTable.numEntries);
    pkg.decRef(state);
}

TEST_F(Statistics, GarbageCollectCountsOnlyRealSweeps) {
    DwPackage pkg(3);
    const auto state = ghz(pkg, 3);
    // Well below the initial gc limit, so this call returns without sweeping.
    EXPECT_FALSE(pkg.garbageCollect(false));
    EXPECT_EQ(pkg.statistics().vUniqueTable.gcRuns, 0U);

    pkg.garbageCollect(true);
    EXPECT_EQ(pkg.statistics().vUniqueTable.gcRuns, 1U);
    pkg.decRef(state);
}

// ---------------------------------------------------------------------
// Memory managers
// ---------------------------------------------------------------------

TEST_F(Statistics, FreshMemoryManagerIsSeededAndUnused) {
    const DwPackage pkg(3);
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.vMemory.numAllocations, 1U);
    EXPECT_EQ(stats.vMemory.numAllocated, MemoryManager<DwVNode>::kInitialAllocationSize);
    EXPECT_EQ(stats.vMemory.numUsed, 0U);
    EXPECT_EQ(stats.vMemory.peakNumUsed, 0U);
    EXPECT_EQ(stats.json()["vector"]["memory_manager"], "unused");
}

TEST_F(Statistics, UsedEntriesMatchUniqueTableEntries) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    // Every get() in makeVEdge/makeMEdge is immediately followed by an
    // emplace, and nothing has been returned yet.
    EXPECT_EQ(stats.vMemory.numUsed, stats.vUniqueTable.numEntries);
    EXPECT_EQ(stats.mMemory.numUsed, stats.mUniqueTable.numEntries);
    pkg.decRef(state);
}

TEST_F(Statistics, AllocatesSecondChunkPastInitialSize) {
    DwPackage pkg(12, NormalizationStrategy::None);
    const auto state = manyDistinctNodes(pkg, 12);
    const auto stats = pkg.statistics();
    ASSERT_GT(stats.vMemory.numUsed, MemoryManager<DwVNode>::kInitialAllocationSize);
    EXPECT_GE(stats.vMemory.numAllocations, 2U);
    EXPECT_GT(stats.vMemory.numAllocated, MemoryManager<DwVNode>::kInitialAllocationSize);
    EXPECT_GT(stats.vMemory.getAllocatedMemoryMiB(), 0.);
    (void)state;
}

TEST_F(Statistics, GarbageCollectReturnsEntriesForReuse) {
    DwPackage pkg(6);
    // Discarding ghz's incRef'd handle keeps its final root alive, but every
    // gate DD and superseded intermediate state along the way is already
    // unreferenced -- plenty for the sweep to reclaim.
    (void)ghz(pkg, 6);
    const auto before = pkg.statistics();
    ASSERT_GT(before.vMemory.numUsed, 0U);

    pkg.garbageCollect(true);
    const auto after = pkg.statistics();
    EXPECT_GT(after.vMemory.numAvailableForReuse, 0U);
    EXPECT_LT(after.vMemory.numUsed, before.vMemory.numUsed);
    EXPECT_GE(after.vMemory.peakNumAvailableForReuse, after.vMemory.numAvailableForReuse);
    EXPECT_EQ(after.vMemory.peakNumUsed, before.vMemory.numUsed);
}

TEST_F(Statistics, NewNodesReuseFreedEntries) {
    DwPackage pkg(6);
    (void)ghz(pkg, 6);
    pkg.garbageCollect(true);
    const auto before = pkg.statistics();
    ASSERT_GT(before.vMemory.numAvailableForReuse, 0U);

    (void)ghz(pkg, 6);
    const auto after = pkg.statistics();
    EXPECT_LT(after.vMemory.numAvailableForReuse, before.vMemory.numAvailableForReuse);
    EXPECT_EQ(after.vMemory.numAllocations, before.vMemory.numAllocations);
}

TEST_F(Statistics, UsageRatioAndMemoryAreConsistent) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto &mem = pkg.statistics().vMemory;
    EXPECT_EQ(mem.getTotalNumAvailable(), mem.numAllocated - mem.numUsed);
    EXPECT_EQ(mem.getNumAvailableFromChunks(), mem.getTotalNumAvailable() - mem.numAvailableForReuse);
    EXPECT_GT(mem.getUsageRatio(), 0.);
    EXPECT_LE(mem.getUsageRatio(), 1.);
    EXPECT_LE(mem.getUsedMemoryMiB(), mem.getAllocatedMemoryMiB());
    EXPECT_LE(mem.getUsedMemoryMiB(), mem.getPeakUsedMemoryMiB());
    pkg.decRef(state);
}

// ---------------------------------------------------------------------
// Compute tables
// ---------------------------------------------------------------------

TEST_F(Statistics, ComputeTableTracksLookupsAndHits) {
    DwPackage pkg(4);
    const auto state = pkg.makeZeroState();
    const auto gate = pkg.makeSingleQubitGateDD(0, gates::h());
    (void)pkg.multiply(gate, state);
    const auto first = pkg.statistics();
    ASSERT_GT(first.mvComputeTable.inserts, 0U);

    (void)pkg.multiply(gate, state);
    const auto second = pkg.statistics();
    EXPECT_GT(second.mvComputeTable.lookups, first.mvComputeTable.lookups);
    EXPECT_GT(second.mvComputeTable.hits, 0U);
}

TEST_F(Statistics, MatrixComputeTableUnusedUntilMatrixMultiply) {
    DwPackage pkg(4);
    const auto state = pkg.makeZeroState();
    (void)pkg.multiply(pkg.makeSingleQubitGateDD(0, gates::h()), state);
    const auto stats = pkg.statistics();
    EXPECT_EQ(stats.mmComputeTable.lookups, 0U);
    EXPECT_EQ(stats.json()["compute_tables"]["matrix_matrix_mult"], "unused");
    EXPECT_NE(stats.json()["compute_tables"]["matrix_vector_mult"], "unused");
}

TEST_F(Statistics, ComputeTableCountersSurviveCacheClear) {
    DwPackage pkg(6);
    (void)ghz(pkg, 6);
    const auto before = pkg.statistics();
    ASSERT_GT(before.mvComputeTable.lookups, 0U);
    ASSERT_GT(before.mvComputeTable.inserts, 0U);

    // Sweeps, and therefore clears mvCache_ -- but deliberately leaves the
    // lifetime counters alone (see DwPackage::garbageCollect).
    ASSERT_TRUE(pkg.garbageCollect(true));
    const auto after = pkg.statistics();
    EXPECT_EQ(after.mvComputeTable.lookups, before.mvComputeTable.lookups);
    EXPECT_EQ(after.mvComputeTable.inserts, before.mvComputeTable.inserts);
    EXPECT_EQ(after.mvComputeTable.numEntries, 0U);
    EXPECT_GE(after.mvComputeTable.peakNumEntries, before.mvComputeTable.numEntries);
}

// ---------------------------------------------------------------------
// Report / JSON
// ---------------------------------------------------------------------

TEST_F(Statistics, JsonHasExpectedTopLevelSections) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto j = pkg.statistics().json();
    EXPECT_TRUE(j.contains("vector"));
    EXPECT_TRUE(j.contains("matrix"));
    EXPECT_TRUE(j.contains("compute_tables"));
    EXPECT_TRUE(j.contains("active_memory_mib"));
    EXPECT_TRUE(j.contains("peak_memory_mib"));
    EXPECT_FALSE(j.contains("statistics_disabled"));
    pkg.decRef(state);
}

TEST_F(Statistics, JsonVectorSectionHasBothSubsections) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto j = pkg.statistics().json();
    EXPECT_TRUE(j["vector"].contains("unique_table"));
    EXPECT_TRUE(j["vector"].contains("memory_manager"));
    EXPECT_TRUE(j["matrix"].contains("unique_table"));
    EXPECT_TRUE(j["matrix"].contains("memory_manager"));
    EXPECT_TRUE(j["vector"]["unique_table"].contains("hits"));
    EXPECT_TRUE(j["vector"]["unique_table"].contains("gc_runs"));
    pkg.decRef(state);
}

TEST_F(Statistics, ToStringIsParseablePrettyJson) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    const std::string text = stats.toString();
    EXPECT_NE(text.find('\n'), std::string::npos); // dump(2), not dump()
    EXPECT_NO_THROW((void)nlohmann::json::parse(text));
    EXPECT_EQ(nlohmann::json::parse(text), stats.json());
    pkg.decRef(state);
}

TEST_F(Statistics, StreamOperatorMatchesToString) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    std::ostringstream oss;
    oss << stats;
    EXPECT_EQ(oss.str(), stats.toString());
    EXPECT_EQ(pkg.statisticsString(), stats.toString());
    pkg.decRef(state);
}

TEST_F(Statistics, ActiveMemoryDoesNotExceedPeakMemory) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto stats = pkg.statistics();
    EXPECT_GE(stats.activeMemoryMiB(), 0.);
    EXPECT_LE(stats.activeMemoryMiB(), stats.peakMemoryMiB());
    pkg.decRef(state);
}

TEST_F(Statistics, ResetClearsRunningCountersButKeepsPeaks) {
    // reset() exists for the Statistics base-class contract and MQT Core
    // parity; nothing inside exact-dd calls it. Operates on the returned
    // snapshot, never on the package.
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    auto stats = pkg.statistics();
    const auto peakEntries = stats.vUniqueTable.peakNumEntries;
    const auto peakActive = stats.vUniqueTable.peakNumActiveEntries;
    const auto lookups = stats.vUniqueTable.lookups;
    ASSERT_GT(peakEntries, 0U);

    stats.reset();
    EXPECT_EQ(stats.vUniqueTable.numEntries, 0U);
    EXPECT_EQ(stats.vUniqueTable.numActiveEntries, 0U);
    EXPECT_EQ(stats.vUniqueTable.peakNumEntries, peakEntries);
    EXPECT_EQ(stats.vUniqueTable.peakNumActiveEntries, peakActive);
    EXPECT_EQ(stats.vUniqueTable.lookups, lookups);
    pkg.decRef(state);
}

TEST_F(Statistics, SnapshotDoesNotMutatePackage) {
    DwPackage pkg(4);
    const auto state = ghz(pkg, 4);
    const auto first = pkg.statistics().json();
    const auto second = pkg.statistics().json();
    EXPECT_EQ(first, second);
    pkg.decRef(state);
}

// ---------------------------------------------------------------------
// Disabled build
// ---------------------------------------------------------------------

TEST(StatisticsDisabled, ReportsMarkerWhenCompiledOut) {
    if constexpr (kStatisticsEnabled) {
        GTEST_SKIP() << "built with EXACT_DD_STATISTICS=ON";
    } else {
        DwPackage pkg(4);
        const auto state = ghz(pkg, 4);
        const auto j = pkg.statistics().json();
        EXPECT_TRUE(j["statistics_disabled"].get<bool>());
        EXPECT_TRUE(j.contains("note"));
        EXPECT_FALSE(j.contains("vector"));
        // Node counting is independent of the option.
        EXPECT_GT(pkg.vNodeCount(), 0U);
        pkg.decRef(state);
    }
}
