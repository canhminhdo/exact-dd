#include "dd/exact/statistics/PackageStatistics.hpp"

#include "dd/exact/DwNode.hpp"
#include "dd/exact/DwPackage.hpp"
#include "dd/exact/statistics/StatisticsConfig.hpp"

#include <nlohmann/json.hpp>

namespace dd::exact {

namespace {
// These are lower bounds for object sizes: a Dw's arbitrary-precision limbs live on the
// heap, outside anything sizeof() can measure. See PackageStatistics' note.
constexpr double kMiB = static_cast<double>(1ULL << 20U);
constexpr double kVNodeMemoryMiB = static_cast<double>(sizeof(DwVNode)) / kMiB;
constexpr double kMNodeMemoryMiB = static_cast<double>(sizeof(DwMNode)) / kMiB;
constexpr double kVEdgeMemoryMiB = static_cast<double>(sizeof(DwEdge<DwVNode>)) / kMiB;
constexpr double kMEdgeMemoryMiB = static_cast<double>(sizeof(DwEdge<DwMNode>)) / kMiB;
} // namespace

double PackageStatistics::activeMemoryMiB() const noexcept {
    const auto v = static_cast<double>(vUniqueTable.numActiveEntries) * (kVNodeMemoryMiB + kVEdgeMemoryMiB);
    const auto m = static_cast<double>(mUniqueTable.numActiveEntries) * (kMNodeMemoryMiB + kMEdgeMemoryMiB);
    return v + m;
}

double PackageStatistics::peakMemoryMiB() const noexcept {
    const auto v = static_cast<double>(vMemory.peakNumUsed) * (kVNodeMemoryMiB + kVEdgeMemoryMiB);
    const auto m = static_cast<double>(mMemory.peakNumUsed) * (kMNodeMemoryMiB + kMEdgeMemoryMiB);
    return v + m;
}

void PackageStatistics::reset() noexcept {
    vUniqueTable.reset();
    mUniqueTable.reset();
    vMemory.reset();
    mMemory.reset();
    mvComputeTable.reset();
    mmComputeTable.reset();
}

nlohmann::json PackageStatistics::json() const {
    nlohmann::json j;

    if constexpr (!kStatisticsEnabled) {
        // Deliberately NOT the usual key tree filled with zeros: a consumer
        // indexing into j["vector"]["unique_table"] should fail loudly here
        // rather than read a zero and believe it.
        j["statistics_disabled"] = true;
        j["note"] = "exact-dd was built with EXACT_DD_STATISTICS=OFF; no counters were collected. "
                    "Reconfigure with -DEXACT_DD_STATISTICS=ON.";
        return j;
    }

    auto &vector = j["vector"];
    vector["unique_table"] = vUniqueTable.json();
    vector["memory_manager"] = vMemory.json();

    auto &matrix = j["matrix"];
    matrix["unique_table"] = mUniqueTable.json();
    matrix["memory_manager"] = mMemory.json();

    auto &computeTables = j["compute_tables"];
    computeTables["matrix_vector_mult"] = mvComputeTable.json();
    computeTables["matrix_matrix_mult"] = mmComputeTable.json();

    j["active_memory_mib"] = activeMemoryMiB();
    j["peak_memory_mib"] = peakMemoryMiB();

    return j;
}

// Defined here rather than in DwPackage.cpp so that the latter never needs
// the full <nlohmann/json.hpp>; a member's definition may live in any
// translation unit, and this one already includes everything required.
PackageStatistics DwPackage::statistics() const {
    PackageStatistics s;
    s.vUniqueTable = vUniqueStats_;
    s.mUniqueTable = mUniqueStats_;
    s.vMemory = vMemory_.statistics();
    s.mMemory = mMemory_.statistics();
    s.mvComputeTable = mvCacheStats_;
    s.mmComputeTable = mmCacheStats_;

    // Refresh the map-derived fields on the COPY -- this is what keeps this
    // method const and leaves the package unmutated. See
    // TableStatistics::snapshot() for what "map-derived" covers and why
    // collisions in particular are computed here rather than tracked.
    s.vUniqueTable.snapshot(vUnique_);
    s.mUniqueTable.snapshot(mUnique_);
    s.mvComputeTable.snapshot(mvCache_);
    s.mmComputeTable.snapshot(mmCache_);
    return s;
}

} // namespace dd::exact
