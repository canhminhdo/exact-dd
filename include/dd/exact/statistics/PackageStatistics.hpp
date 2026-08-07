#ifndef DD_EXACT_PACKAGE_STATISTICS_HPP
#define DD_EXACT_PACKAGE_STATISTICS_HPP

#include "dd/exact/DwNode.hpp"
#include "dd/exact/statistics/MemoryManagerStatistics.hpp"
#include "dd/exact/statistics/Statistics.hpp"
#include "dd/exact/statistics/TableStatistics.hpp"
#include "dd/exact/statistics/UniqueTableStatistics.hpp"

#include <nlohmann/json_fwd.hpp>

namespace dd::exact {

/**
 * A snapshot of every statistic DwPackage collects, returned by value from
 * DwPackage::statistics().
 *
 * 1. A MEMBER, NOT FREE FUNCTIONS.
 *    Every member of DwPackage is private
 *    (vUnique_, mUnique_, mvCache_, mmCache_, vMemory_, mMemory_), and
 *    neither making them public nor befriending a free function buys
 *    anything over simply asking the package for its own report -- which is
 *    also how the rest of DwPackage's introspection surface already reads
 *    (vNodeCount(), vectorToString(), printVectorDiagram()).
 *
 * 2. A SNAPSHOT, NOT REFERENCES. That is what lets DwPackage::statistics()
 *    be const: the map-derived fields (numEntries, numBuckets, collisions --
 *    see TableStatistics::snapshot()) are refreshed on the returned copy, so
 *    nothing in the package is mutated and no member has to be mutable. The
 *    copy is six small structs of counters; it is a report-time cost only.
 *
 * One further reading note: garbageCollect() clears both compute-table
 * caches without resetting their counters, so a compute table's numEntries
 * reflects the cache since the last collection while its lookups/hits/
 * inserts are lifetime totals.
 */
struct PackageStatistics : public Statistics {
    UniqueTableStatistics vUniqueTable;
    UniqueTableStatistics mUniqueTable;
    MemoryManagerStatistics<DwVNode> vMemory;
    MemoryManagerStatistics<DwMNode> mMemory;
    /// Matrix-vector multiplication cache (DwPackage's mvCache_).
    TableStatistics mvComputeTable;
    /// Matrix-matrix multiplication cache (DwPackage's mmCache_).
    TableStatistics mmComputeTable;

    /// Memory held by nodes currently reachable from an external root,
    /// counting node and edge storage. A lower bound -- see the class note.
    [[nodiscard]] double activeMemoryMiB() const noexcept;
    /// Memory corresponding to each memory manager's high-water mark of
    /// entries in use. A lower bound -- see the class note.
    [[nodiscard]] double peakMemoryMiB() const noexcept;

    void reset() noexcept override;

    /// The combined report. When built with EXACT_DD_STATISTICS=OFF this
    /// emits only a "statistics_disabled" marker rather than the usual key
    /// tree full of zeros.
    [[nodiscard]] nlohmann::json json() const override;
};

} // namespace dd::exact

#endif // DD_EXACT_PACKAGE_STATISTICS_HPP
