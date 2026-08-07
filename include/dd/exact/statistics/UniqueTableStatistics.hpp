#ifndef DD_EXACT_UNIQUE_TABLE_STATISTICS_HPP
#define DD_EXACT_UNIQUE_TABLE_STATISTICS_HPP

#include "dd/exact/statistics/StatisticsConfig.hpp"
#include "dd/exact/statistics/TableStatistics.hpp"

#include <algorithm>
#include <cstddef>
#include <nlohmann/json_fwd.hpp>

namespace dd::exact {

/**
 * Statistics for one hash-consing unique table (DwPackage's vUnique_ /
 * mUnique_). Adds the reference-counting and garbage-collection counters
 * that a compute-table cache has no use for.
 *
 * The collision- and memory-reporting caveats documented on TableStatistics
 * apply here unchanged.
 */
struct UniqueTableStatistics : public TableStatistics {
    /// Entries whose reference count is greater than zero, i.e. those
    /// garbageCollect() must not reclaim.
    std::size_t numActiveEntries = 0U;
    std::size_t peakNumActiveEntries = 0U;
    /// Sweeps actually performed; garbageCollect()'s early exit is not one.
    std::size_t gcRuns = 0U;

    /**
     * Records a 0 -> 1 reference-count transition (DwPackage::incRef), the
     * moment a node becomes externally reachable.
     *
     * Note that a node whose count saturates at kMaxRefCount is never
     * decremented again, so it stays counted as active forever.
     * trackGcSweep() re-derives the count after every sweep.
     */
    void trackActiveEntry() noexcept {
        if constexpr (kStatisticsEnabled) {
            ++numActiveEntries;
            peakNumActiveEntries = std::max(peakNumActiveEntries, numActiveEntries);
        }
    }

    /// Records a 1 -> 0 reference-count transition (DwPackage::decRef).
    void untrackActiveEntry() noexcept {
        if constexpr (kStatisticsEnabled) {
            if (numActiveEntries > 0U) {
                --numActiveEntries;
            }
        }
    }

    /**
     * Records one completed sweep of `map` by DwPackage::garbageCollect().
     *
     * After a sweep every surviving entry has ref > 0 by construction, so
     * numActiveEntries is re-derived from the map exactly; this also self-heals
     * any drift accumulated by the incRef/decRef tracking above.
     */
    template <class Map> void trackGcSweep(const Map &map) noexcept {
        if constexpr (kStatisticsEnabled) {
            ++gcRuns;
            numEntries = map.size();
            numActiveEntries = numEntries;
        }
    }

    void reset() noexcept override;
    [[nodiscard]] nlohmann::json json() const override;
};

} // namespace dd::exact

#endif // DD_EXACT_UNIQUE_TABLE_STATISTICS_HPP
