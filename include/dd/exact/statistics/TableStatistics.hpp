#ifndef DD_EXACT_TABLE_STATISTICS_HPP
#define DD_EXACT_TABLE_STATISTICS_HPP

#include "dd/exact/statistics/Statistics.hpp"
#include "dd/exact/statistics/StatisticsConfig.hpp"

#include <algorithm>
#include <cstddef>
#include <nlohmann/json_fwd.hpp>

namespace dd::exact {

/**
 * Statistics for one hash table, used here for
 * DwPackage's two compute-table caches (mvCache_/mmCache_); the unique
 * tables extend it via UniqueTableStatistics.
 *
 * Every MiB figure here is a LOWER BOUND for exact-dd
 * each node holds one Dw per outgoing edge, and
 * a Dw's five arbitrary-precision integers
 * own heap storage that sizeof(T) cannot see.
 */
struct TableStatistics : public Statistics {
    /// Size of a single table entry in bytes; seeded once by DwPackage's ctor.
    std::size_t entrySize = 0U;
    std::size_t numBuckets = 0U;
    std::size_t numEntries = 0U;
    std::size_t peakNumEntries = 0U;
    std::size_t collisions = 0U;
    std::size_t hits = 0U;
    std::size_t lookups = 0U;
    std::size_t inserts = 0U;

    /// Records one table probe, whatever its outcome.
    void trackLookup() noexcept {
        if constexpr (kStatisticsEnabled) {
            ++lookups;
        }
    }

    /// Records a probe that found its key.
    void trackHit() noexcept {
        if constexpr (kStatisticsEnabled) {
            ++hits;
        }
    }

    /**
     * Records a newly inserted entry.
     *
     * Taking the map by reference rather than a precomputed size keeps the
     * call site strictly free when kStatisticsEnabled is false: nothing is
     * evaluated there, only a reference formed.
     */
    template <class Map> void trackInsert(const Map &map) noexcept {
        if constexpr (kStatisticsEnabled) {
            ++inserts;
            numEntries = map.size();
            peakNumEntries = std::max(peakNumEntries, numEntries);
        }
    }

    /**
     * Refreshes every field that is a property of the map's current state
     * rather than of the history of operations on it: numEntries, numBuckets
     * and collisions (see the collision-semantics note on this class).
     */
    template <class Map> void snapshot(const Map &map) noexcept {
        if constexpr (kStatisticsEnabled) {
            numEntries = map.size();
            peakNumEntries = std::max(peakNumEntries, numEntries);
            numBuckets = map.bucket_count();
            collisions = 0U;
            for (std::size_t b = 0; b < numBuckets; ++b) {
                const std::size_t n = map.bucket_size(b);
                if (n > 1U) {
                    collisions += n - 1U;
                }
            }
        }
    }

    /// Zeroes numEntries only, the cumulative hits/lookups/collisions/inserts
    /// history deliberately survives, as do the peak values.
    void reset() noexcept override;

    [[nodiscard]] double hitRatio() const noexcept;
    [[nodiscard]] double colRatio() const noexcept;
    [[nodiscard]] double loadFactor() const noexcept;
    [[nodiscard]] double getEntrySizeMiB() const noexcept;
    [[nodiscard]] double getMemoryMiB() const noexcept;

    /// Returns the JSON string "unused" for a table that was never probed
    [[nodiscard]] nlohmann::json json() const override;
};

} // namespace dd::exact

#endif // DD_EXACT_TABLE_STATISTICS_HPP
