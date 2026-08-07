#ifndef DD_EXACT_MEMORY_MANAGER_STATISTICS_HPP
#define DD_EXACT_MEMORY_MANAGER_STATISTICS_HPP

#include "dd/exact/statistics/Statistics.hpp"
#include "dd/exact/statistics/StatisticsConfig.hpp"

#include <algorithm>
#include <cstddef>
#include <nlohmann/json_fwd.hpp>

namespace dd::exact {

/**
 * Statistics for one pooled node allocator. The only templated record in this module:
 * entrySize and kEntryMemoryMiB come straight from sizeof(T).
 *
 * As with TableStatistics, every MiB figure here is a LOWER BOUND for
 * exact-dd: each node holds one Dw per outgoing edge, and
 * a Dw's five arbitrary-precision integers own heap storage that
 * sizeof(T) cannot see.
 *
 * The non-inline members are defined in MemoryManagerStatistics.cpp with
 * explicit instantiations for DwVNode/DwMNode, which is what keeps
 * <nlohmann/json.hpp> out of MemoryManager.hpp's -- and therefore
 * DwPackage.hpp's -- include chain.
 */
template <typename T> struct MemoryManagerStatistics : public Statistics {
    std::size_t entrySize = sizeof(T);
    /// Number of chunk allocations, including the initial one.
    std::size_t numAllocations = 0U;
    /// Total number of entries carved out across all chunks.
    std::size_t numAllocated = 0U;
    std::size_t numUsed = 0U;
    std::size_t numAvailableForReuse = 0U;
    std::size_t peakNumUsed = 0U;
    std::size_t peakNumAvailableForReuse = 0U;

    static constexpr double kEntryMemoryMiB = static_cast<double>(sizeof(T)) / static_cast<double>(1ULL << 20U);

    /// Records the chunk allocated by MemoryManager's constructor.
    void trackInitialAllocation(std::size_t numEntries) noexcept {
        if constexpr (kStatisticsEnabled) {
            numAllocations = 1U;
            numAllocated = numEntries;
        }
    }

    /// Records a further chunk from MemoryManager::allocateNewChunk().
    void trackAllocation(std::size_t numEntries) noexcept {
        if constexpr (kStatisticsEnabled) {
            ++numAllocations;
            numAllocated += numEntries;
        }
    }

    /// Records entries handed out fresh from a chunk.
    void trackUsedEntries(std::size_t numEntries = 1U) noexcept {
        if constexpr (kStatisticsEnabled) {
            numUsed += numEntries;
            peakNumUsed = std::max(peakNumUsed, numUsed);
        }
    }

    /// Records entries handed out from the free list rather than a chunk.
    void trackReusedEntries(std::size_t numEntries = 1U) noexcept {
        if constexpr (kStatisticsEnabled) {
            numUsed += numEntries;
            peakNumUsed = std::max(peakNumUsed, numUsed);
            numAvailableForReuse -= numEntries;
        }
    }

    /// Records an entry returned to the free list.
    void trackReturnedEntry() noexcept {
        if constexpr (kStatisticsEnabled) {
            ++numAvailableForReuse;
            peakNumAvailableForReuse = std::max(peakNumAvailableForReuse, numAvailableForReuse);
            if (numUsed > 0U) {
                --numUsed;
            }
        }
    }

    /// Entries still untouched in the current chunk, i.e. available without
    /// going through the free list.
    [[nodiscard]] std::size_t getNumAvailableFromChunks() const noexcept;
    /// Every entry not currently in use, from chunks and free list alike.
    [[nodiscard]] std::size_t getTotalNumAvailable() const noexcept;
    [[nodiscard]] double getUsageRatio() const noexcept;
    [[nodiscard]] double getAllocatedMemoryMiB() const noexcept;
    [[nodiscard]] double getUsedMemoryMiB() const noexcept;
    [[nodiscard]] double getPeakUsedMemoryMiB() const noexcept;

    void reset() noexcept override;

    /// Returns the JSON string "unused" for a manager that never handed out
    /// an entry, mirroring MQT Core.
    [[nodiscard]] nlohmann::json json() const override;
};

} // namespace dd::exact

#endif // DD_EXACT_MEMORY_MANAGER_STATISTICS_HPP
