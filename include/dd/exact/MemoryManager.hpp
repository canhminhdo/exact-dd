#ifndef DD_EXACT_MEMORY_MANAGER_HPP
#define DD_EXACT_MEMORY_MANAGER_HPP

#include "dd/exact/DwNode.hpp"
#include "dd/exact/statistics/MemoryManagerStatistics.hpp"

#include <cstddef>
#include <type_traits>
#include <vector>

namespace dd::exact {

/**
 * A pooled bump allocator for DD node objects of type T, ported from MQT
 * Core's dd::MemoryManager<T>. Objects are stored in contiguous chunks that
 * grow by kGrowthFactor when exhausted. Freed objects are recycled via a
 * singly linked free list built on T's own `next` member (the same field
 * MQT reuses for this purpose), avoiding per-node heap allocation and
 * deallocation for the high churn of node creation/collection typical of a
 * decision-diagram package.
 *
 * Pointers previously returned by get() remain valid even after the chunks_
 * vector grows: growing chunks_ moves (never copies) its std::vector<T>
 * elements, and moving a std::vector preserves the validity of pointers
 * into its element storage.
 *
 * Allocation and reuse are counted into a MemoryManagerStatistics<T>, which
 * statistics() exposes and DwPackage folds into its own report. Every one
 * of those tracking calls compiles away when EXACT_DD_STATISTICS is OFF.
 */
template <typename T> class MemoryManager {
    static_assert(std::is_same_v<decltype(T::next), T *>, "T must have a `next` member of type T*");
    static_assert(std::is_same_v<decltype(T::ref), RefCount>, "T must have a `ref` member of type RefCount");

public:
    /// The number of entries allocated in the first chunk.
    static constexpr std::size_t kInitialAllocationSize = 2048U;
    /// The factor by which chunk size grows each time a new chunk is needed.
    static constexpr std::size_t kGrowthFactor = 2U;

    explicit MemoryManager(std::size_t initialAllocationSize = kInitialAllocationSize)
        : chunks_(1, std::vector<T>(initialAllocationSize)), chunkIt_(chunks_[0].begin()),
          chunkEndIt_(chunks_[0].end()) {
        stats_.trackInitialAllocation(initialAllocationSize);
    }

    // available_/chunkIt_/chunkEndIt_ point into this object's own chunks_;
    // copying would leave those pointers/iterators referring to the
    // original's storage. Moving is safe since chunks_ is only ever moved
    // (never copied) along with them.
    MemoryManager(const MemoryManager &) = delete;
    MemoryManager &operator=(const MemoryManager &) = delete;
    MemoryManager(MemoryManager &&) = default;
    MemoryManager &operator=(MemoryManager &&) = default;

    /// Returns a pointer to an entry ready for (re-)use, with `ref` reset to
    /// 0. Reuses a previously returned entry if one is available; otherwise
    /// takes the next unused slot from the current chunk, allocating a new
    /// (larger) chunk first if the current one is exhausted.
    [[nodiscard]] T *get() {
        T *entry = nullptr;
        if (available_ != nullptr) {
            entry = available_;
            available_ = available_->next;
            stats_.trackReusedEntries();
        } else {
            if (chunkIt_ == chunkEndIt_) {
                allocateNewChunk();
            }
            entry = &*chunkIt_;
            ++chunkIt_;
            stats_.trackUsedEntries();
        }
        entry->ref = 0;
        return entry;
    }

    /// Returns an entry (which must have ref == 0) to the free list so a
    /// later get() can reuse it. The entry must not be used after this call
    /// until it is handed out again by get().
    void returnEntry(T *entry) noexcept {
        entry->next = available_;
        available_ = entry;
        stats_.trackReturnedEntry();
    }

    /// Allocation and reuse counters for this manager. All zero when built
    /// with EXACT_DD_STATISTICS=OFF.
    [[nodiscard]] const MemoryManagerStatistics<T> &statistics() const noexcept { return stats_; }

private:
    void allocateNewChunk() {
        const std::size_t nextSize = chunks_.back().size() * kGrowthFactor;
        chunks_.emplace_back(nextSize);
        chunkIt_ = chunks_.back().begin();
        chunkEndIt_ = chunks_.back().end();
        stats_.trackAllocation(nextSize);
    }

    /// Singly linked list (via T::next) of entries available for reuse.
    T *available_{};
    /// Contiguous chunks of storage; growing this never invalidates pointers
    /// into the (moved, not copied) std::vector<T> chunks it holds.
    std::vector<std::vector<T>> chunks_;
    typename std::vector<T>::iterator chunkIt_;
    typename std::vector<T>::iterator chunkEndIt_;

    /// Note that its numAvailableForReuse is the only way to know how long
    /// the free list is: available_ is a bare singly linked list threaded
    /// through T::next, with no size of its own.
    MemoryManagerStatistics<T> stats_{};
};

} // namespace dd::exact

#endif // DD_EXACT_MEMORY_MANAGER_HPP
