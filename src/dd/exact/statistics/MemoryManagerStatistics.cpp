#include "dd/exact/statistics/MemoryManagerStatistics.hpp"

#include "dd/exact/DwNode.hpp"

#include <cstddef>
#include <nlohmann/json.hpp>

namespace dd::exact {

template <typename T> std::size_t MemoryManagerStatistics<T>::getNumAvailableFromChunks() const noexcept {
    return getTotalNumAvailable() - numAvailableForReuse;
}

template <typename T> std::size_t MemoryManagerStatistics<T>::getTotalNumAvailable() const noexcept {
    return numAllocated - numUsed;
}

template <typename T> double MemoryManagerStatistics<T>::getUsageRatio() const noexcept {
    if (numAllocated == 0U) {
        return 0.;
    }
    return static_cast<double>(numUsed) / static_cast<double>(numAllocated);
}

template <typename T> double MemoryManagerStatistics<T>::getAllocatedMemoryMiB() const noexcept {
    return static_cast<double>(numAllocated) * kEntryMemoryMiB;
}

template <typename T> double MemoryManagerStatistics<T>::getUsedMemoryMiB() const noexcept {
    return static_cast<double>(numUsed) * kEntryMemoryMiB;
}

template <typename T> double MemoryManagerStatistics<T>::getPeakUsedMemoryMiB() const noexcept {
    return static_cast<double>(peakNumUsed) * kEntryMemoryMiB;
}

template <typename T> void MemoryManagerStatistics<T>::reset() noexcept {
    numAllocations = 0U;
    numAllocated = 0U;
    numUsed = 0U;
    numAvailableForReuse = 0U;
}

template <typename T> nlohmann::json MemoryManagerStatistics<T>::json() const {
    if (peakNumUsed == 0U) {
        return "unused";
    }
    auto j = Statistics::json();
    j["memory_allocated_MiB"] = getAllocatedMemoryMiB();
    j["memory_used_MiB"] = getUsedMemoryMiB();
    j["memory_used_MiB_peak"] = getPeakUsedMemoryMiB();
    j["num_allocated"] = numAllocated;
    j["num_allocations"] = numAllocations;
    j["num_available_for_reuse"] = numAvailableForReuse;
    j["num_available_for_reuse_peak"] = peakNumAvailableForReuse;
    j["num_available_from_chunks"] = getNumAvailableFromChunks();
    j["num_available_total"] = getTotalNumAvailable();
    j["num_used"] = numUsed;
    j["num_used_peak"] = peakNumUsed;
    j["usage_ratio"] = getUsageRatio();
    return j;
}

// The definitions above live here, rather than in the header, so that
// <nlohmann/json.hpp> stays out of MemoryManager.hpp's -- and therefore
// DwPackage.hpp's -- include chain. The trade-off is that a MemoryManager
// over any node type other than these two will not link until its instantiation is added below.
template struct MemoryManagerStatistics<DwVNode>;
template struct MemoryManagerStatistics<DwMNode>;

} // namespace dd::exact
