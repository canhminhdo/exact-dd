#include "dd/exact/statistics/UniqueTableStatistics.hpp"

#include <nlohmann/json.hpp>

namespace dd::exact {

void UniqueTableStatistics::reset() noexcept {
    TableStatistics::reset();
    numActiveEntries = 0U;
}

nlohmann::json UniqueTableStatistics::json() const {
    if (lookups == 0U) {
        return "unused";
    }
    auto j = TableStatistics::json();
    j["num_active_entries"] = numActiveEntries;
    j["peak_num_active_entries"] = peakNumActiveEntries;
    j["gc_runs"] = gcRuns;
    return j;
}

} // namespace dd::exact
