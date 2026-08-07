#include "dd/exact/statistics/Statistics.hpp"

#include <nlohmann/json.hpp>
#include <string>

namespace dd::exact {

nlohmann::json Statistics::json() const { return nlohmann::json{}; }

std::string Statistics::toString() const { return json().dump(2U); }

} // namespace dd::exact
