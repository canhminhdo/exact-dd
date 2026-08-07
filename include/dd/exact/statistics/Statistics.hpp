#ifndef DD_EXACT_STATISTICS_HPP
#define DD_EXACT_STATISTICS_HPP

#include <nlohmann/json_fwd.hpp>
#include <ostream>
#include <string>

namespace dd::exact {

/**
 * Abstract base of every statistics record in this module
 *
 * JSON is the primary serialization format: subclasses override json() and
 * inherit toString() (a two-space pretty-print of json()) and the stream
 * operator for free. Only <nlohmann/json_fwd.hpp> is included here so that
 * the ~25k-line json.hpp stays out of the include chain of DwPackage.hpp,
 * which reaches these headers via MemoryManager.hpp; the full header is
 * pulled in by the .cpp files that actually build the objects.
 */
struct Statistics {
    Statistics() = default;
    Statistics(const Statistics &) = default;
    Statistics(Statistics &&) = default;
    Statistics &operator=(const Statistics &) = default;
    Statistics &operator=(Statistics &&) = default;
    virtual ~Statistics() = default;

    /// Resets the running statistics, leaving historical peak values intact.
    virtual void reset() noexcept {}

    /// A JSON representation of these statistics.
    [[nodiscard]] virtual nlohmann::json json() const;

    /// A pretty-printed (two-space indented) rendering of json().
    [[nodiscard]] virtual std::string toString() const;

    friend std::ostream &operator<<(std::ostream &os, const Statistics &stats) { return os << stats.toString(); }
};

} // namespace dd::exact

#endif // DD_EXACT_STATISTICS_HPP
