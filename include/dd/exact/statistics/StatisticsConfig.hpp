#ifndef DD_EXACT_STATISTICS_CONFIG_HPP
#define DD_EXACT_STATISTICS_CONFIG_HPP

namespace dd::exact {

/**
 * Compile-time switch for the statistics module (dd/exact/statistics/),
 * driven by the CMake option EXACT_DD_STATISTICS, which sets the PUBLIC
 * compile definition DD_EXACT_STATISTICS
 *
 * Every mutating method on the statistics structs wraps its body in
 * `if constexpr (kStatisticsEnabled)` and is defined inline in its header,
 * so that when the option is OFF the bodies are discarded and the calls at
 * DwPackage's hot sites (the makeVEdge/makeMEdge unique-table probe, the
 * compute-table probe in multiplyRec, incRef/decRef) collapse to nothing.
 * This keeps the call sites free of #ifdefs -- there are over twenty of
 * them -- at no runtime cost.
 *
 * Because the tracking methods are inline, this definition MUST agree
 * between the library and every consumer -- otherwise the inline functions
 * have differing definitions across translation units, which is an ODR
 * violation no compiler will diagnose. That is why the CMake definition is
 * PUBLIC rather than PRIVATE.
 */
#ifdef DD_EXACT_STATISTICS
inline constexpr bool kStatisticsEnabled = true;
#else
inline constexpr bool kStatisticsEnabled = false;
#endif

} // namespace dd::exact

#endif // DD_EXACT_STATISTICS_CONFIG_HPP
