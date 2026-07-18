include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(EXACT_DD_CMAKE_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/exact-dd")

install(
    TARGETS ExactDD
    EXPORT exact-dd-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
    EXPORT exact-dd-targets
    FILE exact-dd-targets.cmake
    NAMESPACE ExactDD::
    DESTINATION ${EXACT_DD_CMAKE_INSTALL_DIR}
)

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/exact-dd-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/exact-dd-config.cmake
    INSTALL_DESTINATION ${EXACT_DD_CMAKE_INSTALL_DIR}
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/exact-dd-config-version.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

install(
    FILES
    ${CMAKE_CURRENT_BINARY_DIR}/exact-dd-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/exact-dd-config-version.cmake
    DESTINATION ${EXACT_DD_CMAKE_INSTALL_DIR}
)
if(EXACT_DD_WITH_GMP)
    install(FILES ${PROJECT_SOURCE_DIR}/cmake/FindGMP.cmake DESTINATION ${EXACT_DD_CMAKE_INSTALL_DIR})
endif()
