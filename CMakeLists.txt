cmake_minimum_required(VERSION 3.18)
project(EDLib LANGUAGES C CXX
        VERSION 2.0)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CXXStandard)
use_cxx17()

option(EDLIB_WITH_ALPSCORE
    "Pull EDLib-alpscore (legacy ALPSCore-based API: alps::params / alps::hdf5 / alps::gf) via FetchContent. Off = header-only edlib:: core only."
    ON)

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Eigen3 REQUIRED NO_MODULE)
include_directories(${Boost_INCLUDE_DIRS})

# cpp-arnoldi: header-only C++ Arnoldi/Lanczos eigensolver (replaces ARPACK).
include(FetchContent)
FetchContent_Declare(cpp_arnoldi
    GIT_REPOSITORY https://github.com/Q-solvers/cpp-arnoldi.git
    GIT_TAG        444304d88d6029f717031ea5b45294a4bd2a3517)
FetchContent_MakeAvailable(cpp_arnoldi)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
    message("setting build type to default, i.e. optimized build")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS -DNDEBUG -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")

set(extlibs
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    arnoldi::arnoldi)

if(USE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI")
        set(parlibs MPI::MPI_CXX)
    endif()
    message(INFO " Using parallel libraries: ${parlibs}")
endif()

add_subdirectory(include/edlib)
add_subdirectory(include/ext)

# Header purity: include/edlib and include/ext must never
# `#include <alps...>`. Enforce at configure time.
file(GLOB _edlib_core_headers
     "${CMAKE_CURRENT_SOURCE_DIR}/include/edlib/*.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/include/ext/*.h")
foreach(_h IN LISTS _edlib_core_headers)
    file(STRINGS "${_h}" _hits REGEX "^[ \t]*#[ \t]*include[ \t]+[<\"]alps")
    if(_hits)
        message(FATAL_ERROR
            "edlib header purity violation: ${_h} contains an alps include:\n${_hits}")
    endif()
endforeach()

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)
add_library(EDLib::EDLib ALIAS ${PROJECT_NAME})
target_include_directories(${PROJECT_NAME}
    INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
              $<INSTALL_INTERFACE:include>)
target_link_libraries(${PROJECT_NAME} INTERFACE Eigen3::Eigen arnoldi::arnoldi)

# Legacy ALPSCore-based shim layer lives in the separate EDLib-alpscore
# project, fetched here on demand. It defines the EDLib::Alpscore INTERFACE
# target and the EDLIB_WITH_ALPSCORE compile definition.
if(EDLIB_WITH_ALPSCORE)
    FetchContent_Declare(edlib_alpscore
        GIT_REPOSITORY https://github.com/Q-solvers/EDLib-alpscore.git
        GIT_TAG        main)
    FetchContent_MakeAvailable(edlib_alpscore)
    list(APPEND extlibs EDLib::Alpscore)

    # Fold EDLib-alpscore into the parent install: a single
    # `find_package(EDLib)` provides both the core and the alpscore shim.
    set_target_properties(EDLib-alpscore PROPERTIES EXPORT_NAME Alpscore)
endif()

set(_edlib_install_targets ${PROJECT_NAME})
if(EDLIB_WITH_ALPSCORE)
    list(APPEND _edlib_install_targets EDLib-alpscore)
endif()

install(TARGETS ${_edlib_install_targets}
        EXPORT ${PROJECT_NAME}_targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file("${PROJECT_NAME}ConfigVersion.cmake"
                                 VERSION ${PROJECT_VERSION}
                                 COMPATIBILITY SameMajorVersion)

configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
    INSTALL_DESTINATION
    ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(EXPORT ${PROJECT_NAME}_targets
        FILE ${PROJECT_NAME}Targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)

install(FILES include/edlib.h DESTINATION include)

# Core API: always installed (header-only, no alpscore).
install(DIRECTORY include/edlib DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN ".ipynb_checkpoints" EXCLUDE)
install(DIRECTORY include/ext   DESTINATION include
        FILES_MATCHING PATTERN "*.h" PATTERN ".ipynb_checkpoints" EXCLUDE)

# Legacy alpscore shim headers come from the FetchContent'd EDLib-alpscore.
if(EDLIB_WITH_ALPSCORE)
    install(DIRECTORY ${edlib_alpscore_SOURCE_DIR}/include/edlib/alpscore
            DESTINATION include/edlib
            FILES_MATCHING PATTERN "*.h" PATTERN ".ipynb_checkpoints" EXCLUDE)
    install(DIRECTORY ${edlib_alpscore_SOURCE_DIR}/include/ext/alpscore
            DESTINATION include/ext
            FILES_MATCHING PATTERN "*.h" PATTERN ".ipynb_checkpoints" EXCLUDE)
endif()

option(Testing "Enable testing" OFF)
if(Testing)
    include(EnableGtests)
    include_directories(${GTEST_INCLUDE_DIR})
    enable_testing()

    # Core-only tests: build and run regardless of EDLIB_WITH_ALPSCORE.
    add_subdirectory(test/core)
    add_test(HubbardCoreTest    test/core/HubbardCoreTest)
    add_test(LanczosCoreTest    test/core/LanczosCoreTest)
    add_test(SzSymmetryCoreTest test/core/SzSymmetryCoreTest)
    add_test(NSymmetryCoreTest  test/core/NSymmetryCoreTest)

    # Legacy alps-using tests only when the option is enabled.
    if(EDLIB_WITH_ALPSCORE)
        add_subdirectory(test)
        add_test(SzSymmetryTest        test/SzSymmetryTest)
        add_test(NSymmetryTest         test/NSymmetryTest)
        add_test(HubbardModelTest      test/HubbardModelTest)
        add_test(LanczosTest           test/LanczosTest)
        add_test(StaticObservablesTest test/StaticObservablesTest)
    endif()
endif()

option(Examples "Enable examples" ON)
if(Examples)
    add_subdirectory(examples)
    if(EDLIB_WITH_ALPSCORE)
        set(SOURCE_FILES main.cpp)
        add_executable(Hubbard ${SOURCE_FILES})
        target_link_libraries(Hubbard ${PROJECT_NAME} ${extlibs})
        if(USE_MPI)
            target_link_libraries(Hubbard ${parlibs})
        endif()
    endif()
endif()
