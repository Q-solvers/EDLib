cmake_minimum_required(VERSION 3.8.2)
project(EDLib LANGUAGES C CXX Fortran
        VERSION 2.0)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(CXXStandard)
use_cxx11()

find_package(ALPSCore COMPONENTS params hdf5 gf REQUIRED)
find_package(ARPACK REQUIRED)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALPSCore_INCLUDES})


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
  message("setting build type to default, i.e. optimized build")
endif()

message("ARPACK" ${ARPACK_LIB})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS -DNDEBUG -g")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")


set(extlibs
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ALPSCore_LIBRARIES}
    ${ARPACK_LIB}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    )

if(USE_MPI)
    if(!ALPS_HAVE_MPI)
        MESSAGE(FATAL_ERROR "Enable MPI support in ALPSCore.")
    endif(!ALPS_HAVE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI")
    endif(MPI_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")
    set(parlibs ${PARPACK_LIB}
                ${MPI_CXX_LIBRARIES}
                ${MPI_Fortran_LIBRARIES})
    message(INFO " Using parallel libraries: ${parlibs}")
endif(USE_MPI)

# Commented due to the bug in CMake (fixed in 3.13)
#add_subdirectory(lib)
add_subdirectory(include/edlib)
add_subdirectory(include/ext)

include(GNUInstallDirs)

add_library(${PROJECT_NAME} INTERFACE)
add_library(common-lib lib/dummy.f90)
add_library(EDLib::EDLib ALIAS ${PROJECT_NAME})  # Define namespaced alias target
add_library(EDLib::common-lib ALIAS common-lib)  # Define namespaced alias target
target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${${PROJECT_NAME}_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

install(TARGETS ${PROJECT_NAME} common-lib
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

file (GLOB_RECURSE TARGET_INSTALL_LIB_HEADERS "${PROJECT_SOURCE_DIR}/include/edlib/*.h")
file (GLOB_RECURSE TARGET_INSTALL_EXT_HEADERS "${PROJECT_SOURCE_DIR}/include/ext/*.h")
install(FILES include/edlib.h
        DESTINATION include)
install(FILES ${TARGET_INSTALL_LIB_HEADERS}
        DESTINATION include/edlib)
install(FILES ${TARGET_INSTALL_EXT_HEADERS}
        DESTINATION include/ext)

option(Testing "Enable testing" OFF)
if(Testing)
    include(EnableGtests)
    #UseGtest()
    #include_directories(${Hubbard_SOURCE_DIR}/test)
    add_subdirectory(test)
    include_directories(${GTEST_INCLUDE_DIR})
    enable_testing()

    add_test(SzSymmetryTest test/SzSymmetryTest)
    add_test(NSymmetryTest test/NSymmetryTest)
    add_test(HubbardModelTest test/HubbardModelTest)
    add_test(LanczosTest test/LanczosTest)
    add_test(StaticObservablesTest test/StaticObservablesTest)

endif (Testing)


option(Examples "Enable examples" ON)
if(Examples)
    add_subdirectory(examples)
    set(SOURCE_FILES main.cpp)
    add_executable(Hubbard ${SOURCE_FILES})
    target_link_libraries(Hubbard common-lib ${PROJECT_NAME} ${extlibs})
    if(USE_MPI)
        target_link_libraries(Hubbard ${parlibs})
    endif(USE_MPI)
endif(Examples)
