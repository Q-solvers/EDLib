cmake_minimum_required(VERSION 3.0)
project(Hubbard C CXX Fortran)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(ALPSCore COMPONENTS params hdf5 gf REQUIRED)
find_package(ARPACK)
find_package(Eigen3)
find_package(BLAS)
find_package(LAPACK)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALPSCore_INCLUDES})
include_directories(${EIGEN3_INCLUDE_DIR})


include(EnableGtests)
UseGtest()
include_directories(${GTEST_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS --std=c++11 -DNDEBUG -mtune=native -O3 -ffast-math")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")


set(extlibs
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ALPSCore_LIBRARIES}
    ${ARPACK_LIB}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    )

include_directories(${Hubbard_SOURCE_DIR}/src)
include_directories(${Hubbard_SOURCE_DIR}/lib)
include_directories(${Hubbard_SOURCE_DIR}/test)

add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(test)

enable_testing()

add_test(SzSymmetryTest test/SzSymmetryTest)

set(SOURCE_FILES main.cpp)

add_executable(Hubbard ${SOURCE_FILES})

target_link_libraries(Hubbard common-lib hubbard-lib ${extlibs})