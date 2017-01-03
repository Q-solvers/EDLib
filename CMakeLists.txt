cmake_minimum_required(VERSION 3.0)
project(Hubbard C CXX Fortran)

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include(CXXStandard)
use_cxx11()

find_package(ALPSCore COMPONENTS params hdf5 gf REQUIRED)
find_package(ARPACK REQUIRED)
find_package(BLAS)
find_package(LAPACK)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALPSCore_INCLUDES})


include(EnableGtests)
UseGtest()
include_directories(${GTEST_INCLUDE_DIR})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DBOOST_DISABLE_ASSERTS -DNDEBUG -g -O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -m64")


set(extlibs
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${ALPSCore_LIBRARIES}
    ${PARPACK_LIB}
    ${ARPACK_LIB}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    )


#include_directories(${Hubbard_SOURCE_DIR}/test)
#add_subdirectory(test)

if(USE_MPI)
    if(!ALPS_HAVE_MPI)
        MESSAGE(FATAL_ERROR "Enable MPI support in ALPSCore.")
    endif(!ALPS_HAVE_MPI)
    find_package(MPI REQUIRED)
    if(MPI_FOUND)
        SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI")
    endif(MPI_FOUND)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAG}")
endif(USE_MPI)

include_directories(${Hubbard_SOURCE_DIR}/test)
add_subdirectory(test)

include_directories(${Hubbard_SOURCE_DIR}/include)
include_directories(${Hubbard_SOURCE_DIR}/lib)

add_subdirectory(lib)
add_subdirectory(include/edlib)
add_subdirectory(examples)

enable_testing()

#add_test(SzSymmetryTest test/SzSymmetryTest)
add_test(HubbardModelTest test/HubbardModelTest)

set(SOURCE_FILES main.cpp)

add_executable(Hubbard ${SOURCE_FILES})

target_link_libraries(Hubbard common-lib ${extlibs})
