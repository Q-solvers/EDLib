cmake_minimum_required(VERSION 3.3)
project(Hubbard)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(ALPSCore COMPONENTS mc params REQUIRED)
find_package(ARPACK)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${ALPSCore_INCLUDES})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11")


set(extlibs
    ${MPI_CXX_LIBRARIES}
    ${Boost_LIBRARIES}
    ${ALPSCore_LIBRARIES}
    )

include_directories(${Hubbard_SOURCE_DIR}/src)
include_directories(${Hubbard_SOURCE_DIR}/lib)
include_directories(${Hubbard_SOURCE_DIR}/test)

add_subdirectory(src)
add_subdirectory(lib)
add_subdirectory(test)

enable_testing()


set(SOURCE_FILES main.cpp)

add_executable(Hubbard ${SOURCE_FILES})

target_link_libraries(Hubbard hubbard-lib ${extlibs})