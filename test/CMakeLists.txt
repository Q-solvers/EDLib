
include_directories(${Hubbard_SOURCE_DIR}/src)
include_directories(${Hubbard_SOURCE_DIR}/lib)

add_executable(SzSymmetryTest SzSymmetry_Test.cpp)

target_link_libraries(SzSymmetryTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})