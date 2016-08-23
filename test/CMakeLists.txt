
include_directories(${Hubbard_SOURCE_DIR}/src)
include_directories(${Hubbard_SOURCE_DIR}/lib)

add_executable(SzSymmetryTest SzSymmetry_Test.cpp)
add_executable(NSymmetryTest NSymmetry_Test.cpp)
add_executable(HubbardModelTest HubbardModel_Test.cpp)
add_executable(HolsteinAndersonModelTest HolsteinAndersonModel_Test.cpp)
add_executable(NSymmetryWithBosonTest NSymmetryWithBoson_Test.cpp)

target_link_libraries(SzSymmetryTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(NSymmetryTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(NSymmetryWithBosonTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(HubbardModelTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(HolsteinAndersonModelTest common-lib hubbard-lib ${extlibs} ${GTEST_LIBRARY})