
include_directories(${Hubbard_SOURCE_DIR}/include)
include_directories(${Hubbard_SOURCE_DIR}/lib)

add_executable(SzSymmetryTest SzSymmetry_Test.cpp)
add_executable(NSymmetryTest NSymmetry_Test.cpp)
add_executable(HubbardModelTest HubbardModel_Test.cpp)
add_executable(SpinResolvedStorage SRS.cpp  SpinResolvedStorage_Test.cpp)

target_link_libraries(SzSymmetryTest common-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(NSymmetryTest common-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(HubbardModelTest common-lib ${extlibs} ${GTEST_LIBRARY})
target_link_libraries(SpinResolvedStorage common-lib ${extlibs} ${GTEST_LIBRARY})
