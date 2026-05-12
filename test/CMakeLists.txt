
add_executable(SzSymmetryTest SzSymmetry_Test.cpp)
add_executable(SzSymmetryWithBosonTest SzSymmetryWithBoson_Test.cpp)
add_executable(NSymmetryTest NSymmetry_Test.cpp)
add_executable(HubbardModelTest HubbardModel_Test.cpp)
add_executable(SpinResolvedStorage SRS.cpp  SpinResolvedStorage_Test.cpp)
add_executable(StaticObservablesTest StaticObservables_Test.cpp)
add_executable(StatisticsTest Statistics_Test.cpp)
add_executable(LanczosTest Lanczos_Test.cpp)


target_link_libraries(SzSymmetryTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(SzSymmetryWithBosonTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(NSymmetryTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(HubbardModelTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(SpinResolvedStorage EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(StaticObservablesTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(StatisticsTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(LanczosTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})

if (USE_MPI)
    target_link_libraries(SzSymmetryTest ${parlibs})
    target_link_libraries(SzSymmetryWithBosonTest ${parlibs})
    target_link_libraries(NSymmetryTest ${parlibs})
    target_link_libraries(HubbardModelTest ${parlibs})
    target_link_libraries(SpinResolvedStorage ${parlibs})
    target_link_libraries(StaticObservablesTest ${parlibs})
    target_link_libraries(StatisticsTest ${parlibs})
    target_link_libraries(LanczosTest ${parlibs})
endif(USE_MPI)


file(COPY input DESTINATION ${CMAKE_BINARY_DIR}/test)
