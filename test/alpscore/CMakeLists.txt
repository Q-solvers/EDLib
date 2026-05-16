
add_executable(SzSymmetryAlpsCoreTest SzSymmetry_Test.cpp)
add_executable(SzSymmetryWithBosonAlpsCoreTest SzSymmetryWithBoson_Test.cpp)
add_executable(NSymmetryAlpsCoreTest NSymmetry_Test.cpp)
add_executable(HubbardModelAlpsCoreTest HubbardModel_Test.cpp)
add_executable(SpinResolvedStorageAlpsCoreTest SRS.cpp SpinResolvedStorage_Test.cpp)
add_executable(StaticObservablesAlpsCoreTest StaticObservables_Test.cpp)
add_executable(StatisticsAlpsCoreTest Statistics_Test.cpp)
add_executable(LanczosAlpsCoreTest Lanczos_Test.cpp)


target_link_libraries(SzSymmetryAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(SzSymmetryWithBosonAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(NSymmetryAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(HubbardModelAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(SpinResolvedStorageAlpsCoreTest  EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(StaticObservablesAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(StatisticsAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})
target_link_libraries(LanczosAlpsCoreTest EDLib ${extlibs} ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY})

if (USE_MPI)
    target_link_libraries(SzSymmetryAlpsCoreTest ${parlibs})
    target_link_libraries(SzSymmetryWithBosonAlpsCoreTest ${parlibs})
    target_link_libraries(NSymmetryAlpsCoreTest ${parlibs})
    target_link_libraries(HubbardModelAlpsCoreTest ${parlibs})
    target_link_libraries(SpinResolvedStorageAlpsCoreTest  ${parlibs})
    target_link_libraries(StaticObservablesAlpsCoreTest ${parlibs})
    target_link_libraries(StatisticsAlpsCoreTest ${parlibs})
    target_link_libraries(LanczosAlpsCoreTest ${parlibs})
endif(USE_MPI)


file(COPY ../input DESTINATION ${CMAKE_BINARY_DIR}/test)
