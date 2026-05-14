add_executable(HubbardCoreTest    HubbardCore_Test.cpp)
add_executable(LanczosCoreTest    LanczosCore_Test.cpp)
add_executable(SzSymmetryCoreTest SzSymmetryCore_Test.cpp)
add_executable(NSymmetryCoreTest  NSymmetryCore_Test.cpp)

set(core_test_libs
    EDLib
    arnoldi::arnoldi
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${GTEST_LIBRARY}
    ${GTEST_MAIN_LIBRARY})

target_link_libraries(HubbardCoreTest    ${core_test_libs})
target_link_libraries(LanczosCoreTest    ${core_test_libs})
target_link_libraries(SzSymmetryCoreTest ${core_test_libs})
target_link_libraries(NSymmetryCoreTest  ${core_test_libs})

if(USE_MPI)
    target_link_libraries(HubbardCoreTest    ${parlibs})
    target_link_libraries(LanczosCoreTest    ${parlibs})
    target_link_libraries(SzSymmetryCoreTest ${parlibs})
    target_link_libraries(NSymmetryCoreTest  ${parlibs})
endif()

# LanczosCoreTest reads test/input/GF_Chi/*.dat from CWD; copy it in.
file(COPY ${PROJECT_SOURCE_DIR}/test/input
     DESTINATION ${CMAKE_BINARY_DIR}/test)
