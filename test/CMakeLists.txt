add_executable(HubbardTest    Hubbard_Test.cpp)
add_executable(LanczosTest    Lanczos_Test.cpp)
add_executable(SzSymmetryTest SzSymmetry_Test.cpp)
add_executable(NSymmetryTest  NSymmetry_Test.cpp)
if(USE_MPI)
    add_executable(MpiTest        Mpi_Test.cpp)
endif()
set(core_test_libs
    EDLib
    arnoldi::arnoldi
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    ${GTEST_LIBRARY}
    ${GTEST_MAIN_LIBRARY})

target_link_libraries(HubbardTest    ${core_test_libs})
target_link_libraries(LanczosTest    ${core_test_libs})
target_link_libraries(SzSymmetryTest ${core_test_libs})
target_link_libraries(NSymmetryTest  ${core_test_libs})

if(USE_MPI)
    target_link_libraries(HubbardTest    ${parlibs})
    target_link_libraries(LanczosTest    ${parlibs})
    target_link_libraries(SzSymmetryTest ${parlibs})
    target_link_libraries(NSymmetryTest  ${parlibs})
    target_link_libraries(MpiTest        ${parlibs} ${core_test_libs})
endif()

# LanczosTest reads test/input/GF_Chi/*.dat from CWD; copy it in.
file(COPY ${PROJECT_SOURCE_DIR}/test/input
     DESTINATION ${CMAKE_BINARY_DIR}/test)
