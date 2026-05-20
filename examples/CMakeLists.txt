
if(EDLIB_WITH_ALPSCORE)
    add_executable(hubbard-example Hubbard.cpp)
    add_executable(anderson-example Anderson.cpp)
    add_executable(holstein-example HolsteinAnderson.cpp)

    target_link_libraries(hubbard-example EDLib ${extlibs})
    target_link_libraries(anderson-example EDLib ${extlibs})
    target_link_libraries(holstein-example EDLib ${extlibs})

    if(USE_MPI)
        target_link_libraries(hubbard-example ${parlibs})
        target_link_libraries(anderson-example ${parlibs})
        target_link_libraries(holstein-example ${parlibs})
    endif(USE_MPI)
endif()

# Core-only example: builds with or without ALPSCore.
add_executable(hubbard-core-example HubbardCore.cpp)
target_link_libraries(hubbard-core-example EDLib ${Boost_LIBRARIES} arnoldi::arnoldi
                                            ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
if(USE_MPI)
    target_link_libraries(hubbard-core-example ${parlibs})
endif(USE_MPI)

# CUDA storage example: built by nvcc, links cuBLAS/cudart via arnoldi.
if(EDLIB_USE_CUDA)
    add_executable(hubbard-cuda-example HubbardCuda.cu)
    target_link_libraries(hubbard-cuda-example EDLib ${Boost_LIBRARIES} arnoldi::arnoldi
                                               ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    set_target_properties(hubbard-cuda-example PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    add_executable(anderson-cuda-bench AndersonCudaBench.cu)
    target_link_libraries(anderson-cuda-bench EDLib ${Boost_LIBRARIES} arnoldi::arnoldi
                                              ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    set_target_properties(anderson-cuda-bench PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
endif()
