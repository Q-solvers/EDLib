
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
