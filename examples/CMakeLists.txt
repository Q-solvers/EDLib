
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