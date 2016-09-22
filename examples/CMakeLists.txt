include_directories(${Hubbard_SOURCE_DIR}/include)

add_executable(hubbard-example Hubbard.cpp)
add_executable(anderson-example Anderson.cpp)

target_link_libraries(hubbard-example common-lib ${extlibs})

target_link_libraries(anderson-example common-lib ${extlibs})