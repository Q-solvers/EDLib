include_directories(${Hubbard_SOURCE_DIR}/include)

add_library(common-lib dummy.f90)

target_link_libraries(common-lib)