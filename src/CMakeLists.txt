include_directories(${Hubbard_SOURCE_DIR}/lib)

add_library(hubbard-lib Hamiltonian.cpp Hamiltonian.h Simmetry.cpp Simmetry.h SzSimmetry.cpp SzSimmetry.h)

target_link_libraries(hubbard-lib ${extlibs})