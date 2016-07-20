include_directories(${Hubbard_SOURCE_DIR}/lib)

add_library(hubbard-lib Hamiltonian.cpp Hamiltonian.h Combination.cpp Combination.h SzCombination.cpp SzCombination.h Sector.cpp Sector.h)

target_link_libraries(hubbard-lib ${extlibs})