include_directories(${Hubbard_SOURCE_DIR}/lib)

add_library(hubbard-lib Hamiltonian.cpp Hamiltonian.h Symmetry.cpp Symmetry.h SzSymmetry.cpp SzSymmetry.h
    EigenPair.h HubbardModel.h Lanczos.h GreensFunction.h NSymmetry.h NSymmetryWithBoson.h Combination.h HolsteinAndersonModel.h SingleImpurityAndersonModel.h)

target_link_libraries(hubbard-lib ${extlibs})