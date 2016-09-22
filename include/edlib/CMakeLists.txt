
project(edlib CXX)

add_custom_target(edlib SOURCES Combination.h
    CRSStorage.h
    EDParams.h
    EigenPair.h
    FermionicModel.h
    fortranbinding.h
    GreensFunction.h
    Hamiltonian.h
    HolsteinAndersonModel.h
    HubbardModel.h
    Lanczos.h
    NSymmetry.h
    NSymmetryWithBoson.h
    SingleImpurityAndersonModel.h
    SOCRSStorage.h
    SpinResolvedStorage.h
    Storage.h
    Symmetry.h
    SzSymmetry.h)