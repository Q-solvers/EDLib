
project(edlib CXX)

add_custom_target(edlib SOURCES Combination.h
    CRSStorage.h
    EDParams.h
    EigenPair.h
    FermionicModel.h
    fortranbinding.h
    GreensFunction.h
    ChiLoc.h
    Hamiltonian.h
    HubbardModel.h
    Lanczos.h
    NSymmetry.h
    SingleImpurityAndersonModel.h
    SOCRSStorage.h
    SpinResolvedStorage.h
    StateDescription.h
    Storage.h
    Symmetry.h
    SzSymmetry.h
    HDF5Utils.h
    MeshFactory.h
    CommonUtils.h
    CRSMatrix.h
    ExecutionStatistic.h)