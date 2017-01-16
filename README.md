[![DOI](https://zenodo.org/badge/63707930.svg)](https://zenodo.org/badge/latestdoi/63707930)

##### Overview
EDLib is a C++ template finite temperature Exact diagonalization solver for quantum electron models.
The central class of the library is `Hamiltonian<Storage, Model>`, that is parametrized by Storage and Model:

- There exists a following set of implementation of models for common purposes:
    - `HubbardModel<precision>`. The finite Hubbard model cluster.
    - `SingleImpurityAndersonModel<precision>`. The single multi-orbital impurity Anderson Model.

- For the Hamiltonian matrix storage there are three implementation of sparse matrix storages:
    - `SpinResolvedStorage<Model>`. A storage that takes into account the case when hopping Hamiltonian can be expressed as Kronecker sum for each spin. This storage is implemented with *MPI* support.
    - `SOCRSStorage<Model>`. A storage that store only fermion signs for each element in Hamiltonian. This storage is implemented with *OpenMP* support.
    - `CRSStorage<Model>`. A simple CRS storage.

The resluting eigenpairs are stored as a set of `EigenPair<precision, SymmetrySectorType>` structures in th hHamiltonian object. 

The following observable can be computed by means of Lanczos continuous fraction (`Lanczos<Hamiltonian>` class template) of the Lehmann representation:
- Single-particle Green's function (`GreensFunction<Hamiltonian>` class template).
- Spin suseptibility (`ChiLoc<Hamiltonian>` class template).
Greens functions are implemented on top of *ALPSCore* Greens functions module.

Look for examples in "examples/" directory for detailed information.

##### Installation ###
The code is is provided as a header-only library with a set of examples and tests.
At least the `edlib/Hamiltonian.h` should be included in any derivative projects.
To compile examples and tests create a build directory and run 

1. `cmake -DARPACK_DIR=<path to ARPACK-ng library dir> -DExamples=ON -DTesting=ON {path_to_edlib}`
2. `make`
3. `make test` (for running tests)
4. example will be build in examples subdirectory

To build with MPI support add `-DUSE_MPI=ON` *CMake* flag. *MPI* library should be installed and *ALPSCore* library
should be compiled with *MPI* support. To build with a specific *ALPSCore* library `-DALPSCore_DIR=<path to ALPSCore>` *CMake* flag.
Since the critical for current library implementation MPI-related *ARPACK-ng* bug was recenlty fixed it is stricly recommended 
to use the latest version of *ARPACK-ng* from github repository.

##### Dependencies 
- c++11-compatible compiler (tested with clang >= 3.1, gcc >= 4.8.2, icpc >= 14.0.2)  
- *ALPSCore* library >= 0.5.6-alpha3
- *arpack-ng* >= 3.5.0
- *MPI* standard >= 2.1 (optional)
- *git* to fetch the code 
- *cmake* to build tests and examples (optional)

##### Authors
- Sergei Iskakov, *iskakoff[at]q-solvers.ru*, 2016-now
- Michael Danilov, 2016-now

##### Distribution
Open-source under MIT License.
