[![DOI](https://zenodo.org/badge/63707930.svg)](https://zenodo.org/badge/latestdoi/63707930)

##### Overview
EDLib is a C++ template finite temperature Exact diagonalization solver for quantum electron models.
The central class of the library is `Hamiltonian<Storage, Model>`, that is parametrized by Storage and Model:

- There exists a following set of implementation of models for common purposes:
    - `HubbardModel<precision>`. The finite Hubbard model cluster.
    - `SingleImpurityAndersonModel<precision>`. The single multi-orbital impurity Anderson Model.

- For the Hamiltonian matrix storage there are three implementation of sparse matrix storages:
    - `SpinResolvedStorage<Model>`. A storage that takes into account the case when hopping Hamiltonian 
    can be expressed as Kronecker sum for each spin. This storage is implemented with *MPI* support.
    - `SOCRSStorage<Model>`. A storage that store only fermion signs for each element in Hamiltonian. 
    This storage is implemented with *OpenMP* support.
    - `CRSStorage<Model>`. A simple CRS storage.

The resluting eigenpairs are stored as a set of `EigenPair<precision, SymmetrySectorType>` structures in 
the Hamiltonian object. 

The following observable can be computed by means of Lanczos continuous fraction 
(`Lanczos<Hamiltonian, Mesh, MeshArguments...>` class template) of the Lehmann representation:
- Single-particle Green's function (`GreensFunction<Hamiltonian, Mesh, MeshArguments...>` class template).
- Spin suseptibility (`ChiLoc<Hamiltonian, Mash, MeshArguments...>` class template).
Greens functions are implemented on top of *ALPSCore* Greens functions module and can use either positive 
Matsubara frequency mesh or Real frequency mesh.

Look for examples in the "examples/" directory for a detailed information.

##### Installation ###
The code is is provided as a header-only library with a set of examples and tests.
At least the `edlib/Hamiltonian.h` should be included in any derivative projects.

The eigensolver uses [cpp-arnoldi](https://github.com/Q-solvers/cpp-arnoldi), a header-only C++17 port
of the ARPACK symmetric driver. Clone it alongside EDLib or point CMake at its location:

```
git clone https://github.com/Q-solvers/cpp-arnoldi ../cpp-arnoldi-main
```

To compile examples and tests create a build directory and run 

1. `cmake -DExamples=ON -DTesting=ON {path_to_edlib}`
2. `make`
3. `make test` (for running tests)
4. example will be build in examples subdirectory

To build with MPI support add `-DUSE_MPI=ON` *CMake* flag. *MPI* library should be installed and *ALPSCore* 
library should be compiled with *MPI* support. To build with a specific *ALPSCore* library 
`-DALPSCore_DIR=<path to ALPSCore>` *CMake* flag.

##### Dependencies 
- c++17-compatible compiler (gcc >= 7, clang >= 5, icpc >= 19)  
- *ALPSCore* library >= 0.5.6-alpha3
- *cpp-arnoldi* (header-only, bundled or cloned alongside EDLib)
- *BLAS* and *LAPACK*
- *MPI* standard >= 2.1 (optional)
- *git* to fetch the code 
- *cmake* >= 3.8.2 to build tests and examples (optional)

##### Authors
- Sergei Iskakov, *iskakoff[at]q-solvers.ru*, 2016-now
- Michael Danilov, 2016-now

##### Distribution
Open-source under MIT License.
