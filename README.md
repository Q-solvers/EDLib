##### Overview
EDLib is a C++ template Exact diagonalization solver for quantum electron models. It provides three basic classes of objects:

- Hamiltonian<precision, Storage, Model>. 
There exists a set of implementation of models for common purposes:
    - `HubbardModel<precision>`. The finite Hubbard model cluster.
    - `SingleImpurityAndersonModel<precision>`. The single multi-orbital impurity Anderson Model.
For the Hamiltonian matrix storage there are three implementation of sparse matrix storages:
- `SpinResolvedStorage<precision, Model>`. A storage that takes into account the case when hopping Hamiltonian can be expressed as Kronecker sum for each spin. This storage is implemented with *MPI* support.
- `SOCRSStorage<precision, Model>`. A storage that store only fermion signs for each element in Hamiltonian. This storage is implemented with *OpenMP* support.
- `CRSStorage<precision, Model>`. A simple CRS storage.

The following observable can be computed by means of Lanczos continuous fraction of the Lehmann representation:
- Single-particle Green's function (`GreensFunction<precision, Hamiltonian>` class template).
- Spin suseptibility (`ChiLoc<precision, Hamiltonian>` class template).
Greens functions are implemented on top of *ALPSCore* Greens functions module.

Look for examples in "examples/" directory for detailed information.

##### Installation ###
The code is is provided as a header-only library with a set of examples and tests.
At least the `edlib/Hamiltonian.h` should be included in any derivative projects.
To compile examples and tests create a build directory and run 

1. `cmake -DExamples=ON -DTesting=ON {path_to_edlib}`
2. `make`
3. `make test` (for running tests)
4. example will be build in examples subdirectory

##### Dependencies 
- c++11-compatible compiler (tested with clang >= 3.1, gcc >= 4.8.2, icpc >= 14.0.2)  
- *ALPSCore* library >= 0.5.6-alpha3
- *arpack-ng* >= 3.5.0
- *MPI* standard >= 2.1 (optional)
- *git* to fetch the code 
- *cmake* to build tests and examples (optional)

##### Author
- Sergei Iskakov, *iskakoff[at]q-solvers.ru*, 2016-now

##### Distribution
Open-source under MIT License.
