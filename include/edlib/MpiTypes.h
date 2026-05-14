#ifndef EDLIB_MPITYPES_H
#define EDLIB_MPITYPES_H

#ifdef USE_MPI
#include <mpi.h>
#include <complex>
#include <cstddef>

namespace edlib {

  /**
   * MPI datatype lookup for the small set of arithmetic types EDLib reduces
   * over. Replaces alps::mpi::detail::mpi_type<T>().
   */
  template <class T>
  MPI_Datatype mpi_type() = delete;

  template <> inline MPI_Datatype mpi_type<float>()  { return MPI_FLOAT; }
  template <> inline MPI_Datatype mpi_type<double>() { return MPI_DOUBLE; }
  template <> inline MPI_Datatype mpi_type<int>()    { return MPI_INT; }
  template <> inline MPI_Datatype mpi_type<long>()   { return MPI_LONG; }
  template <> inline MPI_Datatype mpi_type<long long>() { return MPI_LONG_LONG; }
  template <> inline MPI_Datatype mpi_type<unsigned>()  { return MPI_UNSIGNED; }
  template <> inline MPI_Datatype mpi_type<unsigned long>() { return MPI_UNSIGNED_LONG; }
  template <> inline MPI_Datatype mpi_type<unsigned long long>() { return MPI_UNSIGNED_LONG_LONG; }
  template <> inline MPI_Datatype mpi_type<std::complex<float>>()  { return MPI_C_FLOAT_COMPLEX; }
  template <> inline MPI_Datatype mpi_type<std::complex<double>>() { return MPI_C_DOUBLE_COMPLEX; }

}
#endif // USE_MPI

#endif
