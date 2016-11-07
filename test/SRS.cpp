
#include <gtest/gtest.h>

#ifdef USE_MPI
#include <mpi.h>
#endif


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
#ifdef USE_MPI
  MPI_Init(&argc, &argv);
#endif
  int res = RUN_ALL_TESTS();
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return res;
}
