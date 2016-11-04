
#include <gtest/gtest.h>
#include <mpi.h>


int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  MPI_Init(&argc, &argv);
  int res = RUN_ALL_TESTS();
  MPI_Finalize();
  return res;
}
