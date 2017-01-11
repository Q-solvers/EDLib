//
// Created by iskakoff on 23/08/16.
//

#include <gtest/gtest.h>
#include "edlib/HubbardModel.h"

#include "edlib/SpinResolvedStorage.h"
#include "edlib/SingleImpurityAndersonModel.h"
#include "edlib/CRSStorage.h"
#include "edlib/EDParams.h"
#include <fstream>
#include <iomanip>


class SpinResolvedStorageTestEnv : public ::testing::Environment {
  protected:
  virtual void SetUp() {
    char** argv;
    int argc = 0;
    //int mpiError = MPI_Init(&argc, &argv);
  }

  virtual void TearDown() {
    //MPI_Finalize();
  }

  ~SpinResolvedStorageTestEnv(){};

};

::testing::Environment* const foo_env = AddGlobalTestEnvironment(new SpinResolvedStorageTestEnv);


TEST(SpinResolvedStorageTest, Init) {
  alps::params p;
  EDLib::define_parameters(p);
  p["NSITES"]=8;
  p["storage.MAX_SIZE"] = 80000;
  p["storage.MAX_DIM"] = 4900;
  EDLib::Model::HubbardModel<double> m(p);
  EDLib::Model::HubbardModel<double> m2(p);
  EDLib::Storage::SpinResolvedStorage<EDLib::Model::HubbardModel<double> > storage(p, m
#ifdef USE_MPI
  , MPI_COMM_WORLD
#endif
);
  EDLib::Storage::CRSStorage<EDLib::Model::HubbardModel<double> > storage2(p, m2
#ifdef USE_MPI
  , MPI_COMM_WORLD
#endif
);
  while (m.symmetry().next_sector() && m2.symmetry().next_sector()) {
    std::vector<double> v(m.symmetry().sector().size(), 1.0);
    std::vector<double> w(m.symmetry().sector().size(), 0.0);
    std::vector<double> w2(m.symmetry().sector().size(), 0.0);
    m.symmetry().init();
    m2.symmetry().init();
    storage.fill();
    storage2.fill();
    //storage.prepare_work_arrays(v.data());
    //storage.av(v.data(), w.data(), m.symmetry().sector().size());
    //storage2.av(v.data(), w2.data(), m.symmetry().sector().size());
    //storage.finalize();
    for(int i = 0; i< w.size(); ++i) {
      ASSERT_EQ(w[i], w2[i]);
    }
  }
}

TEST(SpinResolvedStorageTest, av) {
  alps::params p;
  EDLib::define_parameters(p);
  p["NSITES"]=8;
  p["storage.MAX_SIZE"] = 80000;
  p["storage.MAX_DIM"] = 4900;
  EDLib::Model::HubbardModel<double> m(p);
  EDLib::Storage::SpinResolvedStorage<EDLib::Model::HubbardModel<double> > storage(p, m
#ifdef USE_MPI
  , MPI_COMM_WORLD
#endif
  );
  typedef typename EDLib::Symmetry::SzSymmetry::Sector Stype;
  Stype s(4,4,4900);
  m.symmetry().set_sector(s);
  m.symmetry().init();
  storage.fill();
  size_t vs = storage.vector_size(s);
  std::cout<<"size:"<<vs<<std::endl;
  std::vector<double> v(storage.vector_size(s), 1.0);
  std::vector<double> w(storage.vector_size(s), 0.0);
  storage.prepare_work_arrays(v.data());
  storage.av(v.data(), w.data(), vs);
  storage.finalize();
  std::vector<double> vv(s.size());
  std::ifstream ifile("v.dat");
  for(int i = 0; i < s.size(); ++i)
    ifile>>vv[i];
  ifile.close();
  int k =0;
#ifdef USE_MPI
  int rank,size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout<<"RANK:"<<rank<<" of "<<size<<std::endl;
  for(int np = 0; np<size; ++np) {
  MPI_Barrier(MPI_COMM_WORLD);
  if(np==rank) {
#endif
  std::ofstream ofile("w.dat", std::ios_base::app);
  for(int i = 0; i<w.size();++i) {
    ofile<<w[i]<<"\n";
  }
  ofile.close();
#ifdef USE_MPI
  }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0) {
#endif
  std::vector<double> ww(s.size());
  std::ifstream ifilew("w.dat");
  std::cout<<s.size()<<"\n";
  for(int i = 0; i < s.size(); ++i)
    ifilew>>ww[i];
  ifilew.close();
  std::remove("w.dat");
  for(int i = 0; i<ww.size();++i) {
    std::cout<<i<<" "<<vv[i]<<" "<<ww[i]<<"\n";
    ASSERT_EQ(vv[i], ww[i]);
  }
#ifdef USE_MPI
 }
#endif
}

