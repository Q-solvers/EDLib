#ifndef EDLIB_EIGENPAIR_H
#define EDLIB_EIGENPAIR_H

#include <vector>

namespace edlib {

  // Evec defaults to a host std::vector<Precision> so existing
  // EigenPair<prec,Sector> uses are unchanged. Device-resident storages
  // instantiate it with their device buffer type; only eigenvalue()/id are
  // used for set ordering, so the eigenvector type never affects ordering.
  template <class Precision, class SectorType,
            class Evec = std::vector<Precision>>
  class EigenPair {
  public:
    EigenPair(const Precision& eval,
              const Evec& evec,
              int id,
              SectorType sec)
        : _eigenvalue(eval), _sector(sec), _eigenvector(evec), _id(id) {}

    virtual ~EigenPair() = default;

    Precision eigenvalue() const { return _eigenvalue; }
    const Evec& eigenvector() const { return _eigenvector; }
    const SectorType& sector() const { return _sector; }

    bool operator<(const EigenPair& o) const {
      return (_eigenvalue < o._eigenvalue) || (_id < o._id);
    }
    bool operator>(const EigenPair& o) const {
      return (_eigenvalue > o._eigenvalue) || (_id > o._id);
    }

  private:
    Precision  _eigenvalue;
    Evec       _eigenvector;
    int        _id;
    SectorType _sector;
  };

}

#endif
