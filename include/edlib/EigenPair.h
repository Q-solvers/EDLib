#ifndef EDLIB_EIGENPAIR_H
#define EDLIB_EIGENPAIR_H

#include <vector>

namespace edlib {

  template <class Precision, class SectorType>
  class EigenPair {
  public:
    EigenPair(const Precision& eval,
              const std::vector<Precision>& evec,
              int id,
              SectorType sec)
        : _eigenvalue(eval), _sector(sec), _eigenvector(evec), _id(id) {}

    virtual ~EigenPair() = default;

    Precision eigenvalue() const { return _eigenvalue; }
    const std::vector<Precision>& eigenvector() const { return _eigenvector; }
    const SectorType& sector() const { return _sector; }

    bool operator<(const EigenPair& o) const {
      return (_eigenvalue < o._eigenvalue) || (_id < o._id);
    }
    bool operator>(const EigenPair& o) const {
      return (_eigenvalue > o._eigenvalue) || (_id > o._id);
    }

  private:
    Precision              _eigenvalue;
    std::vector<Precision> _eigenvector;
    int                    _id;
    SectorType             _sector;
  };

}

#endif
