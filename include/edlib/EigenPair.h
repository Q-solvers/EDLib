//
// Created by iskakoff on 27/07/16.
//

#ifndef HUBBARD_EIGENPAIR_H
#define HUBBARD_EIGENPAIR_H

#include <vector>

namespace EDLib {
  template<typename precision, class SectorType>
  class EigenPair {
  public:

    EigenPair(const precision &eval, const std::vector < precision > &evec, int id, SectorType sec) : _eigenvalue(eval), _sector(sec),
                                                                                              _eigenvector(evec), _id(id) {
    };

    virtual ~EigenPair() {

    }

    precision eigenvalue() const {
      return _eigenvalue;
    }

    const std::vector < precision > &eigenvector() const {
      return _eigenvector;
    }

    const SectorType &sector() const {
      return _sector;
    }

    bool operator>(const EigenPair &pair) const {
      return (_eigenvalue > pair._eigenvalue) || (_id > pair._id);
    };

    bool operator<(const EigenPair &pair) const {
      return (_eigenvalue < pair._eigenvalue) || (_id < pair._id);
    };
  private:
    precision _eigenvalue;
    std::vector < precision > _eigenvector;
    int _id;
    SectorType _sector;
  };

}

#endif //HUBBARD_EIGENPAIR_H
