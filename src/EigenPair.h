//
// Created by iskakoff on 27/07/16.
//

#ifndef HUBBARD_EIGENPAIR_H
#define HUBBARD_EIGENPAIR_H

#include <memory>
#include <vector>

namespace EDLib {

  template<typename precision, class SectorType>
  class EigenPair {
  public:
    EigenPair(const precision &eval, const std::vector < precision > &evec, SectorType sec) : _eigenvalue(eval), _sector(sec),
                                                                                              _eigenvector(new precision[evec.size()], std::default_delete < precision[] >()) {
      ;
      std::cout<<"evec:"<<evec.size()<<std::endl;
      std::memcpy(_eigenvector.get(), &evec[0], evec.size() * sizeof(precision));
    };

    virtual ~EigenPair() {

    }

    precision eigenvalue() const {
      return _eigenvalue;
    }

    const std::shared_ptr < precision > &eigenvector() const {
      return _eigenvector;
    }

    const SectorType &sector() const {
      return _sector;
    }

    bool operator>(const EigenPair &pair) const {
      return _eigenvalue > pair._eigenvalue;
    };

    bool operator<(const EigenPair &pair) const {
      return _eigenvalue < pair._eigenvalue;
    };
  private:
    precision _eigenvalue;
    std::shared_ptr < precision > _eigenvector;
    SectorType _sector;
  };

}

#endif //HUBBARD_EIGENPAIR_H
