//
// Created by iskakoff on 23/08/16.
//

#ifndef HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
#define HUBBARD_SINGLEIMPURITYANDERSONMODEL_H

template<typename precision>
class SingleImpurityAndersonModel {
public:
  SingleImpurityAndersonModel(EDParams &p) : _symmetry(p) {}

  SzSymmetry& symmetry() {
    return _symmetry;
  }

private:
  SzSymmetry _symmetry;
};


#endif //HUBBARD_SINGLEIMPURITYANDERSONMODEL_H
