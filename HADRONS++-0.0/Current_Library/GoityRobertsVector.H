#ifndef GoityRobertsVector_H
#define GoityRobertsVector_H

#include "Current_Base.H"

namespace HADRONS {

  class GoityRobertsVector : public Current_Base {

    double m_Vcb;
  public:
    GoityRobertsVector(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}

#endif
