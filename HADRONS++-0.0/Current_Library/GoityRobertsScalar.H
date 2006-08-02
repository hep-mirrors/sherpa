#ifndef GoityRobertsScalar_H
#define GoityRobertsScalar_H

#include "Current_Base.H"

namespace HADRONS {

  class GoityRobertsScalar : public Current_Base {

    double m_Vcb;
  public:
    GoityRobertsScalar(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}

#endif
