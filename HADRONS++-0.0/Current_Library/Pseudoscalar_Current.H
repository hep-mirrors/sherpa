#ifndef Pseudoscalar_Current_H
#define Pseudoscalar_Current_H

#include "Current_Base.H"


namespace HADRONS {
  class Pseudoscalar_Current : public Current_Base {
    double m_fxx, m_Vxx;
    int m_pionmode;
  public:
    Pseudoscalar_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc();
  };
}
// indices for Pseudoscalar_Current: 0 must be B, 1 HeavyScalar.
#endif
