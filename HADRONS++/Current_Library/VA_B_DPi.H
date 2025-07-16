#ifndef HADRONS_Current_Library_VA_B_DPi_H
#define HADRONS_Current_Library_VA_B_DPi_H

#include "HADRONS++/Current_Library/Current_Base.H"

namespace HADRONS {

  class VA_B_DPi : public Current_Base {

    double m_Vxx;
  public:
    VA_B_DPi(const ATOOLS::Flavour_Vector& flavs,
             const std::vector<int>& indices, const std::string& name) :
      Current_Base(flavs, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
  };
}

#endif
