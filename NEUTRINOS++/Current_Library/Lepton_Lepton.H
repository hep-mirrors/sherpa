#ifndef NEUTRINOS_Current_Library_Lepton_Lepton_H
#define NEUTRINOS_Current_Library_Lepton_Lepton_H

#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"


namespace NEUTRINOS {
  class Lepton_Lepton : public Scatter_Current_Base {
    bool    m_anti;
    Complex QED_coupling, QED_cR, QED_cL;
    Complex Weak_coupling, Weak_cR, Weak_cL;
  public:
    Lepton_Lepton(const ATOOLS::Flavour_Vector& flavs,
		  const std::vector<int>& indices,
		  const std::string& name);
    ~Lepton_Lepton() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F);
  };
}
#endif
