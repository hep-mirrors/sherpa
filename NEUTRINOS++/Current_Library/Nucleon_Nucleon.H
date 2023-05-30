#ifndef NEUTRINOS_Current_Library_Nucleon_Nucleon_H
#define NEUTRINOS_Current_Library_Nucleon_Nucleon_H

#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"


namespace NEUTRINOS {
  class Nucleon_Nucleon : public Scatter_Current_Base {
    double  m_massin, m_massout;
    Complex QED_coupling, QED_cR, QED_cL;
    Complex Weak_NC_coupling, Weak_NC_cR, Weak_NC_cL;
    Complex Weak_CC_coupling, Weak_CC_cR, Weak_CC_cL;
  public:
    Nucleon_Nucleon(const ATOOLS::Flavour_Vector& flavs,
		    const std::vector<int>& indices,
		    const std::string& name);
    ~Nucleon_Nucleon() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F);
  };
}
#endif
