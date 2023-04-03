#ifndef NEUTRINOS_Current_Library_Nucleon_Nucleon_H
#define NEUTRINOS_Current_Library_Nucleon_Nucleon_H

#include "NEUTRINOS++/Current_Library/Scatter_Current_Base.H"


namespace NEUTRINOS {
  class Nucleon_Nucleon : public Scatter_Current_Base {
    Complex m_cR, m_cL;
  public:
    Nucleon_Nucleon(const ATOOLS::Flavour_Vector& flavs,
		    const std::vector<int>& indices,
		    const std::string& name);
    ~Nucleon_Nucleon() {}
    void Calc(const ATOOLS::Vec4D_Vector& moms,METOOLS::XYZFunc * F);
  };
}
#endif
