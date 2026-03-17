#ifndef METOOLS_Current_Library_VA_0_PiPi_H
#define METOOLS_Current_Library_VA_0_PiPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>


namespace METOOLS {  
  class VA_0_PiPi : public Current_Base {    
  private:
    double  m_norm, m_deltaM2;

    FormFactor_Base * p_fplus, * p_fzero;
  public:
    VA_0_PiPi(const ATOOLS::Flavour_Vector& flavs,
	      const std::vector<int>& indices,
	      const std::string& name);
    ~VA_0_PiPi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
