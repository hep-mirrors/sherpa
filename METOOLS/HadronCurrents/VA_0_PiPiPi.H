#ifndef METOOLS_Current_Library_VA_0_PiPi_H
#define METOOLS_Current_Library_VA_0_PiPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>


namespace METOOLS {  
  class VA_0_PiPiPi : public Current_Base {    
  private:
    double  m_norm;

    FormFactor_Base * p_f1, * p_f2, * p_f3, * p_fS;
  public:
    VA_0_PiPiPi(const ATOOLS::Flavour_Vector& flavs,
		const std::vector<int>& indices,
		const std::string& name);
    ~VA_0_PiPiPi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
