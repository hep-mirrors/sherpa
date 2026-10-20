#ifndef METOOLS_Current_Library_V_Dirac_Dirac_H
#define METOOLS_Current_Library_V_Dirac_Dirac_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {
  class V_Dirac_Dirac : public Current_Base {
    private:
    Complex m_cR, m_cL;
    FormFactor_Base* p_f1;
    
    public:
      void CheckTopology() const override;
    
      V_Dirac_Dirac(
          const ATOOLS::Flavour_Vector& flavs,
          const std::vector<int>& indices,
          const std::string& name);
    
      ~V_Dirac_Dirac();
    
      void SetModelParameters(struct GeneralModel model);
    
      void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    };
}

#endif


