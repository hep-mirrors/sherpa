#ifndef METOOLS_Current_Library_VA_0_KPiPi_H
#define METOOLS_Current_Library_VA_0_KPiPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>


namespace METOOLS {
  // Generic current for the "KpiPi" isospin family of
  // Finkemeier & Mirkes, hep-ph/9503474: K^-pi^-pi^+, pi^-K0bar pi^0.
  // See VA_0_KPiK.H for the general design note (identical structure,
  // different flavour content / getter tag).
  class VA_0_KPiPi : public Current_Base {
  private:
    double  m_norm;

    FormFactor_Base * p_f1, * p_f2, * p_f3, * p_fS;
  public:
    VA_0_KPiPi(const ATOOLS::Flavour_Vector& flavs,
	       const std::vector<int>& indices,
	       const std::string& name);
    ~VA_0_KPiPi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
