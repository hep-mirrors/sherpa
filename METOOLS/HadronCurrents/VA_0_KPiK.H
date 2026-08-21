#ifndef METOOLS_Current_Library_VA_0_KPiK_H
#define METOOLS_Current_Library_VA_0_KPiK_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>


namespace METOOLS {
  // Generic current for the "KpiK" isospin family of
  // Finkemeier & Mirkes, hep-ph/9503474: K^-pi^-K^+, K^0pi^-K0bar,
  // K_Spi^-K_S, K_Spi^-K_L, K_Lpi^-K_L, K^-pi^0K^0. Structurally
  // identical to VA_0_PiPiPi (same V1=(p1-p3)_T, V2=(p2-p3)_T basis,
  // same F1*v1+F2*v2+F3*q+FS*v4 assembly), but for 3 fully distinct
  // pseudoscalars (no Bose symmetry) - kept as a separate class/tag
  // because the ordering convention (0=q1,1=q2,2=q3 in Finkemeier-
  // Mirkes' own per-channel labeling, fixed by FF_0_PPP_Base::FixMode())
  // is channel-specific and unrelated to VA_0_KKPi's own (K+,K-,pi-)
  // convention for the physically similar K+K-pi- final state.
  class VA_0_KPiK : public Current_Base {
  private:
    double  m_norm;

    FormFactor_Base * p_f1, * p_f2, * p_f3, * p_fS;
  public:
    VA_0_KPiK(const ATOOLS::Flavour_Vector& flavs,
	      const std::vector<int>& indices,
	      const std::string& name);
    ~VA_0_KPiK();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
