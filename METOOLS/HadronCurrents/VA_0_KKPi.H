#ifndef METOOLS_Current_Library_VA_0_KKPi_H
#define METOOLS_Current_Library_VA_0_KKPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>

///////////////////////////////////////////////////////////////////////////
//
// Skeleton for tau^- -> K^+ K^- pi^- nu_tau (and charge conjugate),
// cf. arXiv:1509.09140, Sec.3.3 / Sec.3.1 for the general 3-meson
// current, Eq.(3.1).
//
// KKpi convention (as discussed - happy to change, this is just what
// we picked to get started):
//   Order: 0 = K^+, 1 = K^-, 2 = pi^-
// i.e. p1 = K^+, p2 = K^-, p3 = pi^- in the notation of Eq.(3.1), so
//   s1 = (p2+p3)^2 = M^2(K^- pi^-),  s2 = (p1+p3)^2 = M^2(K^+ pi^-)
//
// Naming carried over unchanged from VA_0_PiPiPi for consistency:
//   p_f1, p_f2: the two independent axial-vector form factors
//               (F1, F2 of Eq.(3.1), in the same transverse-projected
//               two-vector basis VA_0_PiPiPi already uses).
//   p_f3:       pseudoscalar form factor, multiplies q^mu directly
//               (this is F4 of Eq.(3.1); suppressed ~m_pi^2/q^2).
//   p_fS:       vector / odd-intrinsic-parity form factor, multiplies
//               the Levi-Civita contraction of p1,p2,p3 (this is F5 of
//               Eq.(3.1)). For pi pi pi this vanishes by G-parity
//               (F5^3pi=0, cf. paper Sec.3.1) - for KKpi it does NOT,
//               which is the one real structural difference from
//               VA_0_PiPiPi's Calc().
//
// STATUS: this is a structural skeleton only. Calc() implements the
// full Lorentz decomposition of Eq.(3.1); the individual form factors
// (CPC / RChL, cf. Sec.3.3) still need to be implemented in the
// FF_0_PPP getter (FF_0_PPP.C) under the "F1_0_KKPi" etc. names used
// in SetModelParameters() below. Currently those return 0.
//
///////////////////////////////////////////////////////////////////////////

namespace METOOLS {
  class VA_0_KKPi : public Current_Base {
  private:
    double  m_norm;

    FormFactor_Base * p_f1, * p_f2, * p_f3, * p_fS;
  public:
    VA_0_KKPi(const ATOOLS::Flavour_Vector& flavs,
	      const std::vector<int>& indices,
	      const std::string& name);
    ~VA_0_KKPi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
