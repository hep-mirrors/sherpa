#ifndef METOOLS_Current_Library_VA_0_EtaPiPi_H
#define METOOLS_Current_Library_VA_0_EtaPiPi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/Tools.H"
#include <list>


namespace METOOLS {
  // eta pi^- pi^0 / eta' pi^- pi^0 (tau_two_meson_currents_KS_RChiT.tex,
  // Sec."eta pi pi and eta' pi pi"). Same generic v1/v2/v4 Calc()
  // structure as VA_0_KPiK/VA_0_KPiPi - see VA_0_KPiK.H for the
  // general design note - but F1/F2/F3 are always identically zero
  // here (G-parity forbids axial/scalar contributions in the isospin
  // limit): only FS (the anomalous vector current) is ever nonzero.
  // Order required: 0=eta(or eta'), 1=pi^-, 2=pi^0 - NOT
  // interchangeable (unlike VA_0_KPiK/VA_0_KPiPi's Bose-symmetric-ish
  // channels), since eta(') is not Bose-symmetric with either pion.
  class VA_0_EtaPiPi : public Current_Base {
  private:
    double  m_norm;

    FormFactor_Base * p_f1, * p_f2, * p_f3, * p_fS;
  public:
    VA_0_EtaPiPi(const ATOOLS::Flavour_Vector& flavs,
		 const std::vector<int>& indices,
		 const std::string& name);
    ~VA_0_EtaPiPi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
