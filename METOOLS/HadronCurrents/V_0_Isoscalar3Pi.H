#ifndef METOOLS_HadronCurrents_V_0_Isoscalar3Pi_H
#define METOOLS_HadronCurrents_V_0_Isoscalar3Pi_H

#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/FF_0_Isoscalar3Pi.H"

namespace METOOLS {

  ///////////////////////////////////////////////////////////////////////////
  //
  //  gamma* -> pi+ pi- pi0, the isoscalar electromagnetic current.
  //
  //  The sibling of VA_0_PiPiPi for a photon rather than a W.  A photon has no
  //  axial coupling to three pions, so of the four Lorentz structures there
  //  only the anomalous one survives and the current is
  //
  //      H^mu = eps^{mu al be ga} p_+al p_-be p_0ga  F(q^2; s_ij)
  //
  //  i.e. VA_0_PiPiPi with F1 = F2 = F3 = 0 and FS the isoscalar form factor.
  //  It is a separate class rather than a form-factor option of that one
  //  because the epsilon structure is not written out here at all: it comes
  //  from the declared rho pi chains inside FF_0_Isoscalar3Pi.
  //
  //  Being a Current_Base means the contraction against the leptonic side is
  //  the generic helicity sum, the same one Current_ME performs for the tau,
  //  rather than a tensor written out by hand in the matrix element.
  //
  ///////////////////////////////////////////////////////////////////////////

  class V_0_Isoscalar3Pi : public Current_Base {
  private:
    GeneralModel        m_model;
    FF_0_Isoscalar3Pi * p_ff;
  public:
    V_0_Isoscalar3Pi(const ATOOLS::Flavour_Vector& flavs,
		     const std::vector<int>& indices,
		     const std::string& name);
    ~V_0_Isoscalar3Pi();
    void Calc(const ATOOLS::Vec4D_Vector& moms, bool anti);
    void SetModelParameters(struct GeneralModel model);
  };
}

#endif
