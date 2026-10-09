/*!
  \file ee_aa_QED_VV.C

  VV_Generator: Internal for e+e- -> gamma gamma, see ee_aa_QED_VV.H.
*/
#include "EXTRA_XS/NNLO/ee_aa_QED_VV.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_TwoLoop.H"
#include "EXTRA_XS/NNLO/ee_aa_QED_Poles.H"
#include "PHASIC++/Process/Process_Info.H"
#include "ATOOLS/Math/MathTools.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

namespace {

  //! e- e+ (either order) -> gamma gamma
  bool IsEEToAA(const Flavour_Vector &fl)
  {
    return fl.size() == 4 && fl[0].Kfcode() == kf_e && fl[1] == fl[0].Bar()
      && fl[2].IsPhoton() && fl[3].IsPhoton();
  }

}

ee_aa_QED_VV::ee_aa_QED_VV(const Process_Info &pi, const Flavour_Vector &flavs):
  Virtual_ME2_Base(pi, flavs)
{
  // Catani-subtracted: finite by construction
  m_providespoles = false;
}

void ee_aa_QED_VV::Calc(const Vec4D_Vector &p)
{
  namespace AGT = ee_aa_QED_TwoLoop;
  m_res = METOOLS::DivArrD(0.);
  m_born = 0.;
  m_ok = false;
  if (p.size() < 4) return;
  const double s((p[0] + p[1]).Abs2()), t((p[0] - p[2]).Abs2());
  if (!(s > 0.)) return;
  // symmetric under t <-> u, so the photon labelling does not matter
  m_x = Max(1e-12, Min(1. - 1e-12, -t/s));
  m_s = s;
  // the pieces at mu^2 = s (S = 0), for scheme conversions downstream
  m_fin02 = AGT::Fin02(m_x);
  m_fin11 = AGT::Fin11(m_x);
  m_r1 = AGT::R1(m_x);
  // the provider's own result at the renormalisation scale (S = ln s/mu^2)
  const double S(m_mur2 > 0. ? std::log(s/m_mur2) : 0.);
  m_res.Finite() = AGT::Fin02(m_x, S) + AGT::Fin11(m_x, S);
  // spin-averaged tree, 2 e^4 (t/u + u/t) = e^4 T0/4
  const double e2(4.*M_PI*AlphaQED());
  m_born = e2*e2*AGT::T0(m_x)/4.;
  m_ok = true;
}

void ee_aa_QED_VV::MasslessSeries(double n02[5], double n11[5], double m1[5]) const
{
  ee_aa_QED_Poles::MasslessSeries(m_x, n02, n11, m1);
}

void ee_aa_QED_VV::MassiveSeries(double me, double m1m[3], double m2m[3]) const
{
  ee_aa_QED_Poles::MassiveSeries(m_x, std::log(m_s/(me*me)), m1m, m2m, true);
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::ee_aa_QED_VV, "ee_aa_QED_VV")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base, PHASIC::Process_Info, EXTRAXS::ee_aa_QED_VV>::
operator()(const Process_Info &pi) const
{
  if (pi.m_vvgenerator != "Internal") return NULL;
  if (pi.m_fi.m_nlotype != nlo_type::vv) return NULL;
  const Flavour_Vector fl(pi.ExtractFlavours());
  if (!IsEEToAA(fl)) return NULL;
  return new EXTRAXS::ee_aa_QED_VV(pi, fl);
}
