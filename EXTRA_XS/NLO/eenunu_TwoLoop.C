#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace MODEL;

namespace EXTRAXS {
class eenunu_TwoLoop : public PHASIC::Virtual_ME2_Base {
  double m_eps2, m_eps, m_fin;
  double ME, ME2, m_alpha;
  double m_u, m_t, m_s, m_photon_mass, m_irscale;

public:
  eenunu_TwoLoop(const Process_Info &pi, const Flavour_Vector &flavs,
                 const double &ep2, const double &ep);

  ~eenunu_TwoLoop() {}

  void Calc(const ATOOLS::Vec4D_Vector &momenta);
};
} // namespace EXTRAXS

using namespace EXTRAXS;

eenunu_TwoLoop::eenunu_TwoLoop(const Process_Info &pi,
                               const Flavour_Vector &flavs, const double &ep2,
                               const double &ep)
    : Virtual_ME2_Base(pi, flavs), m_eps2(ep2), m_eps(ep) {
  Scoped_Settings s{Settings::GetMainSettings()["YFS"]};
  m_photon_mass = s["PHOTON_MASS"].Get<double>();
  ME = Flavour(kf_e).Mass();
  ME2 = ME * ME;
  m_mode = 0;
  m_irscale = 100;
}

void eenunu_TwoLoop::Calc(const Vec4D_Vector &momenta) {
  double factor(1.);
  if (m_stype & sbt::qcd)
    factor = 2. * M_PI / AlphaQCD();
  else if (m_stype & sbt::qed)
    factor = 2. * M_PI / AlphaQED();
  else
    THROW(fatal_error, "Unknown coupling.");
  m_s = (momenta[2] + momenta[3]).Abs2();
  m_t = (momenta[0] - momenta[2]).Abs2();
  m_u = (momenta[0] - momenta[3]).Abs2();
  m_alpha = (*aqed)(m_s);
  DivArrD massph(0, -1, 0, 0, 0, 0);
  DivArrD massph2(0, 0, -1, 0, 0, 0);
  const double L = log(m_s / ME2);
  DivArrD res;
  m_res.Finite() = 1;
  res = 2.*(massph - log(4. * M_PI * sqr(m_irscale) / 4. / ME2 / 4. / M_PI)) * (L-1);
        // (-1. / 8. * L * L * L + 0.5 * L * L +
        //  (-7. / 8. + 5. / 2. * M_PI * M_PI / 6.) * L + 0.5 -
        //  13. / 4. * M_PI * M_PI / 6);
  res += 2.*(massph2 - 2. * log(4. * M_PI * sqr(m_irscale) / 4. / ME2 / 4. / M_PI)) *(L-1)*(L-1);
//                      (L * L / 8. - 0.25 * L + 1. / 8 - 3/4 * M_PI * M_PI / 6.);
//   res*=pow(m_alpha/M_PI, 2);
  m_res=res;
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::eenunu_TwoLoop, "eenunu_TwoLoop")
Virtual_ME2_Base *
ATOOLS::Getter<PHASIC::Virtual_ME2_Base, PHASIC::Process_Info,
               EXTRAXS::eenunu_TwoLoop>::operator()(const Process_Info &pi)
    const {
  if (pi.m_vvgenerator.find("Internal") != 0)
    return NULL;
  if (!pi.Has(nlo_type::vv))
    return NULL;
  Flavour_Vector fl(pi.ExtractFlavours());
  if (fl.size() != 4)
    return NULL;
  if (fl[2].Kfcode() == 12 || fl[2].Kfcode() == -12)
    return NULL;
  if ((fl[0].Kfcode() == Flavour(kf_e) && fl[1] == fl[0].Bar()) ||
      (fl[1].Kfcode() == Flavour(kf_e) && fl[1] == fl[0].Bar()) &&
          (fl[2].IsNeutrino() && fl[3].IsNeutrino()) &&
          (fl[2] == fl[3].Bar())) {
    return new eenunu_TwoLoop(pi, fl, 0.0, 0);
  }
  return NULL;
}
