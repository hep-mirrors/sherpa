#include "Griffin_Virtual.H"

#include "AddOns/Griffin/Griffin_Interface.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace ATOOLS;

Griffin::Griffin_Virtual::Griffin_Virtual(const Process_Info& pi,
        const Flavour_Vector& flavs, const GriffinProcess& proc,
        const bool higherorder):
  Virtual_ME2_Base(pi, flavs), m_process(proc), m_higherorder(higherorder)
{
  m_mode=1;
  /*
     GRIFFIN's helicity sum X of   DiffXSecTerm is dsigma/dcos(theta) = N_c s X/(32 pi) (GRIFFIN
    examples/test_ap2.cc, N_c = 3 for d dbar), and Sherpa's squared ME has
    dsigma/dcos(theta) = |M|^2/(32 pi s), so |M|^2 = N_c s^2 X. 
    The default normalisation is 1/pi, so the |M|^2 that Griffin_Virtual
  */
  const bool fix(ATOOLS::Settings::GetMainSettings()["GRIFFIN"]["FIX_NORMALISATION"]
                 .SetDefault(false).Get<bool>());
  const double nc(flavs.size() > 2 && flavs[2].Strong() ? 3. : 1.);
  m_norm = fix ? nc : 1./M_PI;
}

void Griffin::Griffin_Virtual::Calc(const Vec4D_Vector& momenta) {
  m_lastmom = momenta;
  if (m_higherorder)
    Griffin_Interface::Instance().EvaluateHigherOrder(momenta, m_process, m_res);
  else
    Griffin_Interface::Instance().EvaluateLoop(momenta, m_process, m_res);
  const double s = (momenta[0]+momenta[1]).Abs2();
  m_res.Finite() *= s*s*m_norm;
}



bool Griffin::Griffin_Virtual::HelicityAmplitudes(std::vector<Complex> &a0,
                                                  std::vector<Complex> &a1,
                                                  std::vector<std::vector<int> > &hel)
{
  if (!m_higherorder || m_lastmom.size() != 4) return false;
  return Griffin_Interface::Instance().HigherOrderHelicity(m_lastmom, m_process,
                                                           a0, a1, hel);
}


DECLARE_VIRTUALME2_GETTER(Griffin::Griffin_Virtual,"Griffin_Virtual")
Virtual_ME2_Base *ATOOLS::Getter<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,Griffin::Griffin_Virtual>::
operator()(const PHASIC::Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  const bool ho(pi.m_loopgenerator=="GriffinHO");
  if (pi.m_loopgenerator!="Griffin" && !ho) return nullptr;
  if (pi.m_fi.m_nlotype!=nlo_type::loop) return nullptr;
  Flavour_Vector flavs = pi.ExtractFlavours();
  if (flavs.size()!=4) return nullptr;
  return new Griffin::Griffin_Virtual(pi, flavs,
                                       Griffin::Griffin_Interface::MakeProcess(pi), ho);
}
