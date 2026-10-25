#include "YFS/Main/YFS_Weight_Builder.H"
#include "YFS/Main/YFS_Handler.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include <sstream>
#include <iomanip>

using namespace std;
using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;
using namespace PHASIC;
using namespace METOOLS;

void YFS_Weight_Builder::ProbeSign(const double w_full)
{
  static const bool sp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["SIGN_PROBE"].SetDefault(false).Get<bool>());
  if (sp) {
    double xi(0.), xf(0.);
    for (const Vec4D &k : m_h.m_ev.m_ISRPhotons) xi = Max(xi, 2.*k[0]/sqrt(m_h.m_s));
    for (const Vec4D &k : m_h.m_ev.m_FSRPhotons) xf = Max(xf, 2.*k[0]/sqrt(m_h.m_s));
    const double mll(m_h.m_ev.m_plab.size() > 3 ? (m_h.m_ev.m_plab[2] + m_h.m_ev.m_plab[3]).Mass() : -1.);
    const double weex(!IsZero(m_h.m_ev.m_corr_nominal) ? w_full*m_h.m_ev.m_corr_eex/m_h.m_ev.m_corr_nominal : 0.);
    std::ostringstream o;   // one write per line: MPI ranks share the stream
    o<<std::setprecision(10)<<"@@@ SIGN nisr="<<m_h.m_ev.m_ISRPhotons.size()
     <<" nfsr="<<m_h.m_ev.m_FSRPhotons.size()<<" xisr="<<xi<<" xfsr="<<xf
     <<" mll="<<mll<<" wnom="<<w_full<<" weex="<<weex<<"\n";
    std::cerr<<o.str();
  }
}
