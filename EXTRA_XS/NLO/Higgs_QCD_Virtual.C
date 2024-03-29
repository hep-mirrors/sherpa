#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#define CF 1.33333333333333333
#define CA 3.
#define TR 0.5

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class Higgs_QCD_Virtual : public Virtual_ME2_Base {
    double m_fac, m_pij, m_b0, m_finiteconst;
  public:
    Higgs_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs, int con) :
      Virtual_ME2_Base(pi, flavs)
    {
      Flavour lq(kf_quark);
      double nlf = double(lq.Size())/2.0;
      m_b0 = (11.0/3.0*CA-4.0/3.0*TR*nlf)/4.0/M_PI;
      if (flavs[0].IsQuark()) {
        m_fac = CF;
        m_pij = 1.5;
        /// @todo finite part
        THROW(not_implemented, "qq -> h virtual not implemented.");
      }
      else if (flavs[1].IsGluon()) {
        m_fac = CA;
        m_pij = 2.0*M_PI*m_b0/CA;
        m_finiteconst = sqr(M_PI) + (con?0.0:11.0/3.0);
      }
      else THROW(fatal_error, "Internal Error.");
    }

    ~Higgs_QCD_Virtual() {
    }

    void Calc(const ATOOLS::Vec4D_Vector& momenta);

  };
}

using namespace EXTRAXS;

void Higgs_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  double p2  = 2.*momenta[0]*momenta[1];
  // 1/epsIR
  m_res.IR()=(-2.*m_pij+2.*log(p2/m_mur2))*m_fac;
  // 1/epsIR2
  m_res.IR2()=-2.*m_fac;
  // finite
  m_res.Finite()=(m_finiteconst - sqr(log(p2/m_mur2)))*m_fac;
}

DECLARE_VIRTUALME2_GETTER(EXTRAXS::Higgs_QCD_Virtual,"Higgs_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<PHASIC::Virtual_ME2_Base,PHASIC::Process_Info,EXTRAXS::Higgs_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nlotype==nlo_type::loop) {
    if (pi.m_fi.m_nlocpl[1]!=0.) return NULL;
    if (pi.m_fi.m_asscontribs!=asscontrib::none) {
      msg_Error()<<"Higgs_QCD_Virtual(): Error: cannot provide requested "
                 <<"associated contributions "<<pi.m_fi.m_asscontribs<<std::endl;
      return NULL;
    }
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl[0].IsGluon() && fl[1].IsGluon() &&
        pi.m_fi.m_ps.size()==1 && pi.m_fi.m_ps[0].m_fl.Kfcode()==kf_h0) {
      if (pi.m_maxcpl[0]==3 && (pi.m_maxcpl[1]==1 || (pi.m_maxcpl.size()>2 && pi.m_maxcpl[2]==1)) &&
          pi.m_mincpl[0]==3 && (pi.m_mincpl[1]==1 || (pi.m_mincpl.size()>2 && pi.m_mincpl[2]==1))) {
        for (size_t i=2; i<fl.size(); ++i) {
          if (fl[i].Strong()) return NULL;
        }
        Settings& s = Settings::GetMainSettings();
        int con = s["HNNLO_KF_MODE"].Get<int>();
        return new Higgs_QCD_Virtual(pi, fl, con);
      }
    }
  }
  return NULL;
}
