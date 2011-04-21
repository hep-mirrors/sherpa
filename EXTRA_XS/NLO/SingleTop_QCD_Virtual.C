#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"

#define CF 1.33333333333333333
#define CA 3.
#define TR 0.5

using namespace PHASIC;
using namespace ATOOLS;

namespace EXTRAXS {
  class SingleTop_QCD_Virtual : public PHASIC::Virtual_ME2_Base {
    double m_fac;
  public:
    SingleTop_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs) :
      Virtual_ME2_Base(pi, flavs)
    {
      m_fac = CF;
    }

    ~SingleTop_QCD_Virtual() {
    }

    void   Calc(const ATOOLS::Vec4D_Vector& momenta);
    double DoublePole();
    double SinglePole(const ATOOLS::Vec4D_Vector& momenta);
    double Finite(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
  };
}

using namespace EXTRAXS;

void SingleTop_QCD_Virtual::Calc(const Vec4D_Vector& momenta) {
  // 1/epsIR
  m_res.IR()=SinglePole(momenta);
  // 1/epsIR2
  m_res.IR2()=DoublePole();
  // finite
  m_res.Finite()=Finite(momenta);
}

double SingleTop_QCD_Virtual::DoublePole() {
  return -3.*m_fac;
}

double SingleTop_QCD_Virtual::SinglePole(const Vec4D_Vector& momenta) {
  double s03 = (momenta[1]+momenta[3]).Abs2();
  double s2t = (momenta[0]+momenta[2]).Abs2();
  double mt  = Flavour(kf_t).Mass();
  return (2.*log(s03/m_mur2) - 3. + 2.*log((s2t-mt*mt)/(mt*sqrt(m_mur2))) - 1. -3./2.)*m_fac;
}

double SingleTop_QCD_Virtual::Finite(const Vec4D_Vector& momenta) {
  double s13 = (momenta[1]+momenta[3]).Abs2();
  double s2t = (momenta[0]+momenta[2]).Abs2();
  double s12 = (momenta[0]+momenta[1]).Abs2();
  double s23 = (momenta[0]+momenta[3]).Abs2();
  double s34 = (momenta[2]+momenta[3]).Abs2();
  double mt  = Flavour(kf_t).Mass();
  double mt2 = mt*mt;
  double mur = sqrt(m_mur2);
  double v13 = -0.5*sqr(log(s13/m_mur2)) +1.5*log(s13/m_mur2) -4. +M_PI*M_PI/12.;
  double v2t = DiLog(1.-mt2/(s2t-mt2)) -2. - M_PI*M_PI/24. 
             - 0.5*sqr(log((s2t-mt2)/mt/mur))+1./8.*sqr(log(mt2/m_mur2))
             +(s2t-mt2)/(4.*(2.*mt2-s2t))*log(mt2/m_mur2) 
             + log((s2t-mt2)/(mt*mur))*(1.-(s2t-mt2)/(2.*(2.*mt2-s2t)) 
             -0.5*log(mt2/m_mur2));
  double dmt = -2. + 1.5*log(mt2/m_mur2);
  double mw = Flavour(kf_Wplus).Mass();
  double nonfac = mt2*s12*s23*log((s2t-mt2)/mt2)/((2.*mt2-s2t)*sqr(s13+mw*mw));
  double born = s12*(s34-mt2)/sqr(s13+mw*mw);

  return 2.*(v13+v2t) + dmt + nonfac/born;
}

double SingleTop_QCD_Virtual::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

DECLARE_VIRTUALME2_GETTER(SingleTop_QCD_Virtual_Getter,"SingleTop_QCD_Virtual")
Virtual_ME2_Base *SingleTop_QCD_Virtual_Getter::operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=4) return NULL;
    if (fl[0].IsQuark() && fl[1].IsQuark() &&
	fl[2].IsQuark() && fl[3].IsQuark() &&
	(((fl[2] == Flavour(kf_t)) && (fl[3].Bar() != Flavour(kf_t))) || 
         ((fl[3] == Flavour(kf_t)) && (fl[2].Bar() != Flavour(kf_t)))) ) {
      if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
//         return new SingleTop_QCD_Virtual(pi, fl);
      }
    }
  }
  return NULL;
}
