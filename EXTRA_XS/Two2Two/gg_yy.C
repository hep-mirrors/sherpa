#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


namespace EXTRAXS {

  class gg_yy : public ME2_Base {
  private:
    double m_fac;
  public:

    gg_yy(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };

  gg_yy::gg_yy(const Process_Info& pi, const Flavour_Vector& fl)
    : ME2_Base(pi, fl)
  {
    rpa.gen.AddCitation
      (1,"The box diagram implementation for $gg \\to \\gamma \\gamma$ is\
 documented in \\cite{Hoeche:2009xc} and references therein.");
    m_oew=2;
    m_oqcd=2;
    m_sintt=1;
    double alphaqed2 = sqr(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa.gen.Ecms())));
    double alphas2 = sqr(MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_S",sqr(rpa.gen.Ecms())));

    for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;

    p_colours[0][0] = p_colours[1][1] = 500;
    p_colours[0][1] = p_colours[1][0] = 501;
    
    
    double qcharge[6] = {-1., 2, -1., 2., -1., 2.};
    double sumc = 0;
    for (size_t i=0; i<5; ++i) {
      sumc += sqr(qcharge[i]);
    }
    double sumc4 = sqr(sumc)/81.0;
    
    m_fac = 1.0;
    // sum over loop charges
    m_fac *= sumc4;

    // couplings at E_cms
    m_fac *= alphas2*alphaqed2;

    // average over initial state heli/colour
    m_fac /= 4.0*64;

    // sym factor for two photons
    m_fac /= 2.0;
  }
  
  double gg_yy::operator()(const ATOOLS::Vec4D_Vector& mom)
  {
    double s=(mom[0]+mom[1]).Abs2();
    double t=(mom[0]-mom[2]).Abs2();
    double u=(mom[0]-mom[3]).Abs2();

    double s2 = sqr(s);
    double t2 = sqr(t);
    double u2 = sqr(u);
    double pi2 = sqr(M_PI);

    double box = 
      256.0*sqr(-1.0-(t-u)/s*log(t/u)-(t2+u2)/s2*(sqr(log(t/u))+pi2)/2.0)+
      256.0*sqr(-1.0-(s-u)/t*log(-s/u)-(s2+u2)/t2*sqr(log(-s/u))/2.0)+
      256.0*pi2*sqr((s-u)/t+(s2+u2)/t2*log(-s/u))+
      256.0*sqr(-1.0-(s-t)/u*log(-s/t)-(s2+t2)/u2*sqr(log(-s/t))/2.0)+
      256.0*pi2*sqr((s-t)/u+(s2+t2)/u2*log(-s/t))+
      1280.0;

    // global factor see above
    return box*m_fac*CouplingFactor(2,2);
  }
}

DECLARE_TREEME2_GETTER(gg_yy_Getter,"gg_yy")
Tree_ME2_Base *gg_yy_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born)
    return NULL;
  if (pi.m_loopgenerator!="gg_yy") return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if (fl[0].Kfcode()==kf_gluon && fl[1].Kfcode()==kf_gluon &&
      fl[2].Kfcode()==kf_photon && fl[3].Kfcode()==kf_photon)
  {
    return new gg_yy(pi, fl);
  }
  return NULL;
}
