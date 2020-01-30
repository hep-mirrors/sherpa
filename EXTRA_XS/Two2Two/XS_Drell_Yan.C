#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace std;


/* 
   In all the differential cross sections the factor 1/16 Pi is cancelled
   by the factor 4 Pi for each alpha. Hence one Pi remains in the game.
*/

namespace EXTRAXS {

  class XS_ee_ffbar : public ME2_Base {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
  private:

    bool   barred;
    double qe,qf,ae,af,ve,vf,mass;
    double kappa,sin2tw,MZ2,GZ2,alpha;
    double chi1,chi2,term1,term2;
    double colfac;
    int    kswitch;
    double fac,fin;

  public:

    XS_ee_ffbar(const External_ME_Args& args);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
  };
}

XS_ee_ffbar::XS_ee_ffbar(const External_ME_Args& args)
  : ME2_Base(args)
{
  DEBUG_INFO("now entered EXTRAXS::XS_ee_ffbar ...");
  m_sintt=1;
  m_oew=2;
  m_oqcd=0;
  MZ2    = sqr(ATOOLS::Flavour(kf_Z).Mass());
  GZ2    = sqr(ATOOLS::Flavour(kf_Z).Width());
 
  alpha  = MODEL::s_model->ScalarConstant("alpha_QED");
  sin2tw = std::abs(MODEL::s_model->ComplexConstant(string("csin2_thetaW")));
  if (ATOOLS::Flavour(kf_Z).IsOn()) 
    kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  else
    kappa  = 0.;

  mass     = m_flavs[2].Mass();
  qe       = m_flavs[0].Charge();
  qf       = m_flavs[2].Charge();
  ae       = m_flavs[0].IsoWeak();      
  af       = m_flavs[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  colfac   = 1.;

  kswitch  = 0;  
  fac      = 2./(3.*M_PI);
  fin      = 2.*M_PI/9. - 7./(3.*M_PI) + 9./(3.*M_PI);

  for (short int i=0;i<4;i++) m_colours[i][0] = m_colours[i][1] = 0;
  if (m_flavs[0].IsLepton() && m_flavs[1].IsLepton()) {
    barred = m_flavs[2].IsAnti();
    m_colours[2][barred] = m_colours[3][1-barred] = 500;
    colfac = 3.;
  }

  if (m_flavs[0].IsQuark() && m_flavs[1].IsQuark())  {
    barred = m_flavs[0].IsAnti();
    m_colours[0][barred] = m_colours[1][1-barred] = 500;
    colfac  = 1./3.;
    kswitch = 1;
  }

  m_cfls[3].push_back(kf_photon);
  m_cfls[3].push_back(kf_Z);
  m_cfls[12].push_back(kf_photon);
  m_cfls[12].push_back(kf_Z);
}

double XS_ee_ffbar::operator()(const ATOOLS::Vec4D_Vector& momenta) {
  double s(0.),t(0.);
  if (kswitch == 0 || kswitch==1) {
    s=(momenta[0]+momenta[1]).Abs2();
    t=(momenta[barred]-momenta[2]).Abs2();
  }
  else if (kswitch==2) { // meaning of t and s interchanged in DIS
    t=(momenta[0]+momenta[1]).Abs2();
    s=(momenta[0]-momenta[2]).Abs2();
  }
  else THROW(fatal_error,"Internal error.")

  //if (s<m_threshold) return 0.;
  chi1  = kappa * s * (s-MZ2)/(sqr(s-MZ2) + GZ2*MZ2);
  chi2  = sqr(kappa * s)/(sqr(s-MZ2) + GZ2*MZ2);

  term1 = (1+sqr(1.+2.*t/s)) * (sqr(qf*qe) + 2.*(qf*qe*vf*ve) * chi1 +
				(ae*ae+ve*ve) * (af*af+vf*vf) * chi2);
  term2 = (1.+2.*t/s) * (4. * qe*qf*ae*af * chi1 + 8. * ae*ve*af*vf * chi2);

  // Divide by two ????
  return sqr(4.*M_PI*alpha) * CouplingFactor(0,2) * colfac * (term1+term2); 
}

DECLARE_TREEME2_GETTER(XS_ee_ffbar,"XS_ee_ffbar")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,XS_ee_ffbar>::
operator()(const External_ME_Args& args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=4) return NULL;
  if ((fl[2].IsLepton() && fl[3]==fl[2].Bar() && fl[0].IsQuark() && 
       fl[1]==fl[0].Bar()) ||   
      (fl[0].IsLepton() && fl[1]==fl[0].Bar() && fl[2].IsQuark() && 
       fl[3]==fl[2].Bar())) {
    if (args.m_orders[0]==0 && args.m_orders[1]==2) {
      return new XS_ee_ffbar(args);
    }
  }
  return NULL;
}
