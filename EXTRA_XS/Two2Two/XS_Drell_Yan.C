#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Interaction_Models/Interaction_Model_Base.H"
#include "ATOOLS/Org/Message.H"

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

    XS_ee_ffbar(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    bool   SetColours(double,double,double);
  };
}

XS_ee_ffbar::XS_ee_ffbar(const Process_Info& pi, const Flavour_Vector& fl)
  : ME2_Base(pi, fl)
{
  DEBUG_INFO("now entered EXTRAXS::XS_ee_ffbar ...");
  m_sintt=1;
  m_oew=2;
  m_oqcd=0;
  MZ2    = sqr(ATOOLS::Flavour(kf_Z).Mass());
  GZ2    = sqr(ATOOLS::Flavour(kf_Z).Width());
 
  alpha  = MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()));
  sin2tw = MODEL::s_model->ScalarConstant(string("sin2_thetaW"));
  if (ATOOLS::Flavour(kf_Z).IsOn()) 
    kappa  = 1./(4.*sin2tw*(1.-sin2tw));
  else
    kappa  = 0.;

  mass     = fl[2].Mass();
  qe       = fl[0].Charge();
  qf       = fl[2].Charge();
  ae       = fl[0].IsoWeak();      
  af       = fl[2].IsoWeak();
  ve       = ae - 2.*qe*sin2tw;
  vf       = af - 2.*qf*sin2tw;
  colfac   = 1.;

  kswitch  = 0;  
  fac      = 2./(3.*M_PI);
  fin      = 2.*M_PI/9. - 7./(3.*M_PI) + 9./(3.*M_PI);

  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (fl[0].IsLepton() && fl[1].IsLepton()) {
    barred = fl[2].IsAnti();
    p_colours[2][barred] = p_colours[3][1-barred] = 500;
    colfac = 3.;
  }

  if (fl[0].IsQuark() && fl[1].IsQuark())  {
    barred = fl[0].IsAnti();
    p_colours[0][barred] = p_colours[1][1-barred] = 500;
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
    t=(momenta[0]-momenta[2]).Abs2();
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

bool XS_ee_ffbar::SetColours(double s,double t,double u) 
{ 
  /*
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] = 
    m_scale[PHASIC::stp::sis] = 
    s+dabs(p_momenta[0].Abs2())+dabs(p_momenta[1].Abs2());
  m_scale[PHASIC::stp::sfs] = s;
  */
  return 1; 
}


DECLARE_TREEME2_GETTER(DY_Getter,"DY")
Tree_ME2_Base *DY_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  if ((fl[2].IsLepton() && fl[3]==fl[2].Bar() && fl[0].IsQuark() && 
       fl[1]==fl[0].Bar()) ||   
      (fl[0].IsLepton() && fl[1]==fl[0].Bar() && fl[2].IsQuark() && 
       fl[3]==fl[2].Bar())) {
    if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
      return new XS_ee_ffbar(pi, fl);
    }
  }
  return NULL;
}


namespace EXTRAXS {

  class XS_Charged_Drell_Yan : public ME2_Base {  // == XS_ffbar_ee but not XS_ffbar_f'fbar' !
  private:

    bool   barred;
    double kappa,sin2tw,MW2,GW2,alpha;
    double term;
    double colfac;

  public:

    XS_Charged_Drell_Yan(const Process_Info& pi, const Flavour_Vector& fl);

    double operator()(const ATOOLS::Vec4D_Vector& mom);
    bool   SetColours(double,double,double);
  };
}

XS_Charged_Drell_Yan::XS_Charged_Drell_Yan(const Process_Info& pi, const Flavour_Vector& fl)
  : ME2_Base(pi, fl)
{
  DEBUG_INFO("now entered EXTRAXS::XS_Charged_Drell_Yan ...");

  m_sintt=1;
  m_oew=2;
  m_oqcd=0;

  MW2    = sqr(ATOOLS::Flavour(kf_Wplus).Mass());
  GW2    = sqr(ATOOLS::Flavour(kf_Wplus).Width());

  alpha  = MODEL::s_model->GetInteractionModel()->ScalarFunction("alpha_QED",sqr(rpa->gen.Ecms()));
  sin2tw = MODEL::s_model->ScalarConstant(string("sin2_thetaW"));
  kappa  = 1./(2.*sin2tw);

  for (short int i=0;i<4;i++) p_colours[i][0] = p_colours[i][1] = 0;
  if (fl[0].IsLepton() && fl[1].IsLepton()) {
    barred = fl[2].IsAnti();
    p_colours[2][barred] = p_colours[3][1-barred] = 500;
    colfac = 3.;
  }

  if (fl[0].IsQuark() && fl[1].IsQuark())  {
    barred = fl[0].IsAnti();
    p_colours[0][barred] = p_colours[1][1-barred] = 500;
    colfac  = 1./3.;
  }
  if (fl[2].IntCharge()+fl[3].IntCharge()>0) {
    m_cfls[12].push_back(kf_Wplus);
    m_cfls[3].push_back(Flavour(kf_Wplus).Bar());
  }
  else {
    m_cfls[3].push_back(kf_Wplus);
    m_cfls[12].push_back(Flavour(kf_Wplus).Bar());
  }
}

double XS_Charged_Drell_Yan::operator()(const ATOOLS::Vec4D_Vector& mom) {
  double s=(mom[0]+mom[1]).Abs2();
  double u=(mom[0]-mom[3]).Abs2();

  term = sqr(kappa)/((sqr(s-MW2)+MW2*GW2)) * (u*u);
  return sqr(4.*M_PI*alpha) * CouplingFactor(0,2) * colfac * term;
}

bool XS_Charged_Drell_Yan::SetColours(double s,double t,double u)
{
  /*
  m_scale[PHASIC::stp::fac] = m_scale[PHASIC::stp::ren] =
    m_scale[PHASIC::stp::sis] =
    s+dabs(p_momenta[0].Abs2())+dabs(p_momenta[1].Abs2());
  m_scale[PHASIC::stp::sfs] = s;
  */
  return 1;
}


DECLARE_TREEME2_GETTER(CDY_Getter,"CDY") // Charged Drell-Yan
Tree_ME2_Base *CDY_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_fi.NLOType()!=nlo_type::lo && pi.m_fi.NLOType()!=nlo_type::born) return NULL;
  Flavour_Vector fl=pi.ExtractFlavours();
  if (fl.size()!=4) return NULL;
  for (size_t i(0);i<fl.size();++i) if (fl[i].Mass()!=0.) return NULL;
  if (!ATOOLS::Flavour(kf_Wplus).IsOn()) return NULL;
  if ((fl[2].IsLepton() && fl[2]!=fl[3].Bar() &&
       fl[2].LeptonFamily()==fl[3].LeptonFamily() &&
       !fl[2].IsAnti() && fl[3].IsAnti() &&
       fl[0].IsQuark() && fl[0]!=fl[1].Bar() &&
       fl[0].QuarkFamily()==fl[1].QuarkFamily())) {
    if ((pi.m_oqcd==0 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
      return new XS_Charged_Drell_Yan(pi, fl);
    }
  }
  return NULL;
}
