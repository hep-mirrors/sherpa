#include "EXTRA_XS/NLO/Virtual_ME2_Base.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Process/Process_Base.H"

using namespace EXTRAXS;
using namespace PHASIC;
using namespace ATOOLS;

#define CF 1.33333333333333333
#define CA 3.
#define TR 0.5

namespace EXTRAXS {
  class ee3jet_QCD_Virtual : public Virtual_ME2_Base {
    double m_cpl;
    int m_nf;
     PHASIC::Process_Base* p_tree;

    double Q2,y12,y13,y23;
    double Ff(double,double,double);
    double Rf(double,double);
    double Tf(double,double,double);
    double Kf(double,double,double);

  public:
    ee3jet_QCD_Virtual(const Process_Info& pi, const Flavour_Vector& flavs,
                       PHASIC::Process_Base* tree) :
      Virtual_ME2_Base(pi, flavs), p_tree(tree)
    {
      DEBUG_INFO("ee3jet virtual opened ...");
      m_cpl=MODEL::s_model->ScalarFunction(std::string("alpha_S"),
                                           sqr(rpa.gen.Ecms()));
      m_cpl/=2.*M_PI;

      Flavour hfl(kf_quark);
      m_nf = hfl.Size()/2;
    }

    ~ee3jet_QCD_Virtual() {
       if (p_tree) delete p_tree;
    }

    void Calc(const ATOOLS::Vec4D_Vector& mom);

    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom) {
      return 2.*M_PI*m_mur2/(mom[0]*mom[2]);
    }
  };
}


void ee3jet_QCD_Virtual::Calc(const Vec4D_Vector& mom) {
  m_born=p_tree->Differential(mom); // todo: divide out ISR for hadron-hadron
  
  m_born*=m_cpl;
  Q2  = mom[0]*mom[1];
  y12 = mom[3]*mom[4]/Q2;
  y13 = mom[3]*mom[2]/Q2;
  y23 = mom[4]*mom[2]/Q2;

  // finite
  double dren = 0.5*(11.-2./3.*(double)m_nf)*log(m_mur2/Q2);
  double ep = 0.5*((2.*CF-CA)*sqr(log(y12))+CA*(sqr(log(y13))+sqr(log(y23))));
  m_res.Finite()=(dren-ep+0.5*sqr(M_PI)*(2.*CF+CA)-8.*CF+Ff(y12,y13,y23))*m_born;

  // eps_IR^1
  double epIR = (2.*CF-CA)*log(y12)+CA*(log(y13)+log(y23));
  m_res.IR()=(epIR-3.*CF-11./6.*CA+2./3.*TR*m_nf)*m_born;

  // eps_IR^2
  m_res.IR2()=-(2.*CF+CA)*m_born;
}

double ee3jet_QCD_Virtual::Rf(double y1,double y2)
{
  return log(y1)*log(y2)-log(y1)*log(1.-y1)-log(y2)*log(1.-y2)+sqr(M_PI)/6.-DiLog(y1)-DiLog(y2);
}

double ee3jet_QCD_Virtual::Tf(double y1,double y2,double y3)
{
  return (sqr(y1+y2)+sqr(y1+y3))/(y2*y3);
}

double ee3jet_QCD_Virtual::Kf(double y1,double y2,double y3)
{
  return log(y2)*(CF*(y1*(4.*y1+2.*y2+4.*y3)+y2*y3)/sqr(y1+y3)+CA*y2/(y1+y3));
}

double ee3jet_QCD_Virtual::Ff(double y12,double y13,double y23)
{ //F/T in ref ERT
  double h = (sqr(y12)+sqr(y12+y13))/(y13*y23)*Rf(y12,y23);
  h += (sqr(y12)+sqr(y12+y23))/(y13*y23)*Rf(y12,y13);
  h += (sqr(y13)+sqr(y23))/(y13*y23*(y13+y23));
  h -= 2.*log(y12)*(1./sqr(y13+y23)-1.);
  double res = CF*(y12/(y12+y13)+y12/(y12+y23)+(y12+y23)/y13+(y12+y13)/y23);
  res += Kf(y12,y13,y23) + Kf(y12,y23,y13);
  res -= (2.*CF-CA)*h;
  res /= Tf(y12,y13,y23);
  res -= CA*Rf(y13,y23);
  return res;
}



DECLARE_VIRTUALME2_GETTER(ee3jet_QCD_Virtual_Getter,"ee3jet_QCD_Virtual")
Virtual_ME2_Base *ee3jet_QCD_Virtual_Getter::operator()(const Process_Info &pi) const
{
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo) return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl=pi.ExtractFlavours();
    if (fl.size()!=5) return NULL;
    if (fl[0]==Flavour(kf_e) && fl[1]==Flavour(kf_e).Bar() &&
        fl[2].IsGluon() && fl[3].IsQuark() && fl[4].IsQuark()) {
      if ((pi.m_oqcd==1 || pi.m_oqcd==99) && (pi.m_oew==2 || pi.m_oew==99)) {
        Process_Info tree_pi(pi);
        tree_pi.m_fi.m_nloqcdtype=nlo_type::lo;
        Process_Base* tree = pi.p_gens->InitializeProcess(tree_pi,true);
        if (!tree) {
          DEBUG_INFO("no corresponding born found");
          return NULL;
        }
        std::string scale=ToString(sqr(rpa.gen.Ecms()));
        tree->SetScale("VAR["+scale+"]");
        tree->SetKFactor("QCD",1,2);
        return new ee3jet_QCD_Virtual(pi, fl, tree);
      }
    }
  }
  return NULL;
}
