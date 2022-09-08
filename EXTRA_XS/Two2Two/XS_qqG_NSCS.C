#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "EXTRA_XS/Two2Two/NSCS_Tools.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace NSCS;

namespace EXTRAXS {

  class XS_qqG_NSCS : public ME2_Base {
  private:

    double qe, qf;
    double alpha, alpha_S, colfac;

  public:

    XS_qqG_NSCS(const External_ME_Args& args): ME2_Base(args)
    {
      DEBUG_INFO("now entered EXTRAXS::XS_qqG_NSCS ...");
      m_sintt=1;
      m_oew=2;
      m_oqcd=2;
      alpha = MODEL::s_model->ScalarConstant("alpha_QED");
      alpha_S = MODEL::s_model->ScalarConstant("alpha_S");
      qe = m_flavs[0].Charge();
      qf = m_flavs[3].Charge();
      colfac=3.*4./3.;
      m_cfls[3].push_back(kf_photon);
      m_cfls[3].push_back(kf_Z);
      m_cfls[60].push_back(kf_photon);
      m_cfls[60].push_back(kf_Z);
    }

#ifndef CHECK_LIMIT
    double operator()(const Vec4D_Vector &p)
    {
      NSCS_Args C(p[3],p[4],p[2],Vec4D(1.e-12,0.,0.,1.e-12));
      ClusterNNLO(C);
#else
    double operator()(const Vec4D_Vector &ip)
    {
      static double x1(-1.), x3;
      static int sec(2), mode(0), pts(1000);
      if (x1<0.) {
	Settings& s = Settings::GetMainSettings();
	Scoped_Settings nscs{ s["NSCS"] };
	nscs["x1"].SetDefault(1.);
	nscs["x3"].SetDefault(1.);
	nscs["pts"].SetDefault(10000000);
	nscs["mode"].SetDefault(0);
	nscs["sec"].SetDefault(0);
	x1=nscs["x1"].Get<double>();
	x3=nscs["x3"].Get<double>();
	pts=nscs["pts"].Get<int>();
	mode=nscs["mode"].Get<int>();
	sec=nscs["sec"].Get<int>();
	DEBUG_VAR(x1<<" "<<x3<<" "<<mode<<" "<<sec);
      }
      NSCS_Args C(ip[3],ip[4],ip[2],Vec4D(1.e-12,0.,0.,1.e-12));
      C.m_sec=sec;
      ClusterNNLO(C);
      if (x1<1.) C.m_x1=x1;
      if (x3<1.) C.m_x3=x3;
      ConstructNNLO(C);
      if (IsBad(C.m_p[1][0])) return 0.;
      Vec4D_Vector p(ip);
      p[3]=C.m_p[0];
      p[4]=C.m_p[1];
      p[2]=C.m_p[2];
      DEBUG_VAR(p[3]);
      DEBUG_VAR(p[4]);
      DEBUG_VAR(p[2]);
#endif
      double s12(2.*p[3]*p[4]), s13(2.*p[3]*p[2]), s23(2.*p[4]*p[2]);
      double Q2(s12+s13+s23);
      double R=s23/s13+s13/s23+2.*Q2*s12/(s13*s23);
      DEBUG_VAR(R);
      // double y132=s13/(s12+s13+s23), y231=s23/(s12+s13+s23);
      // if (y132<1.e-8 || y231<1.e-8) return 0.;
      // double z1=s12/(s12+s23), z2=s12/(s12+s13);
      // double D132 = 1./y132 * (2./(1.-z1*(1.-y132))-(1.+z1));
      // double D231 = 1./y231 * (2./(1.-z2*(1.-y231))-(1.+z2));
      // double S = D132+D231;
      double S=Q2*2.*s12/(s13*s23);
      DEBUG_VAR(S);
      {
	C.m_x3=0.;
	ConstructNNLO(C);
	double x(C.m_p[0][0]/(C.m_p[0][0]+C.m_p[2][0]));
	double C1=Q2/s13*(1.-x);
	double SC1=Q2/s13*2.*x/(1.-x);
	DEBUG_VAR(x/(s12/(s12+s23))-1.);
	DEBUG_VAR(SC1+C1<<" "<<SC1);
	S+=C1;
      }
      {
	C.m_x3=0.;
	ConstructNNLO(C);
	double x(C.m_p[1][0]/(C.m_p[1][0]+C.m_p[2][0]));
	double C2=Q2/s23*(1.-x);
	double SC2=Q2/s23*2.*x/(1.-x);
	DEBUG_VAR(x/(s12/(s12+s13))-1.);
	DEBUG_VAR(SC2+C2<<" "<<SC2);
	S+=C2;
      }
      DEBUG_VAR(R<<" "<<S<<" "<<1.-S/R);
      double cpl=sqr(4.*M_PI*alpha)*(8.*M_PI*alpha_S)
	*CouplingFactor(1,2)*colfac*sqr(qe)/Q2;
#ifdef CHECK_LIMIT
      static double sumr(0.), sums(0.), sumrs(0.), sumrsonr(0.), n(0.);
      double psf(1.0);
      if (mode) {
	psf=x1;
	if (mode&2) psf*=x1*x3;
      }
      sumr+=psf*cpl*(qf*qf*R);
      // sums+=psf*cpl*(qf*qf*(D132+D231));
      // sumrs+=psf*cpl*(qf*qf*(RS-D132-D231));
      sums+=psf*cpl*(qf*qf*(S));
      sumrs+=psf*cpl*(qf*qf*dabs(R-S));
      sumrsonr+=dabs(1.-S/R);
      n+=1.;
      if ((int(n)%pts==0)) {
	msg_Info()<<sumr/n<<" "<<sums/n<<" "<<sumrs/n<<" "<<sumrsonr/n<<" NSCS\n";
	exit(0);
      }
#endif
      return cpl*(qf*qf*(R-S));
    }

  };// end of class XS_qqG_NSCS

}// end of namespace EXTRAXS

DECLARE_TREEME2_GETTER(XS_qqG_NSCS,"XS_qqG_NSCS")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,XS_qqG_NSCS>::
operator()(const External_ME_Args& args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=5) return NULL;
  if (fl[0].IsLepton() && fl[1]==fl[0].Bar() && fl[2].IsGluon() &&
      fl[3].IsQuark() && fl[4]==fl[3].Bar()) {
    if (args.m_orders[0]==1 && args.m_orders[1]==2) {
      return new XS_qqG_NSCS(args);
    }
  }
  return NULL;
}
