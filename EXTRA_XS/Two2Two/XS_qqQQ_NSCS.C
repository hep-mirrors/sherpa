#include "PHASIC++/Process/External_ME_Args.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "MODEL/UFO/UFO_Model.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/Random.H"

#include "EXTRA_XS/Main/ME2_Base.H"
#include "EXTRA_XS/Two2Two/NSCS_Tools.H"

using namespace EXTRAXS;
using namespace MODEL;
using namespace ATOOLS;
using namespace PHASIC;
using namespace NSCS;

namespace EXTRAXS {

  class XS_qqQQ_NSCS : public ME2_Base {
  private:

    double qe, qf1, qf2;
    double alpha, alpha_S, colfac;

  public:

    XS_qqQQ_NSCS(const External_ME_Args& args): ME2_Base(args)
    {
      DEBUG_INFO("now entered EXTRAXS::XS_qqQQ_NSCS ...");
      m_sintt=1;
      m_oew=2;
      m_oqcd=2;
      alpha = MODEL::s_model->ScalarConstant("alpha_QED");
      alpha_S = MODEL::s_model->ScalarConstant("alpha_S");
      qe = m_flavs[0].Charge();
      qf1 = m_flavs[2].Charge();
      qf2 = m_flavs[3].Charge();
      colfac=3.*4./3.*.5;
      m_cfls[3].push_back(kf_photon);
      m_cfls[3].push_back(kf_Z);
      m_cfls[60].push_back(kf_photon);
      m_cfls[60].push_back(kf_Z);
    }

    double CosPhi12ij(const double &z1,const double &z2,const double &z3,
		      const double &s12,const double &s13,const double &s23)
    {
      return sqr(z1+z2)/(4.*z1*z2)*sqr(z1*s23-z2*s13)/(s12*(s13+s23)*(z1+z2)*z3);
    }

    double DCqgDCgq(const double &z1,const double &z2,const double &z3,const double &cp2)
    {
      return 2.*z3/(1.-z3)*(1.-(4.*z1*z2)/sqr(z1+z2)*cp2)
	+(1.-z3)*(1.-2.*z1*z2/pow(z1+z2,2));
    }

    double DSDCqgDCgq(const double &z1,const double &z2,const double &z3,const double &cp2)
    {
      return 2.*z3/(1.-z3)*(1.-(4.*z1*z2)/sqr(z1+z2)*cp2);
    }

    double TCqQQb(const double &z1,const double &z2,const double &z3,
		 const double &s12,const double &s13,const double &s23) const
    {
      double t123(2.*(z1*s23-z2*s13)/(z1+z2)+s12*(z1-z2)/(z1+z2));
      return 0.5*((4.*z3+sqr(z1-z2))/(z1+z2)-sqr(t123)/(s12*(s12+s13+s23))
		  +(z1+z2-s12/(s12+s13+s23)));
    }

    double DSTCqQQb(const double &z1,const double &z2,const double &z3,
		   const double &s12,const double &s13,const double &s23) const
    {
      double t123(2.*(z1*s23-z2*s13)/(z1+z2));
      return 0.5*(s12+s13+s23)/(s13+s23)*(4.*z3/(z1+z2)-sqr(t123)/(s12*(s13+s23)));
    }

    double DS(const Vec4D &p1,const Vec4D &p2,
	      const Vec4D &ip3,const Vec4D &ip4)
    {
      Vec4D p3(ip3), p4(ip4);
      if (p3[0]<p4[0]) std::swap<Vec4D>(p3,p4);
      double s34(2.*p3*p4), s14(2.*p1*p4), s13(2.*p1*p3);
      double s12(2.*p1*p2), s23(2.*p2*p3), s24(2.*p2*p4);
      double Q2(s12+s13+s14+s23+s24+s34);
      NSCS_Args DS(p1,p2,p3,p4);
      ClusterNNLO(DS);
      double PS(PSWeight(DS,2));
      DS.m_x1=0.;
      double DSPS(PSWeight(DS,2)/PS);
      double DSW1314(W1314(DS));
      double DSW1324(W1324(DS));
      double SB(s12/(s34*(s13+s14)*(s23+s24)));
      double SA(SB*sqr(s13*s24-s23*s14)/(s12*s34*(s13+s14)*(s23+s24)));
      double DSME(2.*(SB-SA));
      DEBUG_VAR(s12*DSME<<" "<<DSW1314+DSW1324<<" "<<DSPS<<" "<<DS.m_sec);
      return s12*DSME*DSPS*(DSW1314+DSW1324);
    }

    double TC(const Vec4D &p1,const Vec4D &p2,
	      const Vec4D &ip3,const Vec4D &ip4)
    {
      Vec4D p3(ip3), p4(ip4);
      if (p3[0]<p4[0]) std::swap<Vec4D>(p3,p4);
      NSCS_Args CC1(p1,p2,p3,p4);
      ClusterNNLO(CC1);
      double PS(PSWeight(CC1,2));
      CC1.m_x3=0.;
      ConstructNNLO(CC1);
      double CC1PS(PSWeight(CC1,2)/PS);
      CC1.m_x1=0.;
      double DSCC1PS(PSWeight(CC1,2)/PS);
      double Q((p1+p2+p3+p4).Mass()), x1(CC1.m_p[0][0]/(Q/2.));
      double x3(CC1.m_p[2][0]/(Q/2.)), x4(CC1.m_p[3][0]/(Q/2.));
      if (x1<0. || x3<0. || x4<0. || CC1.m_p[1][0]<0.) return 0.;
      double s12(2.*p1*p2), s34(2.*p3*p4), s134((p1+p3+p4).Abs2());
      double s13(s34/CC1.m_s34s13), s14(s13*CC1.m_s14s13);
      double TC1(TCqQQb(x3,x4,x1,s34,s13,s14)/(s34*s134));
      double DSTC1(DSTCqQQb(x3,x4,x1,s34,s13,s14)/(s34*s134));
      DEBUG_VAR(Q*Q*TC1<<" "<<s12*DSTC1<<" "<<CC1PS<<" "<<DSCC1PS<<" "<<CC1.m_sec);
      return Q*Q*TC1*CC1PS-s12*DSTC1*DSCC1PS;
    }

    double TCDC(const Vec4D &p1,const Vec4D &p2,
		const Vec4D &ip3,const Vec4D &ip4)
    {
      Vec4D p3(ip3), p4(ip4);
      if (p3[0]<p4[0]) std::swap<Vec4D>(p3,p4);
      NSCS_Args CC1C34(p1,p2,p3,p4);
      ClusterNNLO(CC1C34);
      double PS(PSWeight(CC1C34,2));
      if (!(CC1C34.m_sec==2 || CC1C34.m_sec==4)) return 0.;
      NSCS_Args C34(CC1C34);
      C34.m_x4=0.;
      ConstructNNLO(C34);
      CC1C34.m_x3=0.;
      ConstructNNLO(CC1C34);
      double CC1C34PS(PSWeight(CC1C34,2)/PS);
      CC1C34.m_x1=0.;
      double DSCC1C34PS(PSWeight(CC1C34,2)/PS);
      double Q((p1+p2+p3+p4).Mass()), x1(CC1C34.m_p[0][0]/(Q/2.));
      double x3(CC1C34.m_p[2][0]/(Q/2.)), x4(CC1C34.m_p[3][0]/(Q/2.));
      if (x1<0. || x3<0. || x4<0. || CC1C34.m_p[1][0]<0.) return 0.;
      Vec4D p1t(C34.m_p[0]), p34t(C34.m_p[2]+C34.m_p[3]);
      double s12(2.*p1*p2), s34(2.*p3*p4), s134(2.*p1t*p34t);
      double s13(s34/CC1C34.m_s34s13), s14(s13*CC1C34.m_s14s13);
      double cp2(CosPhi12ij(x3,x4,x1,s34,s13,s14));
      double TC1DC34(DCqgDCgq(x3,x4,x1,cp2)/(s34*s134));
      double DSTC1DC34(DSDCqgDCgq(x3,x4,x1,cp2)/(s34*s134));
      DEBUG_VAR(Q*Q*TC1DC34<<" "<<s12*DSTC1DC34<<" "<<CC1C34PS<<" "<<DSCC1C34PS<<" "<<CC1C34.m_sec);
      return Q*Q*TC1DC34*CC1C34PS-s12*DSTC1DC34*DSCC1C34PS;
    }

    double DC(const Vec4D &p1,const Vec4D &p2,
	      const Vec4D &ip3,const Vec4D &ip4)
    {
      Vec4D p3(ip3), p4(ip4);
      if (p3[0]<p4[0]) std::swap<Vec4D>(p3,p4);
      NSCS_Args C34(p1,p2,p3,p4);
      ClusterNNLO(C34);
      if (!(C34.m_sec==2 || C34.m_sec==4)) return 0.;
      double s34(2.*p3*p4), Q2((p1+p2+p3+p4).Abs2());
      double PS(PSWeight(C34,2));
      C34.m_x4=0.;
      ConstructNNLO(C34);
      double C34PS(PSWeight(C34,2)/PS);
      C34.m_x1=0.;
      double DSC34PS(PSWeight(C34,2)/PS);
      Vec4D p1t(C34.m_p[0]), p2t(C34.m_p[1]);
      Vec4D p34t(C34.m_p[2]+C34.m_p[3]);
      if (p1t[0]<0. || p2t[0]<0. || p34t[0]<0.) return 0.;
      Vec4D kt(C34.m_b*sqrt(C34.m_lam)+C34.m_a*sqrt(1.-C34.m_lam));
      double sb13(2.*p1t*p34t), sb23(2.*p2t*p34t), sb12(2.*p1t*p2t);
      double z(C34.m_p[2][0]/(C34.m_p[2][0]+C34.m_p[3][0]));
      double DD(4*z*(1-z)/kt.Abs2());
      double C34ME=4*Q2*sb12/(sb13*sb23)
	+4*Q2*DD*sqr((kt*p1t)/sb13-(kt*p2t)/sb23)
        +2*(sb13/sb23+sb23/sb13)*(1.-DD/2.*kt.Abs2());
      double DSC34ME=4*sb12*sb12/(sb13*sb23)
	+4*sb12*DD*sqr((kt*p1t)/sb13-(kt*p2t)/sb23);
      double C34W1314(W1314(C34,1));
      C34ME/=2.*s34;
      DSC34ME/=2.*s34;
      DEBUG_VAR(C34ME<<" "<<DSC34ME<<" "<<C34W1314<<" "<<C34PS<<" "<<DSC34PS<<" "<<C34.m_sec);
      return (C34ME*C34PS-DSC34ME*DSC34PS)*C34W1314;
    }

    double ERTD(const Vec4D &p1,const Vec4D &p2,
		const Vec4D &p3,const Vec4D &p4)
    {
      double s34(2.*p3*p4), s14(2.*p1*p4), s13(2.*p1*p3);
      double s12(2.*p1*p2), s23(2.*p2*p3), s24(2.*p2*p4);
      double s134(s13+s14+s34), s234(s23+s24+s34);
#ifdef SOFT_UNSAFE_ME
      double ERTD = (2*s12*s13*s14 + s13*s14*s23 - pow(s14,2)*s23 - pow(s13,2)*s24
		 + s13*s14*s24 + s12*s13*s34 + s12*s14*s34 + s13*s23*s34 + s14*s24*s34)/
      (pow(s34,2)*pow(s13 + s14 + s34,2)) + (-(s14*pow(s23,2)) + 2*s12*s23*s24
       + s13*s23*s24 + s14*s23*s24 - s13*pow(s24,2) + s12*s23*s34 + s13*s23*s34 + s12*s24*s34 +
        s14*s24*s34)/(pow(s34,2)*pow(s23 + s24 + s34,2)) - (s12*s14*s23 - s13*s14*s23
        + s14*pow(s23,2) + s12*s13*s24 + pow(s13,2)*s24 - s13*s23*s24 - pow(s12,2)*s34 -
        s12*s13*s34 - s12*s23*s34 - s12*pow(s34,2))/(pow(s34,2)*(s13 + s14 + s34)*(s23 + s24 + s34)) -
     (s12*s14*s23 + pow(s14,2)*s23 + s12*s13*s24 - s13*s14*s24 - s14*s23*s24
      + s13*pow(s24,2) - pow(s12,2)*s34 - s12*s14*s34 - s12*s24*s34 - s12*pow(s34,2))/
      (pow(s34,2)*(s13 + s14 + s34)*(s23 + s24 + s34));
#else
      double ERTD=2.*s12*s12/(s134*s234*s34)
	*(1.-sqr(s13*s24-s14*s23)/(s134*s234*s34*s12)
	  -((s14*s23-s13*s24)*(s14-s24-s13+s23)
	    -s34*(s13-s23)*(s14-s24))/(s134*s234*s12));
      ERTD+=(s13/s134+s13/s234)*(s12*s34+s14*s23-s13*s24)/(s134*s34*s34)
	+(s14/s134+s23/s234)*(s12*s34+s13*s24-s14*s23)/(s134*s34*s34)
	+s12/(s134*s234)+(s13*s23+s14*s24)/(s134*s134*s34);
      ERTD+=(s24/s234+s24/s134)*(s12*s34+s23*s14-s24*s13)/(s234*s34*s34)
	+(s23/s234+s14/s134)*(s12*s34+s24*s13-s23*s14)/(s234*s34*s34)
	+s12/(s234*s134)+(s24*s14+s23*s13)/(s234*s234*s34);
      DEBUG_VAR(ERTD);
#endif
      return ERTD;
    }

#ifndef CHECK_LIMIT
    double operator()(const Vec4D_Vector &p)
    {
#else
    double operator()(const Vec4D_Vector &ip)
    {
      static double x1(-1.), x2, x3, x4;
      static int sec(2), mode(0), pts(1000);
      if (x1<0.) {
	Settings& s = Settings::GetMainSettings();
	Scoped_Settings nscs{ s["NSCS"] };
	nscs["x1"].SetDefault(1.);
	nscs["x2"].SetDefault(1.);
	nscs["x3"].SetDefault(1.);
	nscs["x4"].SetDefault(1.);
	nscs["pts"].SetDefault(10000000);
	nscs["mode"].SetDefault(0);
	nscs["sec"].SetDefault(2);
	x1=nscs["x1"].Get<double>();
	x2=nscs["x2"].Get<double>();
	x3=nscs["x3"].Get<double>();
	x4=nscs["x4"].Get<double>();
	pts=nscs["pts"].Get<int>();
	mode=nscs["mode"].Get<int>();
	sec=nscs["sec"].Get<int>();
	DEBUG_VAR(x1<<" "<<x2<<" "<<x3<<" "<<x4<<" "<<mode<<" "<<sec);
      }
      Vec4D_Vector p={ip[0],ip[1],ip[3],ip[5],ip[2],ip[4]};
      if ((p[4][0]>p[2][0] && p[4][0]>p[2][0]) ||
	  (p[5][0]>p[3][0] && p[5][0]>p[3][0])) {
	std::swap<Vec4D>(p[4],p[2]);
	std::swap<Vec4D>(p[5],p[3]);
      }
      NSCS_Args C(p[2],p[3],p[4],p[5]);
      C.m_sec=sec;
      ClusterNNLO(C);
      if (x1<1.) C.m_x1=x1;
      if (x2<1.) C.m_x2=x2;
      if (x3<1.) C.m_x3=x3;
      if (x4<1.) C.m_x4=x4;
      ConstructNNLO(C);
      p[3]=C.m_p[0];
      p[5]=C.m_p[1];
      p[2]=C.m_p[2];
      p[4]=C.m_p[3];
      if (p[2][0]<0. || p[3][0]<0. || p[4][0]<0. || p[5][0]<0.) return 0.;
#endif
      double M24=ERTD(p[3],p[5],p[2],p[4]);
      DEBUG_VAR(M24);
      double DS3=DS(p[3],p[5],p[2],p[4]);
      double DS5=DS(p[5],p[3],p[2],p[4]);
      double TC3=TC(p[3],p[5],p[2],p[4]);
      double TC5=TC(p[5],p[3],p[2],p[4]);
      double DC324=DC(p[3],p[5],p[2],p[4]);
      double DC524=DC(p[5],p[3],p[2],p[4]);
      double TCDC3=TCDC(p[3],p[5],p[2],p[4]);
      double TCDC5=TCDC(p[5],p[3],p[2],p[4]);
      double RS24=M24-DS3-DS5-TC3-TC5+TCDC3+TCDC5-DC324-DC524;
      DEBUG_VAR(RS24<<" "<<RS24/M24);
      double M35=ERTD(p[2],p[4],p[3],p[5]);
      DEBUG_VAR(M35);
      double DS2=DS(p[2],p[4],p[3],p[5]);
      double DS4=DS(p[4],p[2],p[3],p[5]);
      double TC2=TC(p[2],p[4],p[3],p[5]);
      double TC4=TC(p[4],p[2],p[3],p[5]);
      double DC235=DC(p[2],p[4],p[3],p[5]);
      double DC435=DC(p[4],p[2],p[3],p[5]);
      double TCDC2=TCDC(p[2],p[4],p[3],p[5]);
      double TCDC4=TCDC(p[4],p[2],p[3],p[5]);
      double RS35=M35-DS2-DS4-TC2-TC4+TCDC2+TCDC4-DC235-DC435;
      DEBUG_VAR(RS35<<" "<<RS35/M35);
      double cpl=sqr(4.*M_PI*alpha)*sqr(8.*M_PI*alpha_S)
	*CouplingFactor(2,2)*colfac*sqr(qe);
      cpl*=1./(p[0]+p[1]).Abs2();
      if (IsBad(RS24+RS35)) return 0.;
#ifdef CHECK_LIMIT
      static double sumr(0.), sums(0.), sumrs(0.), sumrsonr(0.), n(0.);
      double psf(1.0);
      if (mode) {
	psf=x1*x1*x1*x2*x3*x4;
	if (mode&2) psf*=x1*x2*x3*x4;
      }
      sumr+=psf*cpl*(qf2*qf2*M24+qf1*qf1*M35);
      sums+=psf*cpl*(qf2*qf2*(M24-RS24)+qf1*qf1*(M35-RS35));
      sumrs+=psf*cpl*dabs(qf2*qf2*RS24+qf1*qf1*RS35);
      DEBUG_VAR(psf*cpl*dabs(qf2*qf2*RS24+qf1*qf1*RS35));
      sumrsonr+=dabs((qf2*qf2*RS24+qf1*qf1*RS35)/(qf2*qf2*M24+qf1*qf1*M35));
      n+=1.;
      if ((int(n)%pts==0)) {
	msg_Info()<<sumr/n<<" "<<sums/n<<" "<<sumrs/n<<" "<<sumrsonr/n<<" NSCS\n";
	exit(0);
      }
#endif
      // return cpl*(qf2*qf2*M24+qf1*qf1*M35);
      return cpl*(qf2*qf2*RS24+qf1*qf1*RS35);
    }

  };// end of class XS_qqQQ_NSCS

}// end of namespace EXTRAXS

DECLARE_TREEME2_GETTER(XS_qqQQ_NSCS,"XS_qqQQ_NSCS")
Tree_ME2_Base *ATOOLS::Getter<Tree_ME2_Base,External_ME_Args,XS_qqQQ_NSCS>::
operator()(const External_ME_Args& args) const
{
  if (dynamic_cast<UFO::UFO_Model*>(MODEL::s_model)) return NULL;

  static int mode(-1);
  if (mode==-1) {
    Settings& s = Settings::GetMainSettings();
    Scoped_Settings nscs{ s["NSCS"] };
    nscs["ME"].SetDefault(0);
    mode=nscs["ME"].Get<int>();
  }
  if (mode==1) return NULL;

  const Flavour_Vector fl=args.Flavours();
  if (fl.size()!=6) return NULL;
  if ((fl[0].IsLepton() && fl[1]==fl[0].Bar() && fl[2].IsQuark() &&
       fl[3].IsQuark() && fl[2]!=fl[3] &&
       fl[4]==fl[2].Bar() && fl[5]==fl[3].Bar())) {
    if (args.m_orders[0]==2 && args.m_orders[1]==2) {
      return new XS_qqQQ_NSCS(args);
    }
  }
  return NULL;
}
