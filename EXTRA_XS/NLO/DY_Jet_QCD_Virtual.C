#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "METOOLS/Loops/Divergence_Array.H"
#include "MODEL/Main/Model_Base.H"
#include "EXTRA_XS/NLO/Logarithms.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Scoped_Settings.H"

using namespace PHASIC;
using namespace METOOLS;
using namespace MODEL;
using namespace ATOOLS;

namespace EXTRAXS {

  class VJ_Amplitude {
  private:
    int m_swap;
    double m_s[5][5];
    Complex m_za[5][5], m_zb[5][5];
  public:

    inline const double &s(const int i,const int j) const
    { return m_s[i-1][j-1]; }
    inline const Complex &za(const int i,const int j) const
    { return m_swap?m_zb[i-1][j-1]:m_za[i-1][j-1]; }
    inline const Complex &zb(const int i,const int j) const
    { return m_swap?m_za[i-1][j-1]:m_zb[i-1][j-1]; }

    DivArrC A51(const int j1,const int j2,const int j3,
		const int j4, const int j5,const double &musq)
    {
      Complex A5lom=-sqr(za(j3,j4))/(za(j1,j2)*za(j2,j3)*za(j4,j5));
      Complex l12=LnRat(musq,-s(j1,j2)), l23=LnRat(musq,-s(j2,j3));
      DivArrC Vcc, Vsc;
      Complex Fcc, Fsc;
      Vcc.IR2()=-2.0;
      Vcc.IR()=-2.0-(l12+l23);
      Vcc.Finite()=-0.5*(sqr(l12)+sqr(l23))-2.0*l23-4.0;
      Fcc=sqr(za(j3,j4))/(za(j1,j2)*za(j2,j3)*za(j4,j5))
	*(Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
	  -2.0*za(j3,j1)*zb(j1,j5)*za(j5,j4)/za(j3,j4)
	  *L0(-s(j2,j3),-s(j4,j5))/s(j4,j5));
      Vsc.IR()=0.5;
      Vsc.Finite()=0.5*l23+1.0;
      Fsc=za(j3,j4)*za(j3,j1)*zb(j1,j5)*za(j5,j4)
	/(za(j1,j2)*za(j2,j3)*za(j4,j5))*L0(-s(j2,j3),-s(j4,j5))/s(j4,j5)
	+0.5*sqr(za(j3,j1)*zb(j1,j5))*za(j4,j5)
	/(za(j1,j2)*za(j2,j3))*L1(-s(j2,j3),-s(j4,j5))/sqr(s(j4,j5));
      A5lom*=Complex(0.,-1.);
      Fcc*=Complex(0.,-1.);
      Fsc*=Complex(0.,-1.);
      return (Vcc+Vsc)*A5lom+Fcc+Fsc;
    }

    DivArrC A52(const int j1,const int j2,const int j3,
		const int j4, const int j5,const double &musq)
    {
      Complex A5lom=sqr(za(j2,j4))/(za(j2,j3)*za(j3,j1)*za(j4,j5));
      Complex l12=LnRat(musq,-s(j1,j2)), l45=LnRat(musq,-s(j4,j5));
      DivArrC Vcc, Vsc;
      Complex Fcc, Fsc;
      Vcc.IR2()=-1.0;
      Vcc.IR()=-2.0-l12;
      Vcc.Finite()=-0.5*sqr(l12)-2.0*l45-4.0;
      Fcc=-sqr(za(j2,j4))/(za(j2,j3)*za(j3,j1)*za(j4,j5))
	*Lsm1(-s(j1,j2),-s(j4,j5),-s(j1,j3),-s(j4,j5))
	+za(j2,j4)*(za(j1,j2)*za(j3,j4)-za(j1,j4)*za(j2,j3))
	/(za(j2,j3)*sqr(za(j1,j3))*za(j4,j5))
	*Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
	+2.0*zb(j1,j3)*za(j1,j4)*za(j2,j4)/(za(j1,j3)*za(j4,j5))
	*L0(-s(j2,j3),-s(j4,j5))/s(j4,j5);
      Vsc.IR()=0.5;
      Vsc.Finite()=0.5*l45+0.5;
      Fsc=sqr(za(j1,j4))*za(j2,j3)/(pow(za(j1,j3),3)*za(j4,j5))
	*Lsm1(-s(j1,j2),-s(j4,j5),-s(j2,j3),-s(j4,j5))
	-0.5*sqr(za(j4,j1)*zb(j1,j3))*za(j2,j3)/(za(j1,j3)*za(j4,j5))
	*L1(-s(j4,j5),-s(j2,j3))/sqr(s(j2,j3))
	+sqr(za(j1,j4))*za(j2,j3)*zb(j3,j1)/(sqr(za(j1,j3))*za(j4,j5))
	*L0(-s(j4,j5),-s(j2,j3))/s(j2,j3)
	-za(j2,j1)*zb(j1,j3)*za(j4,j3)*zb(j3,j5)/za(j1,j3)
	*L1(-s(j4,j5),-s(j1,j2))/sqr(s(j1,j2))
	-za(j2,j1)*zb(j1,j3)*za(j3,j4)*za(j1,j4)/(sqr(za(j1,j3))*za(j4,j5))
	*L0(-s(j4,j5),-s(j1,j2))/s(j1,j2)
	-0.5*zb(j3,j5)*(zb(j1,j3)*zb(j2,j5)+zb(j2,j3)*zb(j1,j5))
	/(zb(j1,j2)*zb(j2,j3)*za(j1,j3)*zb(j4,j5));
      A5lom*=Complex(0.,-1.);
      Fcc*=Complex(0.,-1.);
      Fsc*=Complex(0.,-1.);
      return (Vcc+Vsc)*A5lom+Fcc+Fsc;
    }

    DivArrC A53(const int j1,const int j2,const int j3,
		const int j4, const int j5,
		const double &musq,const double &mt)
    {
      DivArrC Fax;
      Fax.Finite()=-zb(j5,j3)*zb(j3,j1)*za(j2,j4)
	*(L1(-s(j1,j2),-s(j4,j5))/sqr(s(j4,j5))
	  -1.0/(12.0*s(j4,j5)*sqr(mt)));
      return Fax;
    }

    void A5NLO(const int j1,const int j2,const int j3,
	       const int j4, const int j5,
	       const double &musq,const double &mt,
	       Complex &A5LOm,METOOLS::DivArrC &A5NLOm)
    {
      A5LOm=-sqr(za(j1,j4))/(za(j2,j5)*za(j5,j1)*za(j4,j3));
      A5LOm*=Complex(0.,-1.);
      A5NLOm=A51(j2,j5,j1,j4,j3,musq)+A52(j2,j1,j5,j4,j3,musq)/9.0;
    }

    DivArrD Virt5(const double &musq,const double &as,const double &nf,
		  const double &gw,const double &mw,const double &ww,
		  const double &mt)
    {
      Complex A5LOm, A5LOp;
      METOOLS::DivArrC A5NLOm, A5NLOp;
      m_swap=1;
      A5NLO(1,2,3,4,5,musq,mt,A5LOm,A5NLOm);
      m_swap=0;
      A5NLO(2,1,4,3,5,musq,mt,A5LOp,A5NLOp);
      METOOLS::DivArrD res, subuv;
      res.IR2()=(std::conj(A5LOp)*A5NLOp.IR2()+std::conj(A5LOm)*A5NLOm.IR2()).real();
      res.IR()=(std::conj(A5LOp)*A5NLOp.IR()+std::conj(A5LOm)*A5NLOm.IR()).real();
      res.Finite()=(std::conj(A5LOp)*A5NLOp.Finite()+std::conj(A5LOm)*A5NLOm.Finite()).real();
      double born((std::conj(A5LOp)*A5LOp+std::conj(A5LOm)*A5LOm).real());
      double prop=sqr(s(3,4))/(sqr(s(3,4)-sqr(mw))+sqr(mw*ww));
      double fac=2.0*4.0/3.0*3.0*sqr(gw*gw)*4.0*M_PI*as*prop;
      res*=3.0*fac;
      born*=fac;
      subuv.IR()=3.0*(11.0-2.0/3.0*nf)/6.0*born;
      subuv.Finite()=-3.0/6.0*born;
      res=res-subuv;
      msg_Debugging()<<"B     = "<<born<<"\n";
      msg_Debugging()<<"V_fin = "<<res.Finite()/born<<"\n";
      msg_Debugging()<<"V_e1  = "<<res.IR()/born<<"\n";
      msg_Debugging()<<"V_e2  = "<<res.IR2()/born<<"\n";
      return res/born;
    }

    void PreCompute(const Vec4D_Vector &p)
    {
      DEBUG_FUNC("");
      for (size_t i(0);i<5;++i)
	msg_Debugging()<<"p["<<i<<"]=Vec4("<<p[i]<<");\n";
      for (size_t i(0);i<5;++i) {
	Spinor<double> spi(1,p[i]), smi(-1,p[i]); 
	for (size_t j(i+1);j<5;++j) {
	  Spinor<double> spj(1,p[j]), smj(-1,p[j]);
	  m_za[j][i]=-(m_za[i][j]=spj*spi);
	  m_zb[j][i]=-(m_zb[i][j]=smi*smj);
	  m_s[j][i]=m_s[i][j]=2.0*p[i]*p[j];
	  msg_Debugging()<<"<"<<i<<","<<j<<"> = "<<m_za[i][j]
			 <<", ["<<i<<","<<j<<"] = "<<m_zb[i][j]
			 <<", s_{"<<i<<j<<"} = "<<m_s[i][j]<<"\n";
	}
      }
    }    
    
  };// end of class VJ_Amplitude
  
  const double NC(3.0), CA(3.0), CF(0.5*(NC*NC-1.0)/NC), TR(0.5);

  class QQZG_QCD_Virtual: public Virtual_ME2_Base,
			  public VJ_Amplitude {
  private:
    bool m_flipq, m_flipl;
    double m_nf, m_gw, m_mw, m_ww;
  public:
    QQZG_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi, flavs), m_nf(Flavour(kf_jet).Size()/2.0)
    {
      m_drmode=1;
      m_flipq=flavs[0].IsAnti();
      m_flipl=flavs[2].IsAnti();
      msg_Info()<<"QQZG"<<flavs<<"\n";
      double gf(1.0/sqrt(2.0)/std::abs(sqr(s_model->ComplexConstant("cvev"))));
      m_mw=Flavour(kf_Wplus).Mass();
      m_ww=Flavour(kf_Wplus).Width();
      m_gw=sqrt(8.0*sqr(m_mw)*gf/sqrt(2.0));
      m_nf=(Flavour(kf_jet).Size()-1)/2;
    }
    double Eps_Scheme_Factor(const Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=-momenta[0+m_flipq];
      p[1]=-momenta[1-m_flipq];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=momenta[4];
      PreCompute(p);
      m_res=Virt5(m_mur2,AlphaQCD(),m_nf,m_gw,m_mw,m_ww,0.0);
    }
  };// end of class QQZG_QCD_Virtual

  class GQZQ_QCD_Virtual: public Virtual_ME2_Base,
			  public VJ_Amplitude {
  private:
    bool m_flipq, m_flipl;
    double m_nf, m_gw, m_mw, m_ww;
  public:
    GQZQ_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi, flavs)
    {
      m_drmode=1;
      m_flipq=flavs[1].IsAnti();
      m_flipl=flavs[2].IsAnti();
      msg_Info()<<"GQZQ"<<flavs<<"\n";
      double gf(1.0/sqrt(2.0)/std::abs(sqr(s_model->ComplexConstant("cvev"))));
      m_mw=Flavour(kf_Wplus).Mass();
      m_ww=Flavour(kf_Wplus).Width();
      m_gw=sqrt(8.0*sqr(m_mw)*gf/sqrt(2.0));
      m_nf=(Flavour(kf_jet).Size()-1)/2;
    }
    double Eps_Scheme_Factor(const Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Calc(const Vec4D_Vector& momenta)
    {
      Vec4D_Vector p(5);
      p[0]=m_flipq?momenta[4]:-momenta[1];
      p[1]=m_flipq?-momenta[1]:momenta[4];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=-momenta[0];
      PreCompute(p);
      m_res=Virt5(m_mur2,AlphaQCD(),m_nf,m_gw,m_mw,m_ww,0.0);
    }
  };// end of class GQZQ_QCD_Virtual

}

using namespace EXTRAXS;

DECLARE_VIRTUALME2_GETTER(QQZG_QCD_Virtual,"DY_Jet_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,QQZG_QCD_Virtual>::
operator()(const Process_Info &pi) const
{
  DEBUG_FUNC(pi);
  if (pi.m_loopgenerator!="Internal") return NULL;
  if (pi.m_fi.m_nlotype&nlo_type::loop) {
    if (pi.m_mincpl[0]!=2. || pi.m_maxcpl[0]!=2.) return NULL;
    if (pi.m_mincpl[1]!=2. || pi.m_maxcpl[1]!=2.) return NULL;
    if (pi.m_fi.m_nlocpl[0]!=1. || pi.m_fi.m_nlocpl[1]!=0.) return NULL;
    if (pi.m_fi.m_ps.size()!=3) return NULL;
    Flavour_Vector fl=pi.ExtractFlavours();
    for (size_t i(0);i<fl.size();++i)
      if (fl[i].Mass()) return NULL;
    if (!(fl[2].IsLepton() && fl[3].IsNeutrino() &&
	  fl[3].LeptonFamily()==fl[2].LeptonFamily())) return NULL;
    if (fl[4].IsGluon() && fl[0].IsQuark() &&
	fl[1].QuarkFamily()==fl[0].QuarkFamily())
      return new QQZG_QCD_Virtual(pi,fl);
    if (fl[0].IsGluon() && fl[1].IsQuark() &&
	fl[4].QuarkFamily()==fl[1].QuarkFamily())
      return new GQZQ_QCD_Virtual(pi,fl);
  }
  return NULL;
}
