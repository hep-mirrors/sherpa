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

    Complex A53(const int j1,const int j2,const int j3,
		const int j4, const int j5,
		const double &musq,const double &mt)
    {
      Complex Fax=-zb(j5,j3)*zb(j3,j1)*za(j2,j4)
	*(L1(-s(j1,j2),-s(j4,j5))/sqr(s(j4,j5))
	  -1.0/(12.0*s(j4,j5)*sqr(mt))/*
	  (1.+(2.*s(j1,j2)+s(j4,j5))/15./sqr(mt)
	   +(2.*s(j1,j2)*s(j4,j5)+3.*sqr(s(j1,j2))
	   +sqr(s(j4,j5)))/140./pow(mt,4))*/);
      Fax*=Complex(0.,-1.);
      return Fax;
    }

    void A5NLO(const int j1,const int j2,const int j3,
	       const int j4, const int j5,
	       const double &musq,const double &mt,
	       Complex &A5LOm,METOOLS::DivArrC &A5NLOm,Complex &A5ax)
    {
      A5LOm=-sqr(za(j1,j4))/(za(j2,j5)*za(j5,j1)*za(j4,j3));
      A5LOm*=Complex(0.,-1.);
      A5NLOm=A51(j2,j5,j1,j4,j3,musq)+A52(j2,j1,j5,j4,j3,musq)/9.0;
      if (mt) A5ax=A53(j2,j1,j5,j4,j3,musq,mt)/3.0;
      else A5ax=0.0;
    }

    DivArrD Virt5(const int j1,const int j2,const int j3,
		  const int j4, const int j5,
		  const double &musq,const double &mt,
		  double &born,Complex &Vax)
    {
      Complex A5LOm, A5LOp, A5axm, A5axp;
      METOOLS::DivArrC A5NLOm, A5NLOp;
      m_swap=1;
      A5NLO(j1,j2,j3,j4,j5,musq,mt,A5LOm,A5NLOm,A5axm);
      m_swap=0;
      A5NLO(j2,j1,j4,j3,j5,musq,mt,A5LOp,A5NLOp,A5axp);
      METOOLS::DivArrD res;
      res.IR2()=(std::conj(A5LOp)*A5NLOp.IR2()+std::conj(A5LOm)*A5NLOm.IR2()).real();
      res.IR()=(std::conj(A5LOp)*A5NLOp.IR()+std::conj(A5LOm)*A5NLOm.IR()).real();
      res.Finite()=(std::conj(A5LOp)*A5NLOp.Finite()+std::conj(A5LOm)*A5NLOm.Finite()).real();
      Vax=3.0*(std::conj(A5LOp)*A5axp+std::conj(A5LOm)*A5axm);
      born=(std::conj(A5LOp)*A5LOp+std::conj(A5LOm)*A5LOm).real();
      return 3.0*res;
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
  
  class QQGW_QCD_Virtual: public Virtual_ME2_Base,
			  public VJ_Amplitude {
  private:
    double m_nf, m_mw, m_ww;
    Complex m_gw, m_vax;
  public:
    QQGW_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi, flavs), m_nf(Flavour(kf_jet).Size()/2.0)
    {
      m_drmode=1;
      msg_Info()<<"QQWG"<<flavs<<"\n";
      m_mw=Flavour(kf_Wplus).Mass();
      m_ww=Flavour(kf_Wplus).Width();
      double g1(sqrt(4.*M_PI*s_model->ScalarConstant("alpha_QED")));
      m_gw=g1/sqrt(s_model->ComplexConstant("csin2_thetaW"));
      m_nf=(Flavour(kf_jet).Size()-1)/2;
    }
    double Eps_Scheme_Factor(const Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Compute(const Vec4D_Vector &p,const double &norm)
    {
      PreCompute(p);
      m_res=Virt5(1,2,3,4,5,m_mur2,0.0,m_born,m_vax);
      double prop=sqr(s(3,4))/(sqr(s(3,4)-sqr(m_mw))+sqr(m_mw*m_ww));
      double fac=2.0*4.0/3.0*3.0*4.0*M_PI*AlphaQCD()*prop;
      fac*=std::abs(sqr(m_gw*m_gw));
      m_res*=fac/norm;
      m_born*=fac/norm;
      DivArrD subuv;
      subuv.IR()=3.0*(11.0-2.0/3.0*m_nf)/6.0*m_born;
      subuv.Finite()=-3.0/6.0*m_born;
      m_res=(m_res-subuv)/m_born;
      msg_Debugging()<<"B     = "<<m_born<<"\n";
      msg_Debugging()<<"V_fin = "<<m_res.Finite()<<"\n";
      msg_Debugging()<<"V_e1  = "<<m_res.IR()<<"\n";
      msg_Debugging()<<"V_e2  = "<<m_res.IR2()<<"\n";
    }
  };// end of class QQGW_QCD_Virtual

  class GWQQ_QCD_Virtual: public QQGW_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    GWQQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGW_QCD_Virtual(pi,flavs)
    {
      m_flipq=!flavs[3].IsAnti();
      m_flipl=flavs[0].IsAnti();
      msg_Info()<<"GWQQ"<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=momenta[3+m_flipq];
      p[1]=momenta[4-m_flipq];
      p[2]=m_flipl?-momenta[0]:momenta[2];
      p[3]=m_flipl?momenta[2]:-momenta[0];
      p[4]=-momenta[1];
      Compute(p,4.0*8.0);
    }
  };// end of class GWQQ_QCD_Virtual

  class QWGQ_QCD_Virtual: public QQGW_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    QWGQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGW_QCD_Virtual(pi,flavs)
    {
      m_flipq=flavs[1].IsAnti();
      m_flipl=!flavs[0].IsAnti();
      msg_Info()<<"QWGQ"<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=m_flipq?momenta[4]:-momenta[1];
      p[1]=m_flipq?-momenta[1]:momenta[4];
      p[2]=m_flipl?momenta[2]:-momenta[0];
      p[3]=m_flipl?-momenta[0]:momenta[2];
      p[4]=momenta[3];
      Compute(p,4.0*3.0);
    }
  };// end of class QWGQ_QCD_Virtual

  class QQWG_QCD_Virtual: public QQGW_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    QQWG_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGW_QCD_Virtual(pi,flavs)
    {
      msg_Info()<<"QQWG"<<flavs<<"\n";
      m_flipq=flavs[0].IsAnti();
      m_flipl=flavs[2].IsAnti();
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=-momenta[0+m_flipq];
      p[1]=-momenta[1-m_flipq];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=momenta[4];
      Compute(p,4.0*9.0);
    }
  };// end of class QQWG_QCD_Virtual

  class GQWQ_QCD_Virtual: public QQGW_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    GQWQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGW_QCD_Virtual(pi,flavs)
    {
      msg_Info()<<"GQWQ"<<flavs<<"\n";
      m_flipq=flavs[1].IsAnti();
      m_flipl=flavs[2].IsAnti();
    }
    void Calc(const Vec4D_Vector& momenta)
    {
      Vec4D_Vector p(5);
      p[0]=m_flipq?momenta[4]:-momenta[1];
      p[1]=m_flipq?-momenta[1]:momenta[4];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=-momenta[0];
      Compute(p,4.0*24.0);
    }
  };// end of class GQWQ_QCD_Virtual

  class QQGZ_QCD_Virtual: public Virtual_ME2_Base,
			  public VJ_Amplitude {
  protected:
    int m_useax;
    double m_nf, m_g1, m_mt, m_mz, m_wz, m_ql, m_qq;
    Complex m_ll, m_rl, m_lq, m_rq;
  public:
    QQGZ_QCD_Virtual(const Process_Info &pi,const Flavour_Vector &flavs):
      Virtual_ME2_Base(pi, flavs), m_nf(Flavour(kf_jet).Size()/2.0) {}
    void Init(const Flavour &qfl,const Flavour &lfl)
    {
      m_drmode=1;
      Settings &s(Settings::GetMainSettings());
      m_useax=s["EXTRAXS"]["ZJ_AXIAL_PARTS"].
	GetScalarWithOtherDefault<int>(0);
      m_mt=Flavour(kf_t).Mass();
      m_mz=Flavour(kf_Z).Mass();
      m_wz=Flavour(kf_Z).Width();
      m_g1=sqrt(4.0*M_PI*s_model->ScalarConstant("alpha_QED"));
      Complex xw=s_model->ComplexConstant("csin2_thetaW");
      Complex sin2w=2.0*sqrt(xw*(1.0-xw));
      m_ql=lfl.Charge();
      m_qq=qfl.Charge();
      m_rq=-2.0*qfl.Charge()*xw/sin2w;
      m_lq=m_rq+2.0*qfl.IsoWeak()/sin2w;
      m_rl=-2.0*lfl.Charge()*xw/sin2w;
      m_ll=m_rl+2.0*lfl.IsoWeak()/sin2w;
      m_nf=(Flavour(kf_jet).Size()-1)/2;
    }
    double Eps_Scheme_Factor(const Vec4D_Vector& mom)
    {
      return 4.*M_PI;
    }
    void Compute(const Vec4D_Vector &p,const double &norm)
    {
      PreCompute(p);
      Complex prop=s(3,4)/Complex(s(3,4)-sqr(m_mz),m_mz*m_wz);
      double fac=8.0*4.0/3.0*3.0*sqr(m_g1*m_g1)*4.0*M_PI*AlphaQCD();
      fac/=norm;
      double BLL, BLR, BRL, BRR;
      Complex VaxLL, VaxRL, VaxLR, VaxRR;
      DivArrD LL=Virt5(1,2,3,4,5,m_mur2,m_mt,BLL,VaxLL);
      DivArrD LR=Virt5(1,2,4,3,5,m_mur2,m_mt,BLR,VaxLR);
      DivArrD RL=Virt5(2,1,3,4,5,m_mur2,m_mt,BRL,VaxRL);
      DivArrD RR=Virt5(2,1,4,3,5,m_mur2,m_mt,BRR,VaxRR);
      Complex cLL(m_qq*m_ql+m_lq*m_ll*prop);
      Complex cLR(m_qq*m_ql+m_lq*m_rl*prop);
      Complex cRL(m_qq*m_ql+m_rq*m_ll*prop);
      Complex cRR(m_qq*m_ql+m_rq*m_rl*prop);
      m_res=fac*(sqr(std::abs(cLL))*LL
		 +sqr(std::abs(cLR))*LR
		 +sqr(std::abs(cRL))*RL
		 +sqr(std::abs(cRR))*RR);
      m_born=fac*(sqr(std::abs(cLL))*BLL
		  +sqr(std::abs(cLR))*BLR
		  +sqr(std::abs(cRL))*BRL
		  +sqr(std::abs(cRR))*BRR);
      Complex caxL(std::conj((m_rq-m_lq)*m_ll*prop));
      Complex caxR(std::conj((m_rq-m_lq)*m_rl*prop));
      if (m_useax) {
	m_res.Finite()+=fac*
	  (caxL*(cLL*VaxLL+cRL*VaxRL)+
	   caxR*(cLR*VaxLR+cRR*VaxRR)).real();
      }
      DivArrD subuv;
      subuv.IR()=3.0*(11.0-2.0/3.0*m_nf)/6.0*m_born;
      subuv.Finite()=-3.0/6.0*m_born;
      m_res=(m_res-subuv)/m_born;
      msg_Debugging()<<"B     = "<<m_born<<"\n";
      msg_Debugging()<<"V_fin = "<<m_res.Finite()<<"\n";
      msg_Debugging()<<"V_e1  = "<<m_res.IR()<<"\n";
      msg_Debugging()<<"V_e2  = "<<m_res.IR2()<<"\n";
    }
  };// end of class QQGZ_QCD_Virtual

  class ZGQQ_QCD_Virtual: public QQGZ_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    ZGQQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGZ_QCD_Virtual(pi,flavs)
    {
      m_flipq=!flavs[3].IsAnti();
      m_flipl=!flavs[0].IsAnti();
      Init(flavs[m_flipq?3:4],flavs[m_flipl?0:1]);
      msg_Info()<<"ZGQQ"<<m_useax<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=momenta[3+m_flipq];
      p[1]=momenta[4-m_flipq];
      p[2]=-momenta[0+m_flipl];
      p[3]=-momenta[1-m_flipl];
      p[4]=momenta[2];
      Compute(p,4.0);
    }
  };// end of class ZGQQ_QCD_Virtual

  class GZQQ_QCD_Virtual: public QQGZ_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    GZQQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGZ_QCD_Virtual(pi,flavs)
    {
      m_flipq=!flavs[3].IsAnti();
      m_flipl=flavs[0].IsAnti();
      Init(m_flipq?flavs[3]:flavs[3].Bar(),
	   m_flipl?flavs[0].Bar():flavs[0]);
      msg_Info()<<"GZQQ"<<m_useax<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=momenta[3+m_flipq];
      p[1]=momenta[4-m_flipq];
      p[2]=m_flipl?-momenta[0]:momenta[2];
      p[3]=m_flipl?momenta[2]:-momenta[0];
      p[4]=-momenta[1];
      Compute(p,4.0*8.0);
    }
  };// end of class GZQQ_QCD_Virtual

  class QZGQ_QCD_Virtual: public QQGZ_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    QZGQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGZ_QCD_Virtual(pi,flavs)
    {
      m_flipq=flavs[1].IsAnti();
      m_flipl=!flavs[0].IsAnti();
      Init(m_flipq?flavs[1].Bar():flavs[1],
	   m_flipl?flavs[0]:flavs[0].Bar());
      msg_Info()<<"QZGQ"<<m_useax<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=m_flipq?momenta[4]:-momenta[1];
      p[1]=m_flipq?-momenta[1]:momenta[4];
      p[2]=m_flipl?momenta[2]:-momenta[0];
      p[3]=m_flipl?-momenta[0]:momenta[2];
      p[4]=momenta[3];
      Compute(p,4.0*3.0);
    }
  };// end of class QZGQ_QCD_Virtual

  class QQZG_QCD_Virtual: public QQGZ_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
  public:
    QQZG_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGZ_QCD_Virtual(pi,flavs)
    {
      m_flipq=flavs[0].IsAnti();
      m_flipl=flavs[2].IsAnti();
      Init(flavs[m_flipq?1:0],flavs[m_flipl?3:2]);
      msg_Info()<<"QQZG"<<m_useax<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector &momenta)
    {
      Vec4D_Vector p(5);
      p[0]=-momenta[0+m_flipq];
      p[1]=-momenta[1-m_flipq];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=momenta[4];
      Compute(p,4.0*9.0);
    }
  };// end of class QQZG_QCD_Virtual

  class GQZQ_QCD_Virtual: public QQGZ_QCD_Virtual {
  private:
    bool m_flipq, m_flipl;
    double m_nf, m_g1, m_mz, m_wz, m_ql, m_qq;
    Complex m_ll, m_rl, m_lq, m_rq;
  public:
    GQZQ_QCD_Virtual(const Process_Info &pi,
		     const Flavour_Vector &flavs):
      QQGZ_QCD_Virtual(pi,flavs)
    {
      m_flipq=flavs[1].IsAnti();
      m_flipl=flavs[2].IsAnti();
      Init(m_flipq?flavs[1].Bar():m_flavs[1],flavs[m_flipl?3:2]);
      msg_Info()<<"GQZQ"<<m_useax<<flavs<<"\n";
    }
    void Calc(const Vec4D_Vector& momenta)
    {
      Vec4D_Vector p(5);
      p[0]=m_flipq?momenta[4]:-momenta[1];
      p[1]=m_flipq?-momenta[1]:momenta[4];
      p[2]=momenta[2+m_flipl];
      p[3]=momenta[3-m_flipl];
      p[4]=-momenta[0];
      Compute(p,4.0*24.0);
    }
  };// end of class GQZQ_QCD_Virtual

}

using namespace EXTRAXS;

DECLARE_VIRTUALME2_GETTER(QQWG_QCD_Virtual,"W_Jet_QCD_Virtual")
Virtual_ME2_Base *ATOOLS::Getter
<Virtual_ME2_Base,Process_Info,QQWG_QCD_Virtual>::
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
    if (((fl[0].IsLepton() && fl[2].IsNeutrino()) ||
	 (fl[2].IsLepton() && fl[0].IsNeutrino())) &&
	fl[2].LeptonFamily()==fl[0].LeptonFamily()) {
      if (fl[1].IsGluon() && fl[3].IsQuark() &&
	  fl[3].QuarkFamily()==fl[4].QuarkFamily())
	return new GWQQ_QCD_Virtual(pi,fl);
      if (fl[3].IsGluon() && fl[1].IsQuark() &&
	  fl[4].QuarkFamily()==fl[1].QuarkFamily())
	return new QWGQ_QCD_Virtual(pi,fl);
    }
    if ((fl[2].IsLepton() && fl[3].IsNeutrino() &&
	 fl[3].LeptonFamily()==fl[2].LeptonFamily())) {
      if (fl[4].IsGluon() && fl[0].IsQuark() &&
	  fl[1].QuarkFamily()==fl[0].QuarkFamily())
	return new QQWG_QCD_Virtual(pi,fl);
      if (fl[0].IsGluon() && fl[1].IsQuark() &&
	  fl[4].QuarkFamily()==fl[1].QuarkFamily())
	return new GQWQ_QCD_Virtual(pi,fl);
    }
  }
  return NULL;
}

DECLARE_VIRTUALME2_GETTER(QQZG_QCD_Virtual,"Z_Jet_QCD_Virtual")
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
    if ((fl[0].IsLepton() && fl[1]==fl[0].Bar()) ||
	(fl[0].IsNeutrino() && fl[1]==fl[0].Bar())) {
      if (fl[2].IsGluon() && fl[3].IsQuark() && fl[4]==fl[3].Bar())
	return new ZGQQ_QCD_Virtual(pi,fl);
    }
    if ((fl[0].IsLepton() && fl[2]==fl[0]) ||
	(fl[0].IsNeutrino() && fl[2]==fl[0])) {
      if (fl[1].IsGluon() && fl[3].IsQuark() && fl[4]==fl[3].Bar())
	return new GZQQ_QCD_Virtual(pi,fl);
      if (fl[3].IsGluon() && fl[1].IsQuark() && fl[4]==fl[1])
	return new QZGQ_QCD_Virtual(pi,fl);
    }
    if ((fl[2].IsLepton() && fl[3]==fl[2].Bar()) ||
	(fl[2].IsNeutrino() && fl[3]==fl[2].Bar())) {
      if (fl[4].IsGluon() && fl[0].IsQuark() && fl[1]==fl[0].Bar())
	return new QQZG_QCD_Virtual(pi,fl);
      if (fl[0].IsGluon() && fl[1].IsQuark() && fl[4]==fl[1])
	return new GQZQ_QCD_Virtual(pi,fl);
    }
  }
  return NULL;
}
