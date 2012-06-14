#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "AddOns/MCFM/MCFM_Wrapper.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Phys/Color.H"
#include "MODEL/Main/Standard_Model.H"

namespace MCFM {

  class MCFM_qqb_zgam: public PHASIC::Virtual_ME2_Base {
  private:
    int                       m_pID;
    double                  * p_p, *p_msqv;
    //MODEL::Running_AlphaS   * p_as;
    double                    m_aem;
    //double                    m_MZ;
    //double                    m_GZ;
    //double                    m_GF;
    //double                    m_Qu;
    //double                    m_Qd;
    //double                    m_CF;
    //double                    m_NC; 
    //double                    m_normcorr;
    //double                    m_sin2_thetaW;
    //double                    Q[5];
    bool                      m_anom;
    double                    m_h_G_1;
    double                    m_h_G_2;
    double                    m_h_G_3;
    double                    m_h_G_4;
    double                    m_h_Z_1;
    double                    m_h_Z_2;
    double                    m_h_Z_3;
    double                    m_h_Z_4;
    double                    m_unitarization_scale;
   
  public:
    MCFM_qqb_zgam(const int & pID,const PHASIC::Process_Info& pi,
		  const Flavour_Vector& flavs, bool anom);
    ~MCFM_qqb_zgam();
    void Calc(const ATOOLS::Vec4D_Vector& momenta);
    double Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom);
    //Complex * zvirtamp(int p1,int p2,int p3,int p4,int p5,
    //			Complex * za,Complex * zb);
    //Complex * zamp(int p1,int p2,int p3,int p4,int p5,
    //		    Complex * za,Complex * zb);
    //void DoMCFM();

  };

}
// declare functions taken from MCFM
extern "C" { 
  //  Complex vpole_(double * s);
  void spinoru_(const int & N, double *p,Complex * za, Complex * zb);
  //Complex zvirtamps_(const int & i1,const int & i2,const int & i3, 
  //		     const int & i4,const int & i5, 
  //		     Complex * za,Complex * zb, Complex * qqnlo,
  //		     Complex * qqlo);
  //Complex fagamma_(const int & i1, const int & i2,const int & i3, const int & i4,
  //	   const int & i5, Complex * za,Complex * zb);
  void qqb_zgam_v_(double * p, double *a);
  //void qqb_zgam_(double * p, double *a);
  //Complex fbgamma_(const int & i1, const int & i2,const int & i3, const int & i4,
  //		   const int & i5, Complex * za,Complex * zb);
  //Complex zamps_(const int * p1, const int * p2, const int * p3, const int * p4,
  //		 const int * p5, Complex * za, Complex * zb, Complex *a);
  //extern struct{
  //  double L[MCFM_NF],R[MCFM_NF],le,ln,re,rn;
  //  double sin2w,q1,l1,r1,q2,l2,r2;
  //}zcouple_;

}

#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace MCFM;
using namespace PHASIC;
using namespace ATOOLS;
using namespace MODEL;
using namespace std;


MCFM_qqb_zgam::MCFM_qqb_zgam(const int & pID,const PHASIC::Process_Info& pi,
			     const Flavour_Vector& flavs, bool anom):
  Virtual_ME2_Base(pi,flavs), m_pID(pID),
  //p_as((Running_AlphaS *)s_model->GetScalarFunction(string("alpha_S"))),
  m_aem(s_model->ScalarFunction(string("alpha_QED"))),
  // m_MZ(Flavour(kf_Z).Mass()), m_GZ(Flavour(kf_Z).Width()),
  //m_GF(1.16639e-5), m_Qu(Flavour(kf_u).Charge()), m_Qd(Flavour(kf_d).Charge()),
  //m_CF(CF), m_NC(NC), m_normcorr(1.),
  // m_sin2_thetaW(s_model->ScalarConstant(string("sin2_thetaW"))),
  //Q({-1./3.,2./3.,-1./3.,2./3.,-1./3.}),
  m_anom(anom),
  m_h_G_1(0), m_h_G_2(0), m_h_G_3(0), m_h_G_4(0), m_h_Z_1(0), m_h_Z_2(0),
  m_h_Z_3(0), m_h_Z_4(0),
  m_unitarization_scale(0)
     
{
  rpa->gen.AddCitation
    (1,"The NLO matrix elements have been taken from MCFM \\cite{}.");
  
  p_p      = new double[4*MCFM_NMX];
  p_msqv   = new double[sqr(2*MCFM_NF+1)];
  m_drmode = m_mode=1;
}

MCFM_qqb_zgam::~MCFM_qqb_zgam()
{
  delete [] p_p;
  delete [] p_msqv;
}
/*
Complex * MCFM_qqb_zgam::zvirtamp(int p1,int p2,int p3,int p4,int p5,
				 Complex * za,Complex * zb)
{
  // make indices that convert from 2d fortran matrix to 1d c++ array
  int p1mx = (p1-1)*MCFM_NMX-1;
  int p2mx = (p2-1)*MCFM_NMX-1;
  int p3mx = (p3-1)*MCFM_NMX-1;
  int p4mx = (p4-1)*MCFM_NMX-1;
  int p5mx = (p5-1)*MCFM_NMX-1;
  Complex vpl;
  Complex prp34;
 
  Complex * qqbnlo = new Complex[16];
  Complex * qqblo  = 0;
  //Complex qqblo[16];
  Complex a[16];
  double s34;
  qqblo = zamp(p1,p2,p3,p4,p5,za,zb);
  //zamps_(&p1,&p2,&p3,&p4,&p5,za,zb,a);
  //std::cout << "qqblo Sherpa:  " << qqblo[0] << std::endl; 
  //std::cout << "qqblo MCFM:    " << a[0]     << std::endl;

  s34=sprods_.s[p3+p4mx];
  prp34=s34/Complex(s34-pow(m_MZ,2),m_MZ*m_GZ);
  vpl=vpole_(&sprods_.s[p1+p2mx]);
  //std::cout << "vpole:  " << vpl << std::endl;
  
  for (int j=0; j<=1; j++){
    qqbnlo[8*j+4*0+2*0+1]=
      (Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp34)
       *(fagamma_(p1,p2,p3,p4,p5,za,zb)+fbgamma_(p1,p2,p3,p4,p5,za,zb)))
      +vpl*qqblo[8*j+4*0+2*0+1];
    
    qqbnlo[8*j+4*0+2*1+1]=
      (Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp34)
       *(fagamma_(p1,p2,p4,p3,p5,za,zb)+fbgamma_(p1,p2,p4,p3,p5,za,zb)))
      +vpl*qqblo[8*j+4*0+2*1+1];
    
    qqbnlo[8*j+4*0+2*0+0]=
      (Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp34)
	*(fagamma_(p2,p1,p4,p3,p5,zb,za)+fbgamma_(p2,p1,p4,p3,p5,zb,za)))
      +vpl*qqblo[8*j+4*0+2*0+0];
    
    qqbnlo[8*j+4*0+2*1+0]=
      (Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp34)
	*(fagamma_(p2,p1,p3,p4,p5,zb,za)+fbgamma_(p2,p1,p3,p4,p5,zb,za)))
      +vpl*qqblo[8*j+4*0+2*1+0];
    
    qqbnlo[8*j+4*1+2*1+0]=      
      -(Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp34)
	*(fagamma_(p1,p2,p3,p4,p5,zb,za)+fbgamma_(p1,p2,p3,p4,p5,zb,za)))
      +vpl*qqblo[8*j+4*1+2*1+0];
    
    qqbnlo[8*j+4*1+2*0+0]=      
      -(Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp34)
	*(fagamma_(p1,p2,p4,p3,p5,zb,za)+fbgamma_(p1,p2,p4,p3,p5,zb,za)))
      +vpl*qqblo[8*j+4*1+2*0+0];
    
    qqbnlo[8*j+4*1+2*1+1]=      
      -(Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp34)
	*(fagamma_(p2,p1,p4,p3,p5,za,zb)+fbgamma_(p2,p1,p4,p3,p5,za,zb)))
      +vpl*qqblo[8*j+4*1+2*1+1]; 
    
    qqbnlo[8*j+4*1+2*0+1]=      
      -(Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp34)
	*(fagamma_(p2,p1,p3,p4,p5,za,zb)+fbgamma_(p2,p1,p3,p4,p5,za,zb)))
      +vpl*qqblo[8*j+4*1+2*0+1];
  }
  delete [] qqblo;
  return qqbnlo;
}


Complex * MCFM_qqb_zgam::zamp(int p1,int p2,int p3,int p4,int p5, 
			       Complex * za,Complex * zb)
{
  int p1mx = (p1-1)*MCFM_NMX-1;
  int p2mx = (p2-1)*MCFM_NMX-1;
  int p3mx = (p3-1)*MCFM_NMX-1;
  int p4mx = (p4-1)*MCFM_NMX-1;
  int p5mx = (p5-1)*MCFM_NMX-1;

  int j;
  double s12,s34,xfac;
  Complex ai[2][2][2],af[2][2][2],prp34,prp12;
  Complex * qqbamp = new Complex[16];
  Complex anomZZ[2][2][2],anomgamZ[2][2][2],anomZgam[2][2][2],cprop;
  Complex im1=Complex(0.,1.);
  s12=sprods_.s[p1+p2mx];
  //std::cout << "s12 S:  " << s12 << std::endl;
  s34=sprods_.s[p3+p4mx];
  //std::cout << "s23 S: " << s34 << std::endl;
  if (m_anom){
    anomcoup_.tevscale = s_model->ScalarConstant(string("UNITARIZATION_SCALE"));
    anomcoup_.h1gam    = s_model->ScalarConstant(string("h1_gamma"));
    anomcoup_.h2gam    = s_model->ScalarConstant(string("h2_gamma"));
    anomcoup_.h3gam    = s_model->ScalarConstant(string("h3_gamma"));
    anomcoup_.h4gam    = s_model->ScalarConstant(string("h4_gamma"));
    anomcoup_.h1Z      = s_model->ScalarConstant(string("h1_Z"));
    anomcoup_.h2Z      = s_model->ScalarConstant(string("h2_Z"));
    anomcoup_.h3Z      = s_model->ScalarConstant(string("h3_Z"));
    anomcoup_.h4Z      = s_model->ScalarConstant(string("h4_Z"));
  }
  xfac=1./(1.+s12/pow((anomcoup_.tevscale),2));
 
  anomcoup_.h1tZ=pow(xfac,3)*anomcoup_.h1Z/pow(m_MZ,2);
  anomcoup_.h2tZ=pow(xfac,4)*anomcoup_.h2Z/pow(m_MZ,4);
  anomcoup_.h3tZ=pow(xfac,3)*anomcoup_.h3Z/pow(m_MZ,2);
  anomcoup_.h4tZ=pow(xfac,4)*anomcoup_.h4Z/pow(m_MZ,4);
  anomcoup_.h1tgam=pow(xfac,3)*anomcoup_.h1gam/pow(m_MZ,2);
  anomcoup_.h2tgam=pow(xfac,4)*anomcoup_.h2gam/pow(m_MZ,4);
  anomcoup_.h3tgam=pow(xfac,3)*anomcoup_.h3gam/pow(m_MZ,2);
  anomcoup_.h4tgam=pow(xfac,4)*anomcoup_.h4gam/pow(m_MZ,4);

  ai[0][0][1]=pow(za[p1+p3mx],2)*zb[p3+p4mx]/
    (za[p1+p5mx]*za[p2+p5mx]*s34);
  ai[0][1][1]=pow(za[p1+p4mx],2)*zb[p4+p3mx]/
    (za[p1+p5mx]*za[p2+p5mx]*s34);
  ai[0][0][0]=pow(zb[p2+p4mx],2)*za[p4+p3mx]/
    (zb[p1+p5mx]*zb[p2+p5mx]*s34); 
  ai[0][1][0]=pow(zb[p2+p3mx],2)*za[p3+p4mx]/
    (zb[p1+p5mx]*zb[p2+p5mx]*s34);
  // std::cout << "za[p1+3mx]  S   " << za[p1+p3mx]  << std::endl;
  // std::cout << "ai[0][0][1] S:  " << ai[0][0][1] << std::endl;
  anomZZ[0][0][1]=0.25/s34*
    ((im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*za[p1+p3mx]*
     zb[p2+p5mx]*zb[p4+p5mx]
     +(im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*za[p1+p2mx]*
     za[p3+p5mx]*zb[p5+p4mx]*pow(zb[p2+p5mx],2));
  
  anomZZ[0][1][1]=0.25/s34*
    ((im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*za[p1+p4mx]*
     zb[p2+p5mx]*zb[p3+p5mx]
     +(im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*za[p1+p2mx]*
     za[p4+p5mx]*zb[p5+p3mx]*pow(zb[p2+p5mx],2));
  
  anomZZ[0][0][0]=0.25/s34*
    ((-im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*zb[p2+p4mx]*
     za[p1+p5mx]*za[p3+p5mx]
     +(-im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*zb[p2+p1mx]*
     zb[p4+p5mx]*za[p5+p3mx]*pow(za[p1+p5mx],2));
  
  anomZZ[0][1][0]=0.25/s34*
    ((-im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*zb[p2+p3mx]*
     za[p1+p5mx]*za[p4+p5mx]
     +(-im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*zb[p2+p1mx]*
     zb[p3+p5mx]*za[p5+p4mx]*pow(za[p1+p5mx],2));

  anomgamZ[0][0][1]=0.25/s34*
    ((im1*anomcoup_.h1tgam+anomcoup_.h3tgam)*2.*za[p1+p3mx]*
     zb[p2+p5mx]*zb[p4+p5mx]+
     (im1*anomcoup_.h2tgam+anomcoup_.h4tgam)*za[p1+p2mx]*
     za[p3+p5mx]*zb[p5+p4mx]*pow(zb[p2+p5mx],2));

  anomgamZ[0][1][1]=0.25/s34*
    ((im1*anomcoup_.h1tgam+anomcoup_.h3tgam)*2.*za[p1+p4mx]*
     zb[p2+p5mx]*zb[p3+p5mx]+
     (im1*anomcoup_.h2tgam+anomcoup_.h4tgam)*za[p1+p2mx]*
     za[p4+p5mx]*zb[p5+p3mx]*pow(zb[p2+p5mx],2));

  anomgamZ[0][0][0]=0.25/s34*
    ((-im1*anomcoup_.h1tgam+anomcoup_.h3tgam)*2.*zb[p2+p4mx]*
     za[p1+p5mx]*za[p3+p5mx]+
     (-im1*anomcoup_.h2tgam+anomcoup_.h4tgam)*zb[p2+p1mx]*
     zb[p4+p5mx]*za[p5+p3mx]*pow(za[p1+p5mx],2));

  anomgamZ[0][1][0]=0.25/s34*
    ((-im1*anomcoup_.h1tgam+anomcoup_.h3tgam)*2.*zb[p2+p3mx]*
     za[p1+p5mx]*za[p4+p5mx]+
     (-im1*anomcoup_.h2tgam+anomcoup_.h4tgam)*zb[p2+p1mx]*
     zb[p3+p5mx]*za[p5+p4mx]*pow(za[p1+p5mx],2));

  anomZgam[0][0][1]=0.25/s12*
    ((im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*za[p3+p1mx]*
     zb[p4+p5mx]*zb[p2+p5mx]+
     (im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*za[p3+p4mx]*
     za[p1+p5mx]*zb[p5+p2mx]*pow(zb[p4+p5mx],2));

  anomZgam[0][1][1]=0.25/s12*
    ((im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*za[p4+p1mx]*
     zb[p3+p5mx]*zb[p2+p5mx]+
     (im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*za[p4+p3mx]*
     za[p1+p5mx]*zb[p5+p2mx]*pow(zb[p3+p5mx],2));

  anomZgam[0][0][0]=0.25/s12*
    ((-im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*zb[p4+p2mx]*
     za[p3+p5mx]*za[p1+p5mx]+
     (-im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*zb[p4+p3mx]*
     zb[p2+p5mx]*za[p5+p1mx]*pow(za[p3+p5mx],2));

  anomZgam[0][1][0]=0.25/s12*
    ((-im1*anomcoup_.h1tZ+anomcoup_.h3tZ)*2.*zb[p3+p2mx]*
     za[p4+p5mx]*za[p1+p5mx]+
     (-im1*anomcoup_.h2tZ+anomcoup_.h4tZ)*zb[p3+p4mx]*
     zb[p2+p5mx]*za[p5+p1mx]*pow(za[p4+p5mx],2));
 
  if (zerowidth_.zerowidth){
    // std::cout << "Sherpa zerowidth" << std::endl;
    af[0][0][0]=Complex(0.,0.);
    af[1][0][0]=Complex(0.,0.);
    af[0][0][1]=Complex(0.,0.);
    af[0][1][0]=Complex(0.,0.);
  }
  else{
    af[0][0][1]=pow(za[p1+p3mx],2)*zb[p1+p2mx]/
      (za[p3+p5mx]*za[p4+p5mx]*s12);
    af[1][0][1]=pow(za[p2+p3mx],2)*zb[p2+p1mx]/
      (za[p3+p5mx]*za[p4+p5mx]*s12);
    af[0][0][0]=pow(zb[p2+p4mx],2)*za[p2+p1mx]/
      (zb[p3+p5mx]*zb[p4+p5mx]*s12);
    af[1][0][0]=pow(zb[p1+p4mx],2)*za[p1+p2mx]/
      (zb[p3+p5mx]*zb[p4+p5mx]*s12);
    //std::cout << "af[0][0][1] S:  " << af[0][0][1] <<std::endl;
  }

  //Others determined by complex conjugation
  ai[1][1][0]=conj(ai[0][0][1]);
  ai[1][0][0]=conj(ai[0][1][1]);
  ai[1][0][1]=conj(ai[0][1][0]);
  ai[1][1][1]=conj(ai[0][0][0]);

  af[1][1][0]=conj(af[0][0][1]);
  af[0][1][0]=conj(af[1][0][1]);
  af[1][1][1]=conj(af[0][0][0]);
  af[0][1][1]=conj(af[1][0][0]);

  anomZZ[1][1][0]=conj(anomZZ[0][0][1]);
  anomZZ[1][0][0]=conj(anomZZ[0][1][1]);
  anomZZ[1][0][1]=conj(anomZZ[0][1][0]);
  anomZZ[1][1][1]=conj(anomZZ[0][0][0]);

  anomgamZ[1][1][0]=conj(anomgamZ[0][0][1]);
  anomgamZ[1][0][0]=conj(anomgamZ[0][1][1]);
  anomgamZ[1][0][1]=conj(anomgamZ[0][1][0]);
  anomgamZ[1][1][1]=conj(anomgamZ[0][0][0]);

  anomZgam[1][1][0]=conj(anomZgam[0][0][1]);
  anomZgam[1][0][0]=conj(anomZgam[0][1][1]);
  anomZgam[1][0][1]=conj(anomZgam[0][1][0]);
  anomZgam[1][1][1]=conj(anomZgam[0][0][0]);

  prp34=s34/Complex((s34-pow(m_MZ,2)),m_MZ*m_GZ);
  //  std::cout << "prp34 S:  " << prp34 << std::endl;
  prp12=s12/Complex((s12-pow(m_MZ,2)),m_MZ*m_GZ);
  //std::cout << "prp12 S:  " << prp12 << std::endl;
  cprop=(s12-s34)/Complex((s12-pow(m_MZ,2)),m_MZ*m_GZ);
  //  std::cout << "cprop S:  " << cprop << std::endl; 

  for(j=0; j<=1; j++){
    qqbamp[8*j+4*0+2*0+0]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp34)*ai[0][0][0]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp12)*af[0][0][0]
      +zcouple_.L[j]*zcouple_.l1*prp34*anomZZ[0][0][0]*cprop
      +Q[j]*zcouple_.l1*prp34*anomgamZ[0][0][0] 
      +zcouple_.L[j]*zcouple_.q1*prp12*anomZgam[0][0][0];
    std::cout << "l1 S " << zcouple_.l1 << std::endl;
    std::cout << "q1 S " << zcouple_.q1 << std::endl;
    std::cout << "000 for " << j << " S:  " << 
     Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1) << std::endl;
    qqbamp[8*j+4*0+2*0+1]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp34)*ai[0][0][1]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.l1*prp12)*af[0][0][1]
      +zcouple_.L[j]*zcouple_.l1*prp34*anomZZ[0][0][1]*cprop
      +Q[j]*zcouple_.l1*prp34*anomgamZ[0][0][1]
      +zcouple_.L[j]*zcouple_.q1*prp12*anomZgam[0][0][1];

    qqbamp[8*j+4*1+2*1+0]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp34)*ai[1][1][0]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp12)*af[1][1][0]
      +zcouple_.R[j]*zcouple_.r1*prp34*anomZZ[1][1][0]*cprop
      +Q[j]*zcouple_.r1*prp34*anomgamZ[1][1][0]
      +zcouple_.R[j]*zcouple_.q1*prp12*anomZgam[1][1][0];
    qqbamp[8*j+4*1+2*1+1]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp34)*ai[1][1][1]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.r1*prp12)*af[1][1][1]
      +zcouple_.R[j]*zcouple_.r1*prp34*anomZZ[1][1][1]*cprop
      +Q[j]*zcouple_.r1*prp34*anomgamZ[1][1][1]
      +zcouple_.R[j]*zcouple_.q1*prp12*anomZgam[1][1][1];

    qqbamp[8*j+4*0+2*1+0]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp34)*ai[0][1][0]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp12)*af[0][1][0]
      +zcouple_.L[j]*zcouple_.r1*prp34*anomZZ[0][1][0]*cprop
      +Q[j]*zcouple_.r1*prp34*anomgamZ[0][1][0]
      +zcouple_.L[j]*zcouple_.q1*prp12*anomZgam[0][1][0];
    qqbamp[8*j+4*0+2*1+1]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp34)*ai[0][1][1]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.L[j]*zcouple_.r1*prp12)*af[0][1][1]
      +zcouple_.L[j]*zcouple_.r1*prp34*anomZZ[0][1][1]*cprop
      +Q[j]*zcouple_.r1*prp34*anomgamZ[0][1][1]
      +zcouple_.L[j]*zcouple_.q1*prp12*anomZgam[0][1][1];

    qqbamp[8*j+4*1+2*0+0]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp34)*ai[1][0][0]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp12)*af[1][0][0]
      +zcouple_.R[j]*zcouple_.l1*prp34*anomZZ[1][0][0]*cprop
      +Q[j]*zcouple_.l1*prp34*anomgamZ[1][0][0]
      +zcouple_.R[j]*zcouple_.q1*prp12*anomZgam[1][0][0];
    qqbamp[8*j+4*1+2*0+1]=
      Q[j]*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp34)*ai[1][0][1]+
      zcouple_.q1*(Q[j]*zcouple_.q1+zcouple_.R[j]*zcouple_.l1*prp12)*af[1][0][1]
      +zcouple_.R[j]*zcouple_.l1*prp34*anomZZ[1][0][1]*cprop
      +Q[j]*zcouple_.l1*prp34*anomgamZ[1][0][1]
      +zcouple_.R[j]*zcouple_.q1*prp12*anomZgam[1][0][1];
			     
  }
  return qqbamp;
}

void MCFM_qqb_zgam::DoMCFM()
{
  double fac;
  double qbq[2];
  double qqb[2];
  Complex * qbqnlo = 0;
  Complex * qqbnlo = 0;
  Complex * qqblo  = 0;
  Complex * qbqlo  = 0;
 
  qqblo  = zamp(1,2,3,4,5,zprods_.za,zprods_.zb);
  qqbnlo = zvirtamp(1,2,3,4,5,zprods_.za,zprods_.zb);
  Complex b[16];
  Complex c[16];
  //zvirtamps_(1,2,3,4,5,zprods_.za,zprods_.zb, b, c);
  //std::cout << "MCFM:  " << b[0] << std::endl;

  fac = m_CF*8.*pow(4*M_PI*m_aem,3)*m_NC;
  if (nproc_.nproc==305) fac=fac/3.;
  int j,h12,h34,h5;
 
  int jj[2*MCFM_NF + 1]={0,1,0,1,0,0,0,1,0,1,0};

  for (j=0; j<2; j++){
    qqb[j]=0.0;
    for (h12=0; h12<2; h12++){
      for (h34=0; h34<2; h34++){
	for (h5=0; h5<2; h5++){
	  qqb[j]=qqb[j] +
	    2.0*real(conj(qqblo[8*j+4*h12+2*h34+h5])*qqbnlo[8*j+4*h12+2*h34+h5]);
	}
      }
    }
    qqb[j]=fac*qqb[j];
  }
  for (int k=MCFM_NF; k<2*MCFM_NF+1; k++){
    if (k==MCFM_NF) p_msqv[k*(2*MCFM_NF+1)+k] = 0.0;
    else if (k>MCFM_NF){
      p_msqv[k*(2*MCFM_NF+1)-k+2*MCFM_NF]=qqb[jj[k]];
     }
  }
  //std::cout << "SHERPA:  " << qqbnlo[0] << std::endl;
  delete [] qqblo;
  delete [] qqbnlo;
}
*/
void MCFM_qqb_zgam::Calc(const Vec4D_Vector &p)
{
  // order partons from q, qb, gamma, e^-, e^+   to
  //                    q, qb, e^-  , e^+, gamma
  // Also make all momenta outgoing
  for (int n(0);n<2;++n)           GetMom(p_p,n,-p[n]);
  
  GetMom(p_p,2,p[3]);
  GetMom(p_p,3,p[4]);  
  GetMom(p_p,4,p[2]);

  // set up helicity amplitudes
  spinoru_(5,p_p,zprods_.za,zprods_.zb);

  long int i(m_flavs[0]), j(m_flavs[1]);
  if (i==21) { i=0; }// correct for gluons
  if (j==21) { j=0; }
  scale_.musq=m_mur2; //set scale for process
  scale_.scale=sqrt(scale_.musq);
  if (m_anom) {
    anomcoup_.h1Z = 10;
    anomcoup_.tevscale = 1000;
  }
  //double a[11][11];
  //qqb_zgam_(p_p, *a);
  // add MCFM_NF to indices such that they run from 0 
  i+=MCFM_NF;
  j+=MCFM_NF;
  double fac = 36*pow(4*M_PI*m_aem/ewcouple_.esq,3)/qcdcouple_.ason2pi;
  if (m_pID == 305) fac=fac/3.;
  // std::cout << "MCFM born  " << 
  //  a[j][i]*36*pow(4*M_PI*m_aem/ewcouple_.esq,3) << std::endl;
  // set residues
  epinv_.epinv=epinv2_.epinv2=0.0;
  qqb_zgam_v_(p_p,p_msqv);
  //DoMCFM();
  double resMCFM(p_msqv[j*(2*MCFM_NF+1) + i]*fac);
  //double res(p_msqv[i*(2*MCFM_NF+1) + j]);
  epinv_.epinv=1.0;
  qqb_zgam_v_(p_p,p_msqv);
  //DoMCFM(); 
  //double res1(p_msqv[i*(2*MCFM_NF+1) + j]);
  double res1MCFM(p_msqv[j*(2*MCFM_NF+1)+i]*fac);
  epinv2_.epinv2=1.0;
  qqb_zgam_v_(p_p,p_msqv);
  //DoMCFM();
  //double res2(p_msqv[i*(2*MCFM_NF+1) + j]);
  double res2MCFM(p_msqv[j*(2*MCFM_NF+1)+i]*fac);
  
  //std::cout << "MCFM double " << (res2MCFM-res1MCFM) 
  //*36*pow(4*M_PI*m_aem/ewcouple_.esq,3) << std::endl;
  //std::cout << "MCFM single " << (res1MCFM-resMCFM) 
  //*36*pow(4*M_PI*m_aem/ewcouple_.esq,3)  << std::endl;
  //m_res.Finite() = res;
  //m_res.IR()     = (res1-res);
  //m_res.IR2()    = (res2-res1);

  m_res.Finite() = resMCFM;
  m_res.IR()     = (res1MCFM-resMCFM);
  m_res.IR2()    = (res2MCFM-res1MCFM);
}

double MCFM_qqb_zgam::Eps_Scheme_Factor(const ATOOLS::Vec4D_Vector& mom)
{
  return 4.*M_PI;
}

extern "C" { void chooser_(); }

DECLARE_VIRTUALME2_GETTER(MCFM_qqb_zgam_Getter,"MCFM_qqb_zgam")
Virtual_ME2_Base *MCFM_qqb_zgam_Getter::operator()(const Process_Info &pi) const
{
  //std::cout <<"---zgam getter---" << std::endl;
  DEBUG_FUNC("");
  if (pi.m_loopgenerator!="MCFM")                       return NULL;
  if (MODEL::s_model->Name()!=std::string("SM")
      && MODEL::s_model->Name()!=std::string("SM+AGC")) return NULL;
  if (pi.m_fi.m_nloewtype!=nlo_type::lo)                return NULL;
  if (pi.m_fi.m_nloqcdtype&nlo_type::loop) {
    Flavour_Vector fl(pi.ExtractFlavours());
    // std::cout << fl[0] << " " << fl[1] << " " << fl[2] << " "
    //	      << fl[3] << " " << fl[4] << std::endl;
    // two incoming strongly interacting particles.
    if (!fl[0].Strong() || !fl[1].Strong())             return NULL;
    if (fl.size()!=5)                                   return NULL;
    // check for fully leptonic FS
    if (!(fl[0].IsQuark() && fl[1].IsQuark()))          return NULL;
    if (!(fl[4].Kfcode()==22 || fl[2].Kfcode()==22
	  || fl[3].Kfcode()==22))                       return NULL;
    if (!((fl[3].IsLepton() && fl[4]==fl[3].Bar()) ||
	(fl[2].IsLepton() && fl[3]==fl[2].Bar())))      return NULL;
    
    int pID(0);
    bool swapped(false);
    bool anom(false);

    if (fl[2].IsUptype() || fl[3].IsUptype()){
      pID=305;}
    else {pID=300; zerowidth_.zerowidth=false; swapped=true;}
    
    if (pi.m_fi.m_ps.size()==2) {
      ATOOLS::Flavour fl1(pi.m_fi.m_ps[0].m_fl[0]);
      ATOOLS::Flavour fl2(pi.m_fi.m_ps[1].m_fl[0]);
      
      if (fl[3].IsLepton() && fl[4]==fl[3].Bar()) {
	if (MODEL::s_model->Name()!=std::string("SM")
	    && MODEL::s_model->Name()!=std::string("SM+AGC")) {
	  msg_Error()<<"Warning in "<<METHOD<<":"<<std::endl
		     <<"   Try to initialise process qqb->Vgamma in MCFM."
		     <<std::endl
		     <<"   Inconsistent setting with Sherpa: "<<std::endl
		     <<"model = "<<MODEL::s_model->Name()<<"(should be 'SM')."
		     <<std::endl<<"   Will exit the run."<<std::endl;
	  exit(1);
	  return NULL;
	}
	
      }
    }
   
    if (MODEL::s_model->Name()==std::string("SM+AGC")) anom=true;
    
    if (pID!=0) {
      if (nproc_.nproc>=0) {
	if (nproc_.nproc!=pID)
	  THROW(not_implemented,
		"Only one process class allowed when using MCFM");
      }
      nproc_.nproc=pID;
      chooser_();
      // std::cout << "chosen l1:  " << zcouple_.l1 << std::endl;
      msg_Info()<<"Initialise MCFM with nproc = "<<nproc_.nproc<<"\n";
      return new MCFM_qqb_zgam(pID,pi,fl,anom);
    }
  }
  return NULL;
}
