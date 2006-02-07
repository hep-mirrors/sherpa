#include "Tau_Decay_MEs.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Histogram.H"
#include "Run_Parameter.H"
#include <stdio.h>

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  leptonic decay  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Lepton::Tau_Lepton( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_nulep(-1),
  m_lep(-1)
{
  m_metype = string("Tau_Lepton");
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == kf::e ||
        p_flavs[i].Kfcode() == kf::mu ) { m_lep = i; break; }           // that's the lepton
  }
  // find the corresponding neutrinos
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
    if( p_flavs[i].Kfcode() == p_flavs[m_lep].Kfcode()+1 ) m_nulep = i;
  }
}
 
void Tau_Lepton::SetModelParameters( GeneralModel _md ) 
{ 
  m_GF2 = sqr( _md("GF", rpa.gen.ScalarConstant(string("GF")) ) ); 
  m_cR1 = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL1 = Complex(0.,_md("a",1.)+_md("b",1.));
  m_cR2 = Complex(0.,_md("a2",1.)-_md("b2",1.));
  m_cL2 = Complex(0.,_md("a2",1.)+_md("b2",1.));
}
 
double Tau_Lepton::Using_Hels( const Vec4D * _p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h1=0; h1<4; h1++ ) for( int h2=0; h2<4; h2++ ) {
    int h = (h1<<2)+h2;     // helicity combination (nutau,tau,lep,nulep)
    ret += norm( F.Z( m_nutau, 0, m_lep, m_nulep, h, m_cR1, m_cL1, m_cR2, m_cL2 ) );
  }
  F.Delete();
  return ret*0.25;
} // its value is 100% identical to traces calculation
 
double Tau_Lepton::operator()( const Vec4D *_p )
{
  double T (1.);
  T = Using_Hels(_p);
  return T*m_GF2;
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  1 pion/kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Pseudo::Tau_Pseudo( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion(-1)
{
  m_metype = string("Tau_Pseudo");
  for( int i=1; i<3; i++ ) {
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
    else m_pion = i;                            // that's the pion
  }
  m_pionmode = (p_flavs[m_pion].Kfcode() == kf::pi_plus) ? 1 : 0;
    // 1=pion mode, 0=kaon mode
}
 
void Tau_Pseudo::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vxx2 = sqr( m_pionmode ? 
      _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real()) : 
      _md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real()) );
  m_fxx2 = sqr( m_pionmode ? _md("fpi", 0.0924) : _md("fK", 0.113) );
  m_GF2  = sqr( _md("GF", rpa.gen.ScalarConstant(string("GF")) ) ); 
  m_cR   = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL   = Complex(0.,_md("a",1.)+_md("b",1.));
}
  
double Tau_Pseudo::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {
    ret += norm( F.X(m_nutau,m_pion,0,h,m_cR,m_cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical traces calculation

double Tau_Pseudo::operator()( const Vec4D *_p )
{
  double T(1.);
  double q2 = _p[m_pion].Abs2(); 
  T = Using_Hels(_p); 
  return T*0.5*m_GF2*m_Vxx2*m_fxx2; 
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  2 pion/kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Two_Pion::Tau_Two_Pion( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion_ch(-1),
  m_pion0(-1)
{
  m_metype = string("Tau_TwoPion/TwoKaon");
  for( int i=1; i<4; i++ ) {
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )  m_nutau = i;
    else {
      if( p_flavs[i].Kfcode() == kf::pi || 
          p_flavs[i].Kfcode() == kf::K ||
          p_flavs[i].Kfcode() == kf::K_S ||
          p_flavs[i].Kfcode() == kf::K_L )                       m_pion0 = i;
      else                                              m_pion_ch = i;
    }
  }
  m_pionmode = (p_flavs[m_pion_ch].Kfcode() == kf::pi_plus) ? 1 : 0;
    // 1 = 2 pion mode, 0 = 2 kaon mode
}

void Tau_Two_Pion::SetModelParameters( GeneralModel _md ) 
{ 
  m_m        = p_masses[m_pion_ch];
  m_m2       = m_m*m_m;
  m_running  = int( _md("RUNNING_WIDTH", 1 ) );
  m_ff       = int( _md("FORM_FACTOR", 1 ) );
  m_Vud2     = sqr( _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real() ) );
  m_fxx      = m_pionmode ? _md("fpi", 0.0924) : _md("fK", 0.113 );
  m_GF2      = sqr( _md("GF", rpa.gen.ScalarConstant(string("GF")) ) ); 
  m_CG       = m_pionmode ? 1. : 0.5;   // Clebsch-Gordon
  m_cR       = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL       = Complex(0.,_md("a",1.)+_md("b",1.));

  m_MR       = _md("Mass_rho_770", Flavour(kf::rho_770_plus).PSMass() );
  m_MRR      = _md("Mass_rho_1450", Flavour(kf::rho_1450_plus).PSMass() );
  m_MRRR     = _md("Mass_rho_1700", Flavour(kf::rho_1700_plus).PSMass() );
  m_GR       = _md("Width_rho_770", Flavour(kf::rho_770_plus).Width() );
  m_GRR      = _md("Width_rho_1450", Flavour(kf::rho_1450_plus).Width() );
  m_GRRR     = _md("Width_rho_1700", Flavour(kf::rho_1700_plus).Width() );

  m_MR2      = m_MR*m_MR;
  m_MRR2     = m_MRR*m_MRR;
  m_MRRR2    = m_MRRR*m_MRRR;

  m_m2_pi    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_m2_K     = sqr( Flavour(kf::K_plus).PSMass() );
   
  // coefficients for KS model
  m_beta     = _md("beta", 0. );
  m_delta    = _md("delta", 0. );
  m_gamma    = _md("gamma", 0. );
  m_lambda   = _md("lambda", 1. );

  // coefficients for RChT model
  m_gammaR   = _md("gamma_rho_770", 1. );
  m_gammaRR  = _md("gamma_rho_1450", 1. );
  m_gammaRRR = _md("gamma_rho_1700", 1. );
}
 
// loop function
Complex Tau_Two_Pion::A( double x, double y )
{
  Complex ret(0.,0.); 
  Complex sigma = csqrt(1.-4.*x);
  ret = log(y) + 8.*x - 5./3. + pow(sigma,3.)*log((sigma+1.)/(sigma-1.)); 
  return ret;
}

Complex Tau_Two_Pion::FormFactor( double s )
{

  Complex ret(1.,0.);
  double MG_R   = m_MR*m_GR;            // mass * width of rho
  double MG_RR  = m_MRR*m_GRR;
  double MG_RRR = m_MRRR*m_GRRR;
  if( m_ff == 1 ) {         // Breit-Wigner-rho
    if (m_running) {
      MG_R   = Tools::OffShellMassWidth( s, m_MR2, m_GR, m_m2, 1. );
      MG_RR  = Tools::OffShellMassWidth( s, m_MRR2, m_GRR, m_m2, 1. );
      MG_RRR = Tools::OffShellMassWidth( s, m_MRRR2, m_GRRR, m_m2, 1. );
    }
    Complex BWr   = Tools::BreitWigner( s, m_MR2, MG_R );
    Complex BWrr  = Tools::BreitWigner( s, m_MRR2, MG_RR );
    Complex BWrrr = Tools::BreitWigner( s, m_MRRR2, MG_RRR );
    ret = ( BWr + m_beta*BWrr + m_gamma*BWrrr )/( 1.+m_beta+m_gamma );
  }
  if( m_ff == 2 ) {         // Resonance Chiral Theory
    Complex AA = A( m_m2_pi/s, m_m2_pi/m_MR2 ) + 0.5*A( m_m2_K/s, m_m2_K/m_MR2 );
    double expon = -1.*s/(96.*sqr(M_PI*m_fxx))*AA.real();
    double MG_R, MG_RR, MG_RRR;
    if (m_running) {
      MG_R   = -m_gammaR  *1.*m_MR2  *s/(96.*sqr(M_PI*m_fxx)) * AA.imag();
      MG_RR  = -m_gammaRR *1.*m_MRR2 *s/(96.*sqr(M_PI*m_fxx)) * AA.imag();
      MG_RRR = -m_gammaRRR*1.*m_MRRR2*s/(96.*sqr(M_PI*m_fxx)) * AA.imag();
    }
    Complex BW_1 = Tools::BreitWigner( s, m_MR2, MG_R );
    Complex BW_2 = Tools::BreitWigner( s, m_MRR2, MG_RR );
    Complex BW_3 = Tools::BreitWigner( s, m_MRRR2, MG_RRR );
    ret = (BW_1+m_beta*BW_2+m_gamma*BW_3)/(1.+m_beta+m_gamma) * exp(expon);
  }
  return ret;
}

double Tau_Two_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {                // only nutau and tau have a helicity !
    ret += norm( F.X(m_nutau,m_pion_ch,0,h,m_cR,m_cL)
                 - F.X(m_nutau,m_pion0,0,h,m_cR,m_cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to traces calculation

double Tau_Two_Pion::operator()( const Vec4D *_p )
{
  double T = Using_Hels(_p);
  double q2 = (_p[m_pion_ch] + _p[m_pion0] ).Abs2();
  Complex FF = FormFactor(q2);
  return T*0.5*m_GF2*m_Vud2*norm(FF)*m_CG;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   pion-kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Pion_Kaon::Tau_Pion_Kaon( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion(-1),
  m_kaon(-1)
{
  m_metype = string("Tau_PionKaon");
  m_ms[0] = sqr(p_flavs[0].PSMass());
  for( int i=1; i<4; i++ ) {
    m_ms[i] = sqr(p_flavs[i].PSMass());
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )          m_nutau = i;    // neutrino
    else {
      if( p_flavs[i].Kfcode() == kf::pi || 
          p_flavs[i].Kfcode() == kf::pi_plus )   m_pion = i;     // pion
      else if (m_kaon<0)                                        m_kaon = i;     // kaon
    }
  }
  m_chpionmode = (p_flavs[m_pion].Kfcode() == kf::pi_plus) ? 1 : 0;
}

void Tau_Pion_Kaon::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vus2     = sqr( _md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real() ) );
  m_fpi2     = sqr( _md("fpi", 0.0924 ) );
  m_GF2      = sqr( _md("GF", rpa.gen.ScalarConstant(string("GF")) ) ); 
  m_cR       = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL       = Complex(0.,_md("a",1.)+_md("b",1.));
  m_Delta_KP = m_ms[m_kaon] - m_ms[m_pion];
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : /*p_ff = new RChT();
             break;        */       // use RChT on own risk: not tested sufficiently
    case 1 : p_ff = new KS();
             break;
  }
  p_ff->SetModelParameters( _md );
  p_ff->SetMasses2( m_ms[m_pion], m_ms[m_kaon], sqr(Flavour(kf::eta).PSMass()) );
}

// Resonance Chiral Theory

void Tau_Pion_Kaon::RChT::SetModelParameters( GeneralModel _md ) 
{
  m_MK2      = sqr( _md("Mass_K*_892", Flavour(kf::K_star_892_plus).PSMass()) );
  m_GK       = _md("Width_K*_892", Flavour(kf::K_star_892_plus).Width());
  m_MK02     = sqr( _md("Mass_K*0_892", Flavour(kf::K_star_892).PSMass()) );
  m_GK0      = _md("Width_K*0_892", Flavour(kf::K_star_892).Width());
  m_fpi2     = sqr( _md("fpi", 0.0924) );
  m_renorm2  = sqr( _md("renorm",_md("Mass_rho_770",Flavour(kf::rho_770_plus).PSMass())));
}

void Tau_Pion_Kaon::RChT::SetMasses2( double _mPi2, double _mK2, double _mEta2 )
{
  m_mPi2     = _mPi2;
  m_mK2      = _mK2;
  m_mEta2    = _mEta2;
  m_mPi      = sqrt(m_mPi2);
  m_mK       = sqrt(m_mK2);
  m_mEta     = sqrt(m_mEta2);
  m_Sigma_KP = m_mK2+m_mPi2;
  m_Delta_KP = m_mK2-m_mPi2;
}
 
double Tau_Pion_Kaon::RChT::MassWidthVector( double s )
{
  double ret (0.);
  if (s>sqr(m_mK+m_mPi))  ret += pow( Tools::Lambda(s,m_mK2,m_mPi2), 1.5 );
  if (s>sqr(m_mK+m_mEta)) ret += pow( Tools::Lambda(s,m_mK2,m_mEta2), 1.5 );
  ret *= m_MK2 /( 128.*M_PI*m_fpi2*sqr(s) );
}

double Tau_Pion_Kaon::RChT::MassWidthScalar( double s )
{
  double G_s = m_GK0*m_MK02/s
    * pow( Tools::Lambda(s,m_mK2,m_mPi2)/Tools::Lambda(m_MK02,m_mK2,m_mPi2), 
           0.5 );
  return sqrt(m_MK02) * G_s;
}

Complex Tau_Pion_Kaon::RChT::JBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double nu = sqrt( sqr(s) + sqr(MP2) + sqr(MQ2) - 2.*s*Sigma - 2.*s*MP2*MQ2 );
  double  J = 2. + Delta/s*log(MQ2/MP2) 
    - Sigma/Delta*log(MQ2/MP2) 
    - nu/s*log( (sqr(s+nu)-sqr(Delta))/(sqr(s-nu)-sqr(Delta)) ); 
  return J/(32.*sqr(M_PI));
}

Complex Tau_Pion_Kaon::RChT::JBarBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double Jp_0 = ( Sigma/sqr(Delta) + 2.*MP2*MQ2/pow(Delta,3)*log(MQ2/MP2) )/( 32.*sqr(M_PI) );
  return JBar(s,MP2,MQ2,Sigma,Delta) - s*Jp_0;
}

Complex Tau_Pion_Kaon::RChT::Mr( double s, double MP2, double MQ2 )
{
  double Sigma = MP2 + MQ2;
  double Delta = MP2 - MQ2;
  Complex Jb   = JBar(s,MP2,MQ2,Sigma,Delta);
  Complex Jbb  = JBarBar(s,MP2,MQ2,Sigma,Delta);
  double mu2   = m_renorm2;
  double k     = ( MP2*log(MP2/mu2) - MQ2*log(MQ2/mu2) ) / ( Delta*32.*sqr(M_PI) );
  Complex M    = (s-2.*Sigma)*Jb/(12.*s) + sqr(Delta)/(3.*sqr(s))*Jbb - k/6. + 1./(288.*sqr(M_PI));
  return M;
}

Complex Tau_Pion_Kaon::RChT::L( double s, double MP2, double MQ2 )
{
  double Sigma = MP2 + MQ2;
  double Delta = MP2 - MQ2;
  Complex Jb   = JBar(s,MP2,MQ2,Sigma,Delta);
  return( sqr(Delta)/(4.*s)*Jb );
}

Complex Tau_Pion_Kaon::RChT::VectorFormFactor( double s )
{
  Complex ret(1.,0.);
  double    MG_K = MassWidthVector(s);
  Complex     BW = Tools::BreitWigner( s, m_MK2, MG_K );
  Complex M_part = Mr(s,m_mK2,m_mPi2) - Mr(s,m_mK2,m_mEta2);
  Complex L_part = L(s,m_mK2,m_mPi2) - L(s,m_mK2,m_mEta2);
  double   expon = 3./2./m_fpi2 *( s*M_part.real() - L_part.real() );
  ret = BW * exp(expon);
  return ret;
}

Complex Tau_Pion_Kaon::RChT::ScalarFormFactor( double s )
{
  Complex ret(1.,0.);
  double  MG_K0 = MassWidthScalar(s);
  Complex    BW = Tools::BreitWigner( s, m_MK02, MG_K0 );
  double    cd2 = sqr(0.032);
  Complex    F4 = 1./(8.*m_fpi2)
    * ( 5.*s - 2.*m_Sigma_KP - 3.*sqr(m_Delta_KP)/s )
    * JBar(s,m_mK2,m_mPi2,m_mK2+m_mPi2,m_mK2-m_mPi2)
    + 1./(24.*m_fpi2)
    * ( 3.*s - 2.*m_Sigma_KP - sqr(m_Delta_KP)/s )
    * JBar(s,m_mK2,m_mEta2,m_mK2+m_mEta2,m_mK2-m_mEta2);
  double  inter = 1. - (1.-m_fpi2/4./cd2)*m_Sigma_KP/m_MK02;
  Complex  expon = Complex( F4.real(), F4.imag()/(1.+sqr(F4.imag())) );
  ret = BW * inter * exp(expon);
  return ret;
}

// Kuehn Santamaria Model

void Tau_Pion_Kaon::KS::SetModelParameters( GeneralModel _md ) 
{
  m_MR      = _md("Mass_K*_892", Flavour(kf::K_star_892_plus).PSMass());
  m_MRR     = _md("Mass_K*_1410", Flavour(kf::K_star_1410).PSMass());

  m_GR      = _md("Width_K*_892", Flavour(kf::K_star_892_plus).Width());
  m_GRR     = _md("Width_K*_1410", Flavour(kf::K_star_1410_plus).Width());

  m_MR2     = sqr(m_MR);
  m_MRR2    = sqr(m_MRR);

  m_m2      = sqr( Flavour(kf::pi_plus).Mass() );
  m_mK2     = sqr( Flavour(kf::K_plus).Mass() );

  m_beta    = _md("beta", 0.);
  m_running = int( _md("RUNNING_WIDTH", 1 ) );
}
 

Complex Tau_Pion_Kaon::KS::VectorFormFactor( double s )
{
  Complex ret(1.,0.);
  double MG_R   = m_MR*m_GR;            // mass * width of K*
  double MG_RR  = m_MRR*m_GRR;
  if (m_running) {
    MG_R   = Tools::OffShellMassWidth( s, m_MR2, m_GR, m_m2, m_mK2, 1. );
    MG_RR  = Tools::OffShellMassWidth( s, m_MRR2, m_GRR, m_m2, m_mK2, 1. );
  }
  Complex BWr   = Tools::BreitWigner( s, m_MR2, MG_R );
  Complex BWrr  = Tools::BreitWigner( s, m_MRR2, MG_RR );
  ret = ( BWr + m_beta*BWrr )/( 1.+m_beta );
  return ret;
}

Complex Tau_Pion_Kaon::KS::ScalarFormFactor( double s )
{
  Complex ret(0.,0.);
  return ret;
}

// general framework

double Tau_Pion_Kaon::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  double q2 = (_p[m_pion]+_p[m_kaon]).Abs2();
  Complex FS = p_ff->ScalarFormFactor(q2);
  Complex FV = p_ff->VectorFormFactor(q2);
  Complex termK = m_Delta_KP/q2*(FS-FV)+FV;
  Complex termP = m_Delta_KP/q2*(FS-FV)-FV;
  for( int h=0; h<4; h++ ) {                // only nutau and tau have a helicity !
    ret += norm( 
          F.X(m_nutau,m_kaon,0,h,m_cR,m_cL) * termK
        + F.X(m_nutau,m_pion,0,h,m_cR,m_cL) * termP
        );
  }
  F.Delete();
  return ret/4.;
} 

double Tau_Pion_Kaon::operator()( const Vec4D *_p )
{
  double T = Using_Hels(_p);
  return T*0.5*m_GF2*m_Vus2;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  3 pseudo mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
Tau_Three_Pseudo::Tau_Three_Pseudo( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pseudo_3(-1),
  m_pseudo_1(-1),
  m_pseudo_2(-1)
{
  m_metype = string("Tau_ThreePseudo");
  int nPion_0(0), nPion_ch(0), nKaon_0(0), nKaon_ch(0);
  m_ms[0] = sqr( p_flavs[0].PSMass() );
  // count number of pions, kaons and calc. mass^2
  for( int i=1; i<5; i++ ) {
    if( p_flavs[i].Kfcode() == kf::pi_plus )   nPion_ch++;
    if( p_flavs[i].Kfcode() == kf::pi )        nPion_0++;
    if( p_flavs[i].Kfcode() == kf::K_plus )    nKaon_ch++;
    if( p_flavs[i].Kfcode() == kf::K ||                   
        p_flavs[i].Kfcode() == kf::K_L ||                     
        p_flavs[i].Kfcode() == kf::K_S )       nKaon_0++;
    m_ms[i] = sqr( p_flavs[i].PSMass() );
  }
  // sanity check
  if (nPion_ch+nPion_0+nKaon_ch+nKaon_0 != 3) {
    msg.Error()<<"ERROR in HADRONS::Tau_Three_Pseudo constructor\n"
               <<"     number of three outgoing pseudoscalars != 3 ?!.\n"
               <<"     Don't know, what to do. Will abort."<<endl;
    abort();           
  }
  // define mode number
  m_mode = nPion_ch*1000 + nPion_0*100 + nKaon_ch*10 + nKaon_0;
  switch( m_mode ) {
    case 1200 : /* pi0 pi0 pi- mode */
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
                  else {
                    if( p_flavs[i].Kfcode() == kf::pi_plus ) m_pseudo_3 = i;
                    else if( m_pseudo_1<0 ) m_pseudo_1 = i;
                    else m_pseudo_2 = i;
                  }
                }
                break;
    case  210 : /* pi0 pi0 K- mode */
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
                  else {
                    if( p_flavs[i].Kfcode() == kf::K_plus ) m_pseudo_3 = i;
                    else if( m_pseudo_1<0 ) m_pseudo_1 = i;
                    else m_pseudo_2 = i;
                  }
                }
                break;
    case   30 : /* K- K- K+ */
    case 1020 : /* K- pi- K+ */
    case 2010 : /* K- pi- pi+ */
    case 3000 : /* pi- pi- pi+ mode */
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
                  else {
                    if ( p_flavs[i].IsAnti() == p_flavs[0].IsAnti() ) m_pseudo_3 = i;
                    else if( m_pseudo_1<0 ) m_pseudo_1 = i;
                    else m_pseudo_2 = i;
                  }
                }
                break;
    case 1002 : /* KS pi- KL */ 
                /* KS KS pi- */ 
                /* pi- KL KL */
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )    m_nutau = i;
                  else {
                    if( p_flavs[i].Kfcode() == kf::K_S ) {
                      if (m_pseudo_1 < 0)                               m_pseudo_1 = i;
                      else                                              m_pseudo_3 = i;
                    }
                    if( p_flavs[i].Kfcode() == kf::pi_plus )            m_pseudo_2 = i;
                    if( p_flavs[i].Kfcode() == kf::K_L ) {
                      if (m_pseudo_3<0)                                 m_pseudo_3 = i;
                      else                                              m_pseudo_1 = i;
                    }
                  }
                }
                break;
    case  111 : /* K- pi0 K0 */ /*K0 = KS, KL*/
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )    m_nutau = i;
                  else {
                    if( p_flavs[i].Kfcode() == kf::K ||
                        p_flavs[i].Kfcode() == kf::K_L ||
                        p_flavs[i].Kfcode() == kf::K_S )                m_pseudo_3 = i;
                    if( p_flavs[i].Kfcode() == kf::pi )                 m_pseudo_2 = i;
                    if( p_flavs[i].Kfcode() == kf::K_plus )             m_pseudo_1 = i;
                  }
                }
                break;
    case 1101 : /* pi- K0 pi0 */ /*K0 = KS, KL*/
                for( int i=1; i<5; i++ ) {
                  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )    m_nutau = i;
                  else {
                    if( p_flavs[i].Kfcode() == kf::pi )                 m_pseudo_3 = i;
                    if( p_flavs[i].Kfcode() == kf::K ||
                        p_flavs[i].Kfcode() == kf::K_S ||
                        p_flavs[i].Kfcode() == kf::K_L )                m_pseudo_2 = i;
                    if( p_flavs[i].Kfcode() == kf::pi_plus )            m_pseudo_1 = i;
                  }
                }
                break;
  }
  // get internal number of outgoing particles (in KS order)
  m_part[m_pseudo_1-1] = 1;
  m_part[m_pseudo_2-1] = 2;
  m_part[m_pseudo_3-1] = 3;
  m_part[m_nutau-1] = 4;
}
 
void Tau_Three_Pseudo::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vud     = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real() );
  m_Vus     = _md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real() );
  m_fpi2    = sqr(_md("fpi", 0.0924 ));
  m_GF2     = sqr( _md("GF", rpa.gen.ScalarConstant(string("GF")) )); 
  m_cR      = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL      = Complex(0.,_md("a",1.)+_md("b",1.));
  SetA123();
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : p_ff = new RChT();
             break;
    case 1 : p_ff = new KS();
             break;
  }
  p_ff->SetOutgoingMasses2( m_ms[m_pseudo_1], m_ms[m_pseudo_2], m_ms[m_pseudo_3]  );
  p_ff->SetPath( m_path );
  p_ff->SetMode( m_mode );
  p_ff->SetInternalNumbers( m_part );
  p_ff->SetModelParameters( _md );
}
 
void Tau_Three_Pseudo::SetA123()
{
  switch( m_mode ) {
    case 1200 : /* pi0 pi0 pi- mode */
                m_A123 = m_Vud;
                break;
    case   30 : /* K- K- K+ */
                m_A123 = m_Vus;
                break;
    case 2010 : /* K- pi- pi+ */
                m_A123 = -0.5*m_Vus;
                break;
    case 1020 : /* K- pi- K+ */
                m_A123 = -0.5*m_Vud;
                break;
    case 3000 : /* pi- pi- pi+ mode */
                m_A123 = m_Vud;
                break;
    case 1002 : /* K0 pi- K0b */
                m_A123 = -0.5*m_Vud;
                break;
    case  111 : /* K- pi0 K0 */
                m_A123 = 1.5*m_Vud*SQRT_05;
                break;
    case  210 : /* pi0 pi0 K- mode */
                m_A123 = m_Vus/4.;
                break;
    case 1101 : /* pi- K0b pi0 */
                m_A123 = 1.5*m_Vus*SQRT_05;
                break;
    default   : msg.Error()<<"Warning in HADRONS::Tau_Decay_MEs.C in Tau_Three_Pseudo::SetA123() :"
                           <<"     Obviously this three pseudoscalar channel (code "<<m_mode<<")\n"
                           <<"     doesn't have a global A123. Maybe it is not implemented yet.\n"
                           <<"     Take A123=1., will continue and hope for the best."<<endl;
                m_A123 = 1.;
                break;
  }
  m_global = sqr(2.*m_A123/3.);
}

// Parameterisation
// DUMM, PICH, PORTOLES hep-ph/0312183 

void Tau_Three_Pseudo::RChT::SetModelParameters( GeneralModel _md ) 
{ 
  // it is this part's job to map all the vector resonances correctly !
  int a[2], b[2], c[2];                         // help variables to do things right
  double MV[2], GV[2];                          // store information from DC file
  if (!m_twoident) {                            // not two identical vector resonances
    for( int i=0; i<2; i++ ) {          // i = 0,1
      a[i] = (i==0) ? ParticleNo(int(_md("vector1_i", 1))-1) 
                    : ParticleNo(int(_md("vector2_i", 2))-1);
      b[i] = (i==0) ? ParticleNo(int(_md("vector1_j", 3))-1) 
                    : ParticleNo(int(_md("vector2_j", 3))-1);
      if ( a[i]!=3 && b[i]!=3 ) {
        ATOOLS::msg.Error()<<ATOOLS::om::red
          <<"ERROR in Tau_Three_Pseudo::RChT::SetModelParameters\n"
          <<"     Resonances aren't set correctly.\n"
          <<"     Make sure you have the right settings under \"Resonances\" vector"<<i+1<<" -> _ _.\n"
          <<"     mode number : "<<m_mode<<" (= #pi #pi0 #K #K0 in out state)\n"
          <<"     Don't know what to, will abort."
          <<ATOOLS::om::reset<<std::endl;
        abort();
      }
      else c[i] = (a[i]==3)? b[i] : a[i];
    }
    // masses and widths of vector resonances (v and v')
    MV[0] = _md("Mass_vector_1", Flavour(kf::rho_770_plus).PSMass() );
    GV[0] = _md("Width_vector_1", Flavour(kf::rho_770_plus).Width() );
    MV[1] = _md("Mass_vector_2", Flavour(kf::rho_770_plus).PSMass() );
    GV[1] = _md("Width_vector_2", Flavour(kf::rho_770_plus).Width() );
  }
  else {
    c[0] = c[1] = 1;
    // masses and widths of vector resonances (v and v')
    MV[0] = _md("Mass_vector", Flavour(kf::rho_770_plus).PSMass() );
    GV[0] = _md("Width_vector", Flavour(kf::rho_770_plus).Width() );
    MV[1] = MV[0];
    GV[1] = GV[0];
  }
   
  // set correct parameters (using corresponding settings under "Resonances" in DC file)
   
  m_MA       = _md("Mass_axial", Flavour(kf::a_1_1260_plus).PSMass());            // mass of axial resonance
  m_MA2      = sqr(m_MA);                           // mass^2 of axial resonance
  m_msV[0]   = sqr( MV[c[0]-1] );                   // mass^2 of vector resonance 13
  m_msV[1]   = sqr( MV[c[1]-1] );                   // mass^2 of vector resonance 23
   
  m_GA_at_MA2 = _md("Width_axial", Flavour(kf::a_1_1260_plus).Width());         // on-shell axial width
  m_widthV[0] = GV[c[0]-1];                         // on-shell vector 13 width
  m_widthV[1] = GV[c[1]-1];                         // on-shell vector 23 width

  m_fpi2   = sqr(_md("fpi", 0.0924));
  m_l0     = _md("lambda0", 1.);                    // fit parameter lambda0
  m_gammaR = _md("gamma_rho_770", 1.);              // global factor for rho width
  m_m      = Flavour( kf::pi_plus ).PSMass();         // pion mass
  m_m2     = sqr(m_m);                              // pion mass^2
  m_mK2    = sqr( Flavour( kf::K_plus ).PSMass() );   // Kaon mass^2
  m_exp_alpha = _md("exp_alpha", 2.45);             // exponent in off-shell GA
  m_l1     = _md("lambda1", 0.5);                   // fit parameter        
  m_l2     = _md("lambda2", 0.);                    // fit parameter
  m_lsum   = m_l1 + m_l2; 
  // constraints due to short-distance behaviour
  m_FV2    = 2*m_fpi2;                              // vector coupling
  m_FV     = sqrt(m_FV2);
  m_GV     = m_FV/2.;
  m_FA2    = m_FV2 - m_fpi2;                        // axial coupling
  m_FA     = sqrt(m_FA2);
  
  m_running  = int( _md("RUNNING_WIDTH", 3 ) );     // running width

  if( m_running&1 ) {                               // if axial width running
    p_phi=CreatePhiHistogram();
    m_Phi_at_MA2 = Phi( m_MA2 );
  }
  else 
    m_Phi_at_MA2 = 1.;
}
 
void Tau_Three_Pseudo::RChT::SetMode( int m ) 
{ 
  m_mode = m;
  m_twoident = false;
  if (m_mode==1200 ||
      m_mode==210 ||
      m_mode==30 ||
      m_mode==3000 ||
      m_mode==12) m_twoident = true;
}

Histogram * Tau_Three_Pseudo::RChT::CreatePhiHistogram()
{
  // create file name
  char fn[256];
  sprintf(fn, "%sPhaseSpaceFunctions/PhiQ2_MV13=%.3f_MV23=%.3f_fpi=%.4f_gR=%.2f_GV13=%.3f_GV23=%.3f_run=%i.dat",
      m_path.c_str(),
      sqrt(m_msV[0]),
      sqrt(m_msV[1]),
      sqrt(m_fpi2),
      m_gammaR,
      m_widthV[0],
      m_widthV[1],
      (m_running&2));

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {                          // if file exists
    // read table and create histogram
    msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::RChT::CreatePhiHistogram : \n"
             <<"     Read phi(q2) from "<<fn<<"."<<endl;
    myHist = new Histogram( fn );
  }
  else {                            // if file does not exist
    // create histogram (i.e. table of values)
    msg.Out()<<"Create necessary phase space function for choosen parameters.\n"
             <<"This may take some time. Please wait..."<<endl;
    msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::RChT::CreatePhiHistogram : \n"
             <<"     Create phi(q2) in "<<fn<<"."<<endl;
    double low (0.),
    up  (3.2);
    int nbins  (50);
    double step = (up-low)/nbins;
    myHist = new Histogram( 0, low, up+step, nbins+1 );
    double q2  (low),
    phi (0.);
    while( q2<=up+step ) {
      phi = IntegralPhi(q2);                // get phi value
      myHist->Insert( q2, phi );            // insert into histogram
      q2 += step;       
    }
    myHist->Output(string(fn));
  }
  return myHist;
}

double Tau_Three_Pseudo::RChT::IntegralPhi( double Q2 )
{
  int Ns=500, Nt=500;                   // number of subintervals
  double sum (0.);
  double s_max = sqr( sqrt(Q2)-m_m );
  double s_min = 4.*m_m2;
  double ds = (s_max-s_min)/Ns;
  double t_max (0.);
  double t_min (0.);
  double dt = (0.);
  double s (s_min), t, u;
  double V12, V22, V1V2;
  
  while ( s<s_max ) {
    t_max = (   sqr( Q2-Mass2(1-1)-Mass2(2-1)+Mass2(3-1) ) 
              - sqr( Sqrt_Lambda(Q2, s, Mass2(2-1)) - Sqrt_Lambda(s,Mass2(1-1),Mass2(3-1)) ) 
            )/(4.*s);
    t_min = (   sqr( Q2-Mass2(1-1)-Mass2(2-1)+Mass2(3-1) ) 
              - sqr( Sqrt_Lambda(Q2, s, Mass2(2-1)) + Sqrt_Lambda(s,Mass2(1-1),Mass2(3-1)) ) 
            )/(4.*s);
    dt = (t_max-t_min)/Nt;
    t = t_min;
    while ( t<t_max ) {
      u    = Q2 - s - t + Mass2(1-1)+Mass2(2-1)+Mass2(3-1);
      V12  = 2.*(Mass2(1-1)+Mass2(2-1)) - s - sqr(u-t+Mass2(1-1)-Mass2(3-1))/(4.*Q2);
      V22  = 2.*(Mass2(1-1)+Mass2(2-1)) - t - sqr(u-s+Mass2(2-1)-Mass2(3-1))/(4.*Q2);
      V1V2 = (u-s-t+4.*Mass2(3-1))/2. - (u-t+Mass2(1-1)-Mass2(3-1))*(u-s+Mass2(2-1)-Mass2(3-1))/(4.*Q2);
      Complex BW_s = Tools::BreitWigner( s, m_msV[0], MassWidthVector(0,s) );
      Complex BW_t = Tools::BreitWigner( t, m_msV[1], MassWidthVector(1,t) );
      sum += ( V12 * norm(BW_s) + V22*norm(BW_t) + 2.*V1V2*real( BW_s*conj(BW_t)) )*ds*dt;
      t += dt;
    }
    s += ds;
  }
  return sum*Q2;
}

double Tau_Three_Pseudo::RChT::Phi( double x )
{ 
  double val;   
  p_phi->Extrapolate( x+p_phi->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  double p = val;
  return p;
}

double Tau_Three_Pseudo::RChT::MassWidthVector( int a, double s )
{ 
  if(m_running & 2) {
    switch(m_mode) {
      case 3000:
      case 1200: return MassWidthVector( s );
      case 1020: return (a==0)? MassWidthVector(s) : sqrt(m_msV[1])*m_widthV[1];
    }
    // default:
    msg.Error()<<"Warning: this form factor (RChT) \n"
      <<"     hasn't been implemented yet. Please use KS model."<<endl;
    return  0.;
  }
  return( sqrt(m_msV[a])*m_widthV[a] );
}

double Tau_Three_Pseudo::RChT::MassWidthVector( double s )
{ 
  // make sure that this belongs only to rho !!!!!
  double MVGV (0.);
    if( s>4.*m_m2 )  MVGV += pow( 1.-4.*m_m2/s, 1.5 );
    if( s>4.*m_mK2 ) MVGV += pow( 1.-4.*m_mK2/s, 1.5 ) / 2.;
    MVGV *= m_gammaR*m_msV[0]*s/(96.*M_PI*m_fpi2); 
  return MVGV;
}

double Tau_Three_Pseudo::RChT::MassWidthAxial( double Q2 )
{
  if( m_running & 1 )
    return(  m_MA * m_GA_at_MA2 * Phi(Q2) / m_Phi_at_MA2 
           * pow( m_MA2/Q2, m_exp_alpha ) );
  return m_MA*m_GA_at_MA2;
}
 
Complex Tau_Three_Pseudo::RChT::FormFactor( int j, double Q2, double s, double t )
{
  switch( m_mode ) {
    case 1200:
    case 3000: { // 3pion mode
                 if (j==1 || j==2) {        // axial contributions
                   double u = Q2-s-t+Mass2(1-1)+Mass2(2-1)+Mass2(3-1);
                   double x = (j==1)? s : t;
                   double y = (j==1)? t : s;
                   double MVGV_x = MassWidthVector(x);
                   double MVGV_y = MassWidthVector(y);
                   double MAGA = MassWidthAxial(Q2);
                   double F_Q2_x = x/2./Q2 - m_l0*m_m2/Q2;
                   double F_Q2_y = y/2./Q2 - m_l0*m_m2/Q2;
                   double MV2 = m_msV[0];
                   Complex alpha = 1. - 3./2.*x/Complex(x-MV2, MVGV_x);
                   Complex beta =  -3./2.*x/Complex(x-MV2, MVGV_x)
                     + F_Q2_x*(2.*Q2+x-u)/Complex(x-MV2, MVGV_x)
                     + F_Q2_y*(u-x)/Complex(y-MV2, MVGV_y);
                   return SQRT_05*( alpha - Q2/Complex(Q2-m_MA2,MAGA)*beta );
                   // the (global) factor of SQRT_05 was added in order to obtain a good expression for BR
                 }
                 else {                     // pseudoscalar, vector
                   return Complex(0.,0.);
                 }
               }
    case 1020: { // K- pi- K+ mode
                 double u = Q2-s-t+Mass2(1-1)+Mass2(2-1)+Mass2(3-1);
                 switch(j) {
                   case 1:
                   case 2: {                // axial contributions
                             double x = (j==1)? s : t;
                             double y = (j==1)? t : s;
                             int    a = j-1;                          // resonance assoc. with x
                             int    b = (a==1)? 0 : 1;                 // resonance assoc. with y
                             Complex TAchi = -sqrt(2.)/3.;
                             Complex TA1r  = -sqrt(2.)/6. * m_FV*m_GV/m_fpi2 * ( 
                                  1./Complex(m_msV[a]-x,MassWidthVector(a,x))*
                                    ( 3.*x + (1.-2.*m_GV/m_FV)*(2.*Q2-2.*x-u+Mass2(a)-Mass2(b)) )
                                 +1./Complex(m_msV[b]-y,MassWidthVector(b,y))*
                                    ( 2.*(Mass2(b)-Mass2(2)) + (1.-2.*m_GV/m_FV)*(u-x+Mass2(2)-Mass2(b)) )    
                                 );
                             Complex TA2r  = 2./3. * m_FA*m_GV/m_fpi2 * Q2/Complex(m_MA2-Q2,MassWidthAxial(Q2)) * (
                                  1./Complex(m_msV[a]-x,MassWidthVector(a,x))*
                                    ( m_lsum*(-3.*x+Mass2(a)-Mass2(2)) + FFunc(x,Mass2(b),Q2)*(2.*Q2+x-u+Mass2(2)-Mass2(b)) )
                                 +1./Complex(m_msV[b]-y,MassWidthVector(b,y))*
                                    ( 2.*m_lsum*(Mass2(2)+Mass2(b)) + FFunc(y,Mass2(a),Q2)*(u-x+Mass2(b)-Mass2(2)) )
                                 );
                             return TAchi + TA1r + TA2r;
                           }
                   case 3: {                // pseudoscalar contribution
                             Complex TAchi = +sqrt(2.)/3.*3./2.* Tools::BreitWigner(Q2,Mass2(1),0.) *
                                                ( 1. + (Mass2(2)-u)/Q2 );
                             Complex TA1r  = -sqrt(2.)/6.*m_FV*m_GV * 3.*m_GV/m_FV * Tools::BreitWigner(Q2,Mass2(1),0.)/Q2 *
                                                (   s*(t-u)/Complex(m_msV[0]-s,MassWidthVector(0,s))
                                                  + (t*(s-u)+(Q2-Mass2(0))*(Mass2(1)-Mass2(2)))/Complex(m_msV[1]-t,MassWidthVector(1,t)) 
                                                );
                             return TAchi + TA1r;
                           }
                   case 4: {                // vector contribution
                             return Complex(0.,0.);
                           }
                 }
               }
  }
  // default:
  msg.Error()<<"Warning: this form factor (RChT) \n"
    <<"     hasn't been implemented yet. Please use KS model."
    <<"     Decay mode: "<<m_mode<<endl;
  return  (j==3)? Complex(0.,0.) : Complex(1.,0.); 
}

double Tau_Three_Pseudo::RChT::FFunc( double a, double b, double c)
{
  return m_l2 + m_l1*a/c - m_l0*b/c;
}

// Parameterisation
// DECKER, FINKEMEIER, MIRKES hep-ph/9310270

void Tau_Three_Pseudo::KS::SetModelParameters( GeneralModel _md ) 
{
  // it is this part's job to map all the vector resonances correctly !
  int a[2], b[2], c[2];                         // help variables to do things right
  double MV[2], GV[2], Mv[2], Gv[2], Beta[2];   // store information from DC file
  if (!m_twoident) {                            // not two identical vector resonances
    for( int i=0; i<2; i++ ) {          // i = 0,1
      a[i] = (i==0) ? ParticleNo(int(_md("vector1_i", 1))-1) 
                    : ParticleNo(int(_md("vector2_i", 2))-1);
      b[i] = (i==0) ? ParticleNo(int(_md("vector1_j", 3))-1) 
                    : ParticleNo(int(_md("vector2_j", 3))-1);
      if ( a[i]!=3 && b[i]!=3 ) {
        ATOOLS::msg.Error()<<ATOOLS::om::red
          <<"ERROR in Tau_Three_Pseudo::KS::SetModelParameters\n"
          <<"     Resonances aren't set correctly.\n"
          <<"     Make sure you have the right settings under \"Resonances\" vector"<<i+1<<" -> _ _.\n"
          <<"     mode number : "<<m_mode<<" (= #pi #pi0 #K #K0 in out state)\n"
          <<"     Don't know what to, will abort."
          <<ATOOLS::om::reset<<std::endl;
        abort();
      }
      else c[i] = (a[i]==3)? b[i] : a[i];
    }
    // masses and widths of vector resonances (v and v')
    MV[0] = _md("Mass_vector_1", Flavour(kf::rho_770_plus).PSMass() );
    GV[0] = _md("Width_vector_1", Flavour(kf::rho_770_plus).Width() );
    MV[1] = _md("Mass_vector_2", Flavour(kf::rho_770_plus).PSMass() );
    GV[1] = _md("Width_vector_2", Flavour(kf::rho_770_plus).Width() );
    Mv[0] = _md("Mass_vector'_1", Flavour(kf::rho_1450_plus).PSMass() );
    Gv[0] = _md("Width_vector'_1", Flavour(kf::rho_1450_plus).Width());
    Mv[1] = _md("Mass_vector'_2", Flavour(kf::rho_1450_plus).PSMass() );
    Gv[1] = _md("Width_vector'_2", Flavour(kf::rho_1450_plus).Width() );
    // relative strength
    Beta[0] = _md("beta_1", 0. );
    Beta[1] = _md("beta_2", 0. );
  }
  else {
    c[0] = c[1] = 1;
    // masses and widths of vector resonances (v and v')
    MV[0] = _md("Mass_vector", Flavour(kf::rho_770_plus).PSMass() );
    GV[0] = _md("Width_vector", Flavour(kf::rho_770_plus).Width() );
    MV[1] = MV[0];
    GV[1] = GV[0];
    Mv[0] = _md("Mass_vector'", Flavour(kf::rho_1450_plus).PSMass() );
    Gv[0] = _md("Width_vector'", Flavour(kf::rho_1450_plus).Width() );
    Mv[1] = Mv[0];
    Gv[1] = Gv[0];
    // relative strength
    Beta[0] = _md("beta", 0. );
    Beta[1] = Beta[0];
  }
   
  // set correct parameters (using corresponding settings under "Resonances" in DC file)
   
  m_MA       = _md("Mass_axial", Flavour(kf::a_1_1260_plus).PSMass());     // mass of axial resonance
  m_MAA2     = _md("Mass_axial'", Flavour(kf::a_1_1260_plus).Width());   // mass of axial resonance'
  m_MA2      = sqr( m_MA );                         // mass^2 of axial resonance
  m_MAA2     = sqr( m_MAA );                        // mass^2 of axial resonance'
  m_msV[0]   = sqr( MV[c[0]-1] );                   // mass^2 of vector resonance 13
  m_msV[1]   = sqr( MV[c[1]-1] );                   // mass^2 of vector resonance 23
  m_msv[0]   = sqr( Mv[c[0]-1] );                   // mass^2 of vector resonance' 13
  m_msv[1]   = sqr( Mv[c[1]-1] );                   // mass^2 of vector resonance' 23
  m_Beta[0]  = Beta[c[0]-1];                        // weight factor for vector resonance' 13
  m_Beta[1]  = Beta[c[1]-1];                        // weight factor for vector resonance' 23
  m_alpha    = _md("alpha", 0.);                    // weight factor for axial resonance'
   
  m_GA        = _md("Width_axial",Flavour(kf::a_1_1260_plus).Width()  );         // on-shell axial width
  m_GAA       = _md("Width_axial'", Flavour(kf::a_1_1260_plus).Width() );        // on-shell axial' width
  m_widthV[0] = GV[c[0]-1];                         // on-shell vector 13 width
  m_widthV[1] = GV[c[1]-1];                         // on-shell vector 23 width
  m_widthv[0] = Gv[c[0]-1];                         // on-shell vector' 13 width
  m_widthv[1] = Gv[c[1]-1];                         // on-shell vector' 23 width

  m_exp_alpha = _md("exp_alpha", 2.45 );            // exponent in off shell GA

  m_running  = int( _md("RUNNING_WIDTH", 3 ) );     // running width

  if( m_running&1 ) {                               // only if axial width running
    p_G = CreateGHistogram();                             
    m_G_at_MA2  = G(m_MA2);                         // G value at axial mass^2
    m_G_at_MAA2 = G(m_MA2);                         // G value at axial' mass^2
  }
  else {
    m_G_at_MA2  = 1.;
    m_G_at_MAA2 = 1.;
  }
} 
 
void Tau_Three_Pseudo::KS::SetMode( int m ) 
{ 
  m_mode = m;
  m_twoident = false;
  if (m_mode==1200 ||
      m_mode==210 ||
      m_mode==30 ||
      m_mode==3000 ||
      m_mode==12) m_twoident = true;
  switch( m_mode ) {
    case 1200 : /* pi0 pi0 pi- mode */
      m_X123  = Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      break;
    case   30 : /* K- K- K+ */
      m_X123  = 2.*Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      break;
    case 2010 : /* K- pi- pi+ */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      break;
    case 1020 : /* K- pi- K+ */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      break;
    case 3000 : /* pi- pi- pi+ mode */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      break;
    case 1002 : /* K0 pi- K0b */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      break;
    case  111 : /* K- pi0 K0 */
      m_X123  = 0.;
      m_ms123 = Mass2(1);
      m_G123  = 0;
      break;
    case  210 : /* pi0 pi0 K- mode */
      m_X123  = -2.*(Mass2(1)+Mass2(2));
      m_ms123 = Mass2(2);
      m_G123  = 1;
      break;
    case 1101 : /* pi- K0b pi0 */
      m_X123  = 0.;
      m_ms123 = Mass2(1);
      m_G123  = 0;
      break;
    case   12 : /* K- K0b K0 */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      break;
  }
}
 
Histogram * Tau_Three_Pseudo::KS::CreateGHistogram()
{
  // create file name
  char fn[512];
  if (m_G123) {
    sprintf(fn, "%s/PhaseSpaceFunctions/GQ2_MV13=%.3f_Mv13=%.3f_beta13=%.3f_MV23=%.3f_Mv23=%.3f_beta23=%.3f_run=%i.dat",
        m_path.c_str(),
        sqrt(m_msV[0]), sqrt(m_msv[0]), m_Beta[0],
        sqrt(m_msV[1]), sqrt(m_msv[1]), m_Beta[1], (m_running&2) );
  }
  else {
    sprintf(fn, "%s/PhaseSpaceFunctions/G(Q2)_noV13_MV23=%.3f_Mv23=%.3f_beta23=%.3f_run=%i.dat",
        m_path.c_str(),
        sqrt(m_msV[1]), sqrt(m_msv[1]), m_Beta[1], (m_running&2) );
  }

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {                          // if file exists
    // read table and create histogram
    msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
             <<"     Read G(q2) for   MV13="<<sqrt(m_msV[0])<<endl
             <<"                      Mv13="<<sqrt(m_msv[0])<<endl
             <<"                    beta13="<<m_Beta[0]<<endl
             <<"                      MV23="<<sqrt(m_msV[1])<<endl
             <<"                      Mv23="<<sqrt(m_msv[1])<<endl
             <<"                    beta23="<<m_Beta[1]<<endl
             <<"      from "<<fn<<"."<<endl;
    myHist = new Histogram( fn );
  }
  else {                            // if file does not exist
    // create histogram (i.e. table of values)
    msg.Out()<<"Create necessary phase space function for choosen parameters.\n"
             <<"This may take some time. Please wait..."<<endl;
    msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
             <<"     Create G(q2)  for MV13="<<sqrt(m_msV[0])<<endl
             <<"                       Mv13="<<sqrt(m_msv[0])<<endl
             <<"                     beta13="<<m_Beta[0]<<endl
             <<"                       MV23="<<sqrt(m_msV[1])<<endl
             <<"                       Mv23="<<sqrt(m_msv[1])<<endl
             <<"                     beta23="<<m_Beta[1]<<endl
             <<"      in "<<fn<<"."<<endl;
    double low (0.), up (3.2);
    int nbins  (50);
    double step = (up-low)/nbins;
    myHist = new Histogram( 0, low, up+step, nbins+1 );
    double q2  (low),
    phi (0.);
    while( q2<=up+step ) {
      phi = IntegralG(q2);              // get G value
      myHist->Insert( q2, phi );        // insert into histogram
      q2 += step;       
    }
    myHist->Output(string(fn));
  }
  return myHist;
}


double Tau_Three_Pseudo::KS::G( double s )
{
  double val;   
  p_G->Extrapolate( s+p_G->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  double p = val;
  return p;
}

double Tau_Three_Pseudo::KS::IntegralG( double Q2 )
{
  int Ns=500, Nt=500;                   // number of subintervals
  double sum (0.);
  double s_max = Q2 + Mass2(2-1) - 2.*sqrt(Q2*Mass2(2-1));
  double s_min = Mass2(1-1)+Mass2(3-1)+2.*sqrt(Mass2(1-1)*Mass2(3-1));
  double ds = (s_max-s_min)/Ns;
  double t_max (0.);
  double t_min (0.);
  double dt (0.);
  double s (s_min), t, u;
  double V12, V22, V1V2;
  
  while ( s<s_max ) {
    t_max = (   sqr( Q2-Mass2(1-1)-Mass2(2-1)+Mass2(3-1) ) 
              - sqr( Sqrt_Lambda(Q2, s, Mass2(2-1)) - Sqrt_Lambda(s,Mass2(1-1),Mass2(3-1)) ) 
            )/(4.*s);
    t_min = (   sqr( Q2-Mass2(1-1)-Mass2(2-1)+Mass2(3-1) ) 
              - sqr( Sqrt_Lambda(Q2, s, Mass2(2-1)) + Sqrt_Lambda(s,Mass2(1-1),Mass2(3-1)) ) 
            )/(4.*s);
    dt = (t_max-t_min)/Nt;
    t = t_min;
    while ( t<t_max ) {
      u    = Q2 - s - t + Mass2(1-1)+Mass2(2-1)+Mass2(3-1);
      V12  = 2.*(Mass2(1-1)+Mass2(2-1)) - s - sqr(u-t+Mass2(1-1)-Mass2(3-1))/(4.*Q2);
      V22  = 2.*(Mass2(1-1)+Mass2(2-1)) - t - sqr(u-s+Mass2(2-1)-Mass2(3-1))/(4.*Q2);
      V1V2 = (u-s-t+4.*Mass2(3-1))/2. - (u-t+Mass2(1-1)-Mass2(3-1))*(u-s+Mass2(2-1)-Mass2(3-1))/(4.*Q2);
      Complex BW_s = (BW_V(1-1,3-1,s)+m_Beta[0]*BW_v(1-1,3-1,s))/(1+m_Beta[0]);
      Complex BW_t = (BW_V(2-1,3-1,t)+m_Beta[1]*BW_v(2-1,3-1,t))/(1+m_Beta[1]);
      sum += ( V12 * norm(BW_s) + V22*norm(BW_t) + 2.*V1V2*real( BW_s*conj(BW_t)) )*ds*dt;
      t += dt;
    }
    s += ds;
  }
  return sum*Q2;
} 
 
Complex Tau_Three_Pseudo::KS::BW_A( double s )
{
  double MG_A (1.);
  double MG_AA (1.);
  if (m_running & 1) {      // axial width running?
    MG_A  = m_MA *m_GA  * G(s)/m_G_at_MA2 * pow( m_MA2/s, m_exp_alpha );
    MG_AA = m_MAA*m_GAA * G(s)/m_G_at_MAA2* pow( m_MAA2/s,m_exp_alpha );
  }
  else {                    // no axial width running
    MG_A  = m_MA *m_GA;
    MG_AA = m_MAA*m_GAA;
  }
  if (m_alpha!=0.) 
    return( 
        (           Tools::BreitWigner(s,m_MA2,MG_A) 
        + m_alpha * Tools::BreitWigner(s,m_MAA2,MG_AA) ) / (1.+m_alpha) 
    );
  return Tools::BreitWigner(s,m_MA2,MG_A);
}


Complex Tau_Three_Pseudo::KS::Tvector1( int a, int b, double x )
{
  Complex ret =           
                BW_V(a,b,x)*( 1.-(Mass2(a)-Mass2(b))/(3.*m_msV[a]) )
  + m_Beta[a]*( BW_v(a,b,x)*( 1.-(Mass2(a)-Mass2(b))/(3.*m_msv[a]) ) );
  return ret/(1.+m_Beta[a]);
}

Complex Tau_Three_Pseudo::KS::Tvector2( int a, int b, double x )
{
  Complex ret =
                BW_V(a,b,x)/m_msV[a] 
    + m_Beta[a]*BW_v(a,b,x)/m_msv[a];
  return ret*( 2.*(Mass2(a)-Mass2(b)) )/( 3.*(1.+m_Beta[a]) );    
}

Complex Tau_Three_Pseudo::KS::TSvector( int a, int b, int c, double Q2, double s, double t )
{
  Complex ret =
                BW_V(a,b,s) * (   Q2-2.*t-s+2.*Mass2(a)+Mass2(c) 
                                - (Mass2(a)-Mass2(b))/m_msV[a]*( Q2+t-Mass2(c)-Q2*(t-m_msV[a]) ) )
    + m_Beta[a]*BW_v(a,b,s) * (   Q2-2.*t-s+2.*Mass2(a)+Mass2(c) 
                                - (Mass2(a)-Mass2(b))/m_msv[a]*( Q2+t-Mass2(c)-Q2*(t-m_msv[a]) ) );
  return ret/(1.+m_Beta[a]);    
}

Complex Tau_Three_Pseudo::KS::Tgen( int a, int b, int c, double s, double t)
{
  if (m_G123) return( Tvector1(a-1,c-1,s) + Tvector2(b-1,c-1,t) );
  if (a==1) return Tvector2(b-1,c-1,t);
  if (a==2) return Tvector1(a-1,c-1,s);
  msg.Error()<<"ERROR in HADRONS::Tau_Three_Pseudo::KS::Tgen(a,b,c,s,t) : \n"
             <<"    Method was called with m_G123==0 and a != 1,2.\n"
             <<"    This must not happen, will abort."<<endl;
  abort();           
}

Complex Tau_Three_Pseudo::KS::Trho( double s )
{
  return (BW_V(1-1,3-1,s)+m_Beta[0]*BW_v(1-1,3-1,s))/(1+m_Beta[0]);
  return m_msV[0]/Complex(m_msV[0]-s,-sqrt(m_msV[0])*m_widthV[0]);
}

Complex Tau_Three_Pseudo::KS::TKstar( double s )
{
  return BW_V(2-1,3-1,s);
  return m_msV[1]/Complex(m_msV[1]-s,-sqrt(m_msV[1])*m_widthV[1]);
}

Complex Tau_Three_Pseudo::KS::FormFactor( int j, double Q2, double s, double t )
{
  Complex FF(0.,0.);
  switch( j ) {
    case 1 : { FF = BW_A(Q2)*Tgen(1,2,3,s,t);
//               Complex ref = 
//               BW_A(Q2)*(Trho(s)+2./3.*(Mass2(1)-Mass2(2))/m_msV[1]*TKstar(t));
//               BW_A(Q2)*(Trho(s)+2./3.*(Mass2(1)-Mass2(2))/m_msV[1]*BW_V(2-1,3-1,t));
//               PRINT_INFO(FF/ref);  
             break; }
    case 2 : { FF = BW_A(Q2)*Tgen(2,1,3,t,s); 
//             Complex ref =  
//               BW_A(Q2)*BW_V(2-1,3-1,t)*(1.-1./3.*(Mass2(1)-Mass2(2))/m_msV[1]);
//             PRINT_INFO(FF/ref);  
             break; }
    case 3 : { FF = m_X123 + BW_A(Q2)*(Q2-m_MA2)*m_ms123/(m_MA2*Q2) 
                         *( TSvector(1-1,3-1,2-1,Q2,s,t) + TSvector(2-1,3-1,1-1,Q2,t,s) );
             FF /= 2.*(Q2-m_ms123);          
//             Complex ref =
//               1./(2.*(Q2-Mass2(1)))*(
//                   Mass2(1)+Mass2(2) + BW_A(Q2)*(Q2-m_MA2)/(m_MA2*Q2)*(
//                     Trho(s)*Mass2(1)*(Q2-2.*t-s+2.*Mass2(2)+Mass2(1)) +
//                     BW_V(2-1,3-1,t)*(
//                       Mass2(1)*(Q2-2.*s-t+2.*Mass2(1)+Mass2(2)) -
//                       (Mass2(1)-Mass2(2))/m_msV[1] * (Mass2(1)*(Q2+t-Mass2(2))-Q2*(t-m_msV[1]))
//                       )
//                     )
//                   );
//             PRINT_INFO(FF/ref);
             break; }
  }
  return FF;
}
// tested: general form of a special channel 
//         = simplified form for this special channel
 
// General framework

Complex Tau_Three_Pseudo::FormFactor( int j, double Q2, double s, double t )
{
  return( p_ff->FormFactor(j,Q2,s,t) );     
}

double Tau_Three_Pseudo::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  Vec4D p1( _p[m_pseudo_1] ),
        p2( _p[m_pseudo_2] ),
        p3( _p[m_pseudo_3] );
  Vec4D Q( p1+p2+p3 );      
  double s = (p1+p3).Abs2(),
         t = (p2+p3).Abs2();
  double Q2 = Q.Abs2();
  double dot1 = Q*(p1-p3),
         dot2 = Q*(p2-p3);
  double d1 = dot1/Q2,
         d2 = dot2/Q2;
  Complex F1 = FormFactor( 1, Q2, s, t );        
  Complex F2 = FormFactor( 2, Q2, s, t );        
  Complex FS = FormFactor( 3, Q2, s, t );
  double ampl (0.);
  for( int h=0; h<4; h++ ) {    // helicity combination (nu,tau)
    ampl += norm(
             F.X( m_nutau, m_pseudo_1, 0, h, m_cR, m_cL ) * ( FS + F1*(1.-d1) - F2*d2 )
           + F.X( m_nutau, m_pseudo_2, 0, h, m_cR, m_cL ) * ( FS - F1*d1 + F2*(1.-d2) )
           + F.X( m_nutau, m_pseudo_3, 0, h, m_cR, m_cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) )
             ); 
  }
  F.Delete();
  return ampl/2.;
} // ME is invariant under exchange of two identical pseudos !


double Tau_Three_Pseudo::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*m_global*m_GF2/m_fpi2; 
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode (3prong)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Four_Pion_3::Tau_Four_Pion_3( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion1(-1),
  m_pion2(-1),
  m_pion3(-1),
  m_pion0(-1)
{
  m_metype = string("Tau_FourPion");
  // calc. mass^2
  m_ms[0] = sqr( p_flavs[0].PSMass() );
  for( int i=1; i<6; i++ ) {
    m_ms[i] = sqr( p_flavs[i].PSMass() );
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )  m_nutau = i;    // neutrino
    else {
      if( p_flavs[i].Kfcode() == kf::pi )               m_pion0 = i;    // pi0
      else {
        if( p_flavs[0].IsAnti()==p_flavs[i].IsAnti() )  m_pion3 = i;    // pi+
        else {
          if( m_pion1==-1 )                             m_pion1 = i;    // pi-
          else                                          m_pion2 = i;    // pi-
        }
      }
    }
  }
  // set internal numbering
  m_inter[0] = 0;
  m_inter[1] = m_pion3;
  m_inter[2] = m_pion0;
  m_inter[3] = m_pion1;
  m_inter[4] = m_pion2;
  m_inter[5] = m_nutau;
}
 
void Tau_Four_Pion_3::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vud2   = sqr(_md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real()));
  m_Vus2   = sqr(_md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real()) );
  m_GF2    = sqr(_md("GF", rpa.gen.ScalarConstant(string("GF"))) );
  m_cR     = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL     = Complex(0.,_md("a",1.)+_md("b",1.));

  Complex sum (0.,0.);
  char helps[20];
  double absol, phase;
  for (int i=0; i<4; i++) {
    sprintf( helps,"alpha_%i",i );
    absol = _md(helps+string("_abs"), 0. );
    phase = _md(helps+string("_phase"), 0. );
    m_Alpha[i]  = Complex( absol*cos(phase), absol*sin(phase) );
    sum += m_Alpha[i];
  }
  m_SumAlpha = sum;


  p_lorenz = new KS();
  p_lorenz->SetModelParameters( _md );
}
 
void Tau_Four_Pion_3::LorenzBase::SetPrivates( Complex * _x, ATOOLS::Vec4D * _p ) 
{
  p_X  = _x;
  p_p  = _p;
  p_p[0] = p_p[1]+p_p[2]+p_p[3]+p_p[4];         // = q
  m_q2 = p_p[0].Abs2();
  for (int i=2; i<=4; i++ ) {
    m_r[i] = p_p[0]-p_p[i];
    m_s[i] = (p_p[1]+p_p[i]).Abs2();
  }
  // redefinition of variables
  m_s[0] = (p_p[2]+p_p[3]).Abs2();          // = t3
  m_s[1] = (p_p[2]+p_p[4]).Abs2();          // = t4
  p_X[0] = p_X[1]+p_X[2]+p_X[3]+p_X[4];     // X(nu,q,0)
  // unused variables yet
  m_r[0] = m_r[1] = ATOOLS::Vec4D(0.,0.,0.,0.);
};
 
// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

void Tau_Four_Pion_3::KS::SetModelParameters( GeneralModel _md )
{
  m_fpi2    = sqr(_md("fpi",0.0924))*2.;        // redefine fpi
  m_grop    = _md("grop", 12.924);
  m_Go3p    = _md("Go3p", 1476.);

  m_mpi2    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_mpi02   = sqr( Flavour(kf::pi).PSMass() );

  m_MR      = _md("Mass_rho_770",  Flavour(kf::rho_770_plus).PSMass()  );
  m_MRR     = _md("Mass_rho_1450", Flavour(kf::rho_1450_plus).PSMass() );
  m_MRRR    = _md("Mass_rho_1700", Flavour(kf::rho_1700_plus).PSMass() );
  m_MO      = _md("Mass_omega",    Flavour(kf::omega_782).PSMass() );
  m_MF      = _md("Mass_f0_980",   Flavour(kf::f_0_980).PSMass() );
  m_MS      = _md("Mass_sigma",    Flavour(kf::f_0_980).PSMass()  );
  m_MA      = _md("Mass_a1_1260",  Flavour(kf::a_1_1260_plus).PSMass());
   
  m_GR      = _md("Width_rho_770",  Flavour(kf::rho_770_plus).Width()  );
  m_GRR     = _md("Width_rho_1450", Flavour(kf::rho_1450_plus).Width() );
  m_GRRR    = _md("Width_rho_1700", Flavour(kf::rho_1700_plus).Width() );
  m_GO      = _md("Width_omega",    Flavour(kf::omega_782).Width() );
  m_GF      = _md("Width_f0_980",   Flavour(kf::f_0_980).Width() );
  m_GS      = _md("Width_sigma",    Flavour(kf::f_0_980).Width()  );
  m_GA      = _md("Width_a1_1260",  Flavour(kf::a_1_1260_plus).Width());
   
  m_MR2     = m_MR*m_MR;
  m_MRR2    = m_MRR*m_MRR;
  m_MRRR2   = m_MRRR*m_MRRR;
  m_MO2     = m_MO*m_MO;
  m_MS2     = m_MS*m_MS;
  m_MF2     = m_MF*m_MF;
  m_MA2     = m_MA*m_MA;
   
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
  m_sigma   = _md("sigma", 0.);
   
  m_Frho    = _md("frho", 0.266)*m_MR2;

  m_R[0]=m_R[1] = 0.;
  m_R[2]        = -2.;
  m_R[3]=m_R[4] = 1.;

  char helps[20];
  double absol, phase;
  for (int i=0; i<4; i++) {
    sprintf( helps, "beta_omega_pi_%i", i );
    absol = _md( helps+string("_abs"), 0. );
    phase = _md( helps+string("_phase"), 0. );
    m_Beta_opi[i] = Complex( absol*cos(phase), absol*sin(phase) );
     
    sprintf( helps, "beta_a1_pi_%i", i );
    absol = _md( helps+string("_abs"), 0. );
    phase = _md( helps+string("_phase"), 0. );
    m_Beta_api[i] = Complex( absol*cos(phase), absol*sin(phase) );
     
    sprintf( helps, "beta_sigma_rho_%i", i );
    absol = _md( helps+string("_abs"), 0. );
    phase = _md( helps+string("_phase"), 0. );
    m_Beta_srh[i] = Complex( absol*cos(phase), absol*sin(phase) );
     
    sprintf( helps, "beta_f0_rho_%i", i );
    absol = _md( helps+string("_abs"), 0. );
    phase = _md( helps+string("_phase"), 0. );
    m_Beta_frh[i] = Complex( absol*cos(phase), absol*sin(phase) );
  }
}

Complex Tau_Four_Pion_3::KS::Fk( double x, Complex * _beta )
{
  Complex BW_R   = Tools::BreitWigner( x, m_MR2, sqrt(x)*m_GR );
  Complex BW_RR  = Tools::BreitWigner( x, m_MRR2, sqrt(x)*m_GRR );
  Complex BW_RRR = Tools::BreitWigner( x, m_MRRR2, sqrt(x)*m_GRRR );
  return( (_beta[0] + _beta[1]*BW_R + _beta[2]*BW_RR + _beta[3]*BW_RRR)/(_beta[0]+_beta[1]+_beta[2]+_beta[3]) );
}

Complex Tau_Four_Pion_3::KS::Trho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, m_MR2, m_GR, m_mpi2, 1. );
  Complex BW_RR  = Tools::BreitWigner( x, m_MRR2, m_GRR, m_mpi2, 1. );
  Complex BW_RRR = Tools::BreitWigner( x, m_MRRR2, m_GRRR, m_mpi2, 1. );
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Complex Tau_Four_Pion_3::KS::TTrho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, m_MR2, m_MR*m_GR )/m_MR2;
  Complex BW_RR  = Tools::BreitWigner( x, m_MRR2, m_MRR*m_GRR )/m_MRR2;
  return( BW_R + m_sigma*BW_RR );
}

double Tau_Four_Pion_3::KS::Dots( int k, int l )        // k=3,4; l=1,2,3,4 (l!=k)
  // special dot products used for anomalous part of J_omega
{
  int pre = l-1;                // predecessor of l
  if (pre==0) pre  = 4;
  if (pre==k) pre -= 1;     
  int suc  = l+1;               // successor of l
  if (suc==k)  suc += 1;
  if (suc==5)  suc  = 1;
  return( (m_r[k]*p_p[pre])*(p_p[suc]*p_p[k]) 
        - (m_r[k]*p_p[suc])*(p_p[pre]*p_p[k]) );
}

Complex Tau_Four_Pion_3::KS::OmegaPi()
{
  // chiral part
  Complex J_chi;
  Complex sum1 (0.,0.), sum2 (0.,0.);
  sum1 = Complex(0.,0.);
  for (int k=2; k<=4; k++) {
    sum2 = Complex(0.,0.);
    for (int l=2; l<=4; l++) if (l!=k) {
      sum2 += (p_X[0]-2.*p_X[l]) * (m_r[l]*(p_p[k]-p_p[1]))/m_r[l].Abs2();
    }
    sum1 += m_R[k]*Trho(m_s[k]) * (p_X[k]-p_X[1]-sum2);
  }
  J_chi = 2.*sqrt(3.)/m_fpi2*Trho(m_q2) * sum1;

  // anomalous part
  Complex J_a;
  sum1 = Complex(0.,0.);
  for( int k=3; k<=4; k++ ) {
    sum2 = Complex(0.,0.);
    for (int l=1; l<=4; l++) if (l!=k) {
      sum2 += p_X[l] * Dots(k,l);
    }
    sum1 += Tools::BreitWigner( m_r[k].Abs2(), m_MO2, m_MO*m_GO )/m_MO2*sum2;
  }
  J_a = m_Go3p*m_Frho*m_grop*TTrho(m_q2) * sum1;

  // total
  return( (J_chi + J_a)*Fk(m_q2,m_Beta_opi) );
}

Complex Tau_Four_Pion_3::KS::AonePi()
{
  Complex term1 (0.,0.), term2 (0.,0.);
  Complex A, B, C;
  Vec4D   P, Q, R;

  //  1st term
  P = p_p[1]-p_p[3];
  Q = p_p[1]-p_p[4];
  R = m_r[2];
  A = Tools::BreitWigner(m_s[3],m_MR2,m_GR,m_mpi2,1.);
  B = Tools::BreitWigner(m_s[4],m_MR2,m_GR,m_mpi2,1.);
  C = Tools::BreitWigner(R.Abs2(),m_MA2,m_GA,m_mpi2,1.);
  term1  = A*(p_X[1]-p_X[3]) 
         + B*(p_X[1]-p_X[4]) 
         - (p_X[0]-p_X[2])*( A*(R*P) + B*(R*Q) )/R.Abs2() 
         - p_X[0]*( A*(p_p[0]*P) 
                  + B*(p_p[0]*Q) 
                  - (p_p[0]*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() 
                  )/m_q2;
  term1 *= C;

  // 2nd and 3rd term
  int ind;
  Complex help;
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    P = p_p[1]-p_p[2];
    Q = p_p[ind]-p_p[2];
    R = m_r[k];
    A = Tools::BreitWigner(m_s[2],m_MR2,m_GR,m_mpi02,m_mpi2,1.);
    B = Tools::BreitWigner(m_s[ind-3],m_MR2,m_GR,m_mpi02,m_mpi2,1.);  //s[0,1]=t[3,4]
    C = Tools::BreitWigner(R.Abs2(),m_MA2,m_GA,m_mpi2,1.);
    help   = A*(p_X[1]-p_X[2]) 
           + B*(p_X[2]-p_X[ind]) 
           - (p_X[0]-p_X[k])*( A*(R*P) + B*(R*Q) )/R.Abs2() 
           - p_X[0]*( A*(p_p[0]*P) 
                    + B*(p_p[0]*Q) 
                    - (p_p[0]*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() 
                    )/m_q2;
    term2 += C*help;
  }

  // total 
  return (term1-term2)*Fk(m_q2,m_Beta_api);
}

Complex Tau_Four_Pion_3::KS::SigmaRho()
{
  int ind;
  Complex term (0.,0.);
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    BW_S = Tools::BreitWigner(m_s[k],m_MS2,m_GS,m_mpi2,1.);
    BW_R = Tools::BreitWigner(m_s[ind-3],m_MR2,m_GR,m_mpi02,m_mpi2,1.); //s[0,1]=t[3,4]
    term += BW_S*BW_R * ( p_X[2] - p_X[ind] + p_X[0]*(p_p[0]*(p_p[ind]-p_p[2])) );
  }
  return term*Fk(m_q2,m_Beta_srh);
}

Complex Tau_Four_Pion_3::KS::FzeroRho()
{
  int ind;
  Complex term (0.,0.);
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    BW_S = Tools::BreitWigner(m_s[k],m_MF2,m_GR,m_mpi2,1.);
    BW_R = Tools::BreitWigner(m_s[ind-3],m_MR2,m_GR,m_mpi02,m_mpi2,1.);
    term += BW_S*BW_R * ( p_X[2] - p_X[ind] + p_X[0]*(p_p[0]*(p_p[ind]-p_p[2])) );
  }
  return term*Fk(m_q2,m_Beta_frh);
}

Complex Tau_Four_Pion_3::KS::operator()( int number )
{
  switch( number ) {
    case 0 : return OmegaPi();
    case 1 : return AonePi();
    case 2 : return SigmaRho();
    case 3 : return FzeroRho();
  }
}

// General framework

double Tau_Four_Pion_3::Using_Hels( const Vec4D * _p )
{
  // internal numeration and convenient variables
  for (int i=1; i<=5; i++ ) m_p[i] = _p[m_inter[i]];

  // summation over helicities
  XYZFunc F(m_nout,_p,p_flavs);
  double ampl (0.);
  Complex help(0.,0.);
    for (int h=0; h<4; h++) {           // helicity combination (nutau, tau)
      for (int k=1; k<=4; k++) 
        m_X[k] = F.X( m_nutau, m_inter[k], 0, h, m_cR, m_cL );
      m_X[0] = m_X[1] + m_X[2] + m_X[3] + m_X[4];
      p_lorenz->SetPrivates( m_X, m_p );
      help = Complex(0.,0.);
      for (int k=0; k<m_ncontrib; k++) {
        help += m_Alpha[k] * (*p_lorenz)(k) / m_SumAlpha;
      }
      ampl += norm( help );
    }
  F.Delete();
  return ampl/2.;
}

double Tau_Four_Pion_3::operator()( const Vec4D * _p )
{
  double T = Using_Hels( _p );
  return m_GF2/2.*m_Vud2 * T;
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode (1prong)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Four_Pion_1::Tau_Four_Pion_1( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion1(-1),
  m_pion2(-1),
  m_pion3(-1),
  m_pion0(-1)
{
  m_metype = string("Tau_FourPion");
  // calc. mass^2
  m_ms[0] = sqr( p_flavs[0].PSMass() );
  for( int i=1; i<6; i++ ) {
    m_ms[i] = sqr( p_flavs[i].PSMass() );
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )  m_nutau = i;    // neutrino
    else {
      if( p_flavs[i].Kfcode() == kf::pi_plus )          m_pion3 = i;    // pi-
      else {
        if( m_pion0==-1 )                               m_pion0 = i;    // pi0
        else {
          if( m_pion1==-1 )                             m_pion1 = i;    // pi0
          else                                          m_pion2 = i;    // pi0
        }
      }
    }
  }
  // set internal numbering
  m_inter[0] = 0;
  m_inter[1] = m_pion3;
  m_inter[2] = m_pion0;
  m_inter[3] = m_pion1;
  m_inter[4] = m_pion2;
  m_inter[5] = m_nutau;
}
 
void Tau_Four_Pion_1::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vud2   = sqr(_md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real()) );
  m_Vus2   = sqr(_md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real()) );
  m_GF2    = sqr(_md("GF", rpa.gen.ScalarConstant(string("GF"))) );
  m_cR     = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL     = Complex(0.,_md("a",1.)+_md("b",1.));

  p_lorenz = new KS();
  p_lorenz->SetModelParameters( _md );
}
 
void Tau_Four_Pion_1::LorenzBase::SetPrivates( Complex * _X, ATOOLS::Vec4D * _p ) 
{
  p_X  = _X;
  p_p  = _p;
  p_p[0] = p_p[1]+p_p[2]+p_p[3]+p_p[4];           // = q
  m_q2 = p_p[0].Abs2();
  for (int i=2; i<=4; i++ ) {
    m_r[i] = p_p[0]-p_p[i];
    m_s[i] = (p_p[1]+p_p[i]).Abs2();
  }
  // redefinition of variables
  m_s[0] = (p_p[2]+p_p[3]).Abs2();            // = t3
  m_s[1] = (p_p[2]+p_p[4]).Abs2();            // = t4
  p_X[0] = p_X[1]+p_X[2]+p_X[3]+p_X[4];       // X(nu,q,0)
  // unused variables yet
  m_r[0] = m_r[1] = ATOOLS::Vec4D(0.,0.,0.,0.);
};
 
// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

void Tau_Four_Pion_1::KS::SetModelParameters( GeneralModel _md )
{
  m_fpi2    = sqr(_md("fpi",0.0924))*2.;        // redefine fpi with factor sqrt(2)

  m_mpi2    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_mpi02   = sqr( Flavour(kf::pi).PSMass() );

  m_MR      = _md("Mass_rho_770",  Flavour(kf::rho_770_plus).PSMass());
  m_MRR     = _md("Mass_rho_1450", Flavour(kf::rho_1450_plus).PSMass() );
  m_MRRR    = _md("Mass_rho_1700", Flavour(kf::rho_1700_plus).PSMass() );
   
  m_GR      = _md("Width_rho_770",  Flavour(kf::rho_770_plus).Width());
  m_GRR     = _md("Width_rho_1450", Flavour(kf::rho_1450_plus).Width() );
  m_GRRR    = _md("Width_rho_1700", Flavour(kf::rho_1700_plus).Width() );
   
  m_MR2     = m_MR*m_MR;
  m_MRR2    = m_MRR*m_MRR;
  m_MRRR2   = m_MRRR*m_MRRR;
   
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
}

Complex Tau_Four_Pion_1::KS::Trho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, m_MR2, m_GR, m_mpi2, 1. );
  Complex BW_RR  = Tools::BreitWigner( x, m_MRR2, m_GRR, m_mpi2, 1. );
  Complex BW_RRR = Tools::BreitWigner( x, m_MRRR2, m_GRRR, m_mpi2, 1. );
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Complex Tau_Four_Pion_1::KS::operator()()
{
  Complex J_chi;
  Complex sum1 (0.,0.), sum2 (0.,0.);
  sum1 = Complex(0.,0.);
  for (int k=2; k<=4; k++) {
    sum2 = Complex(0.,0.);
    for (int l=2; l<=4; l++) if (l!=k) {
      sum2 += (p_X[0]-2.*p_X[l]) * (m_r[l]*(p_p[k]-p_p[1]))/m_r[l].Abs2();
    }
    sum1 += Trho(m_s[k]) * (p_X[k]-p_X[1]-sum2);
  }
  J_chi = 2.*sqrt(3.)/m_fpi2*Trho(m_q2) * sum1;

  // total
  return( J_chi );
}

// General framework

double Tau_Four_Pion_1::Using_Hels( const Vec4D * _p )
{
  // internal numeration and convenient variables
  for (int i=1; i<=5; i++ ) {
    m_p[i] = _p[m_inter[i]];
  }

  // summation over helicities
  XYZFunc F(m_nout,_p,p_flavs);
  double ampl (0.);
  Complex help(0.,0.);
    for (int h=0; h<4; h++) {           // helicity combination (nutau, tau)
      for (int k=1; k<=4; k++) m_X[k] = F.X( m_nutau, m_inter[k], 0, h, m_cR, m_cL );
      m_X[0] = m_X[1] + m_X[2] + m_X[3] + m_X[4];
      p_lorenz->SetPrivates( m_X, m_p );
      help = Complex(0.,0.);
      help = (*p_lorenz)();
      ampl += norm( help );
    }
  F.Delete();
  return ampl/2.;
}

double Tau_Four_Pion_1::operator()( const Vec4D * _p )
{
  double T = Using_Hels( _p );
  return m_GF2/2.*m_Vud2 * T;
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  no formfactor  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// these are ME with same Lorenz structure as above
// BUT the formfactors are identical to one
// it is a first approach to cope with eta modes

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  eta 2pion mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 
Tau_Eta_Two_Pion::Tau_Eta_Two_Pion( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_eta(-1),
  m_pion(-1),
  m_pion0(-1)
{
  m_metype = string("Tau_EtaTwoPion");
  // calc. mass^2
  m_ms[0] = sqr( p_flavs[0].PSMass() );
  for( int i=1; i<5; i++ ) {
    m_ms[i] = sqr( p_flavs[i].PSMass() );
    if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 )  m_nutau = i;    // neutrino
    else {
      if( p_flavs[i].Kfcode() == kf::pi )               m_pion0 = i;    // pi0
      if( p_flavs[i].Kfcode() == kf::pi_plus )          m_pion  = i;    // pi
      if( p_flavs[i].Kfcode() == kf::eta )              m_eta   = i;    // eta
    }
  }
}
 
void Tau_Eta_Two_Pion::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vud2   = sqr(_md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real()) );
  m_GF2    = sqr(_md("GF", rpa.gen.ScalarConstant(string("GF"))) );
  m_cR     = Complex(0.,_md("a",1.)-_md("b",1.));
  m_cL     = Complex(0.,_md("a",1.)+_md("b",1.));
  m_fpi2   = sqr( _md("fpi", 0.0924));
  m_global =  sqr(2./3.);                           // global factor
}

double Tau_Eta_Two_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  Vec4D p1( _p[m_pion] ),
        p2( _p[m_pion0] ),
        p3( _p[m_eta] );
  Vec4D Q( p1+p2+p3 );      
  double s = (p1+p3).Abs2(),
         t = (p2+p3).Abs2();
  double Q2 = Q.Abs2();
  double dot1 = Q*(p1-p3),
         dot2 = Q*(p2-p3);
  double d1 = dot1/Q2,
         d2 = dot2/Q2;
  Complex F1 (1.,0.);
  Complex F2 (1.,0.);
  Complex FS (0.,0.);
  double ampl (0.);
  for( int h=0; h<4; h++ ) {    // helicity combination (nu,tau)
    ampl += norm(
             F.X( m_nutau, m_pion, 0, h, m_cR, m_cL )   * ( FS + F1*(1.-d1) - F2*d2 )
           + F.X( m_nutau, m_pion0, 0, h, m_cR, m_cL )  * ( FS - F1*d1 + F2*(1.-d2) )
           + F.X( m_nutau, m_eta, 0, h, m_cR, m_cL )    * ( FS - F1*(1.+d1) - F2*(1.+d2) )
             ); 
  }
  F.Delete();
  return ampl/2.;
} // ME is invariant under exchange of pi <-> pi0


double Tau_Eta_Two_Pion::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*m_GF2*m_Vud2/m_fpi2*m_global; 
}
 
