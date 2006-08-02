#include "Tau_Decay_MEs.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Histogram.H"
#include "Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

#define DEFINE_ME_GETTER(CLASS,NAME,TAG)                              \
DECLARE_GETTER(NAME,TAG,HD_ME_Base,Flavour_Info);                     \
HD_ME_Base* NAME::operator()(const Flavour_Info &parameters) const    \
{ return new CLASS(parameters.nout, parameters.flavs); }               \
void NAME::PrintInfo(std::ostream &str,const size_t width) const      \
{ str<<"implement me"; }                                              \

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  leptonic decay  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Lepton,Tau_Lepton_Getter,"Tau_Lepton");

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
  double GF = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
  m_global  = GF*SQRT_05;
  m_cR1 = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL1 = Complex(0.,_md("v",1.)+_md("a",1.));
  m_cR2 = Complex(0.,_md("v2",1.)-_md("a2",1.));
  m_cL2 = Complex(0.,_md("v2",1.)+_md("a2",1.));
}
 
void Tau_Lepton::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
  for( int h=0; h<16; h++ ) {      // for all helicity combinations
      _ampls_tensor->push_back( 
          F.Z( m_nutau, 0, m_lep, m_nulep, h, m_cR1, m_cL1, m_cR2, m_cL2 ) * m_global 
      );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2 3)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(m_nulep,1) );
  _indices->push_back( pair<int,int>(m_lep,1) );
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  1 pion/kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Pseudo,Tau_Pseudo_Getter,"Tau_Pseudo");

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
  double Vxx  = m_pionmode ? 
      _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real()) : 
      _md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real());
  double fxx  = m_pionmode ? _md("fpi", 0.0924) : _md("fK", 0.113);
  double GF   = _md("GF", rpa.gen.ScalarConstant(string("GF")) );
  m_global = fxx * GF * Vxx;
  m_cR   = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL   = Complex(0.,_md("v",1.)+_md("a",1.));
}
  
void Tau_Pseudo::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
  for( int h=0; h<4; h++ ) {
    _ampls_tensor->push_back( 
        -1.0 * Complex(0.0,1.0) * F.X(m_nutau,m_pion,0,h,m_cR,m_cL) * m_global
    );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pion does not have spin index
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  2 pion/kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Two_Pion,Tau_Two_Pion_Getter,"Tau_Two_Pion");

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
          p_flavs[i].Kfcode() == kf::K_L )              m_pion0 = i;
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
  m_ff       = int( _md("FORM_FACTOR", 1 ) );
  m_fpi      = _md("fpi", 0.0924 );
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real() );
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
  double CG  = m_pionmode ? 1. : SQRT_05;   // Clebsch-Gordon
  m_global   = GF * CG * Vud;           // GF * V_CKM * CG
  m_cR       = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL       = Complex(0.,_md("v",1.)+_md("a",1.));

  int    running  = int( _md("RUNNING_WIDTH", 1 ) );
  double MR       = _md("Mass_rho(770)+", Flavour(kf::rho_770_plus).PSMass() );
  double MRR      = _md("Mass_rho(1450)+", Flavour(kf::rho_1450_plus).PSMass() );
  double MRRR     = _md("Mass_rho(1700)+", Flavour(kf::rho_1700_plus).PSMass() );
  double GR       = _md("Width_rho(770)+", Flavour(kf::rho_770_plus).Width() );
  double GRR      = _md("Width_rho(1450)+", Flavour(kf::rho_1450_plus).Width() );
  double GRRR     = _md("Width_rho(1700)+", Flavour(kf::rho_1700_plus).Width() );
  m_R   = ResonanceFlavour( kf::rho_770_plus, MR, GR, running );
  m_RR  = ResonanceFlavour( kf::rho_1450_plus, MRR, GRR, running );
  m_RRR = ResonanceFlavour( kf::rho_1700_plus, MRRR, GRRR, running );

  m_m2_pi    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_m2_K     = sqr( Flavour(kf::K_plus).PSMass() );
   
  // coefficients for KS model
  m_beta     = _md("beta", 0. );
  m_gamma    = _md("gamma", 0. );

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
  if( m_ff == 1 ) {         // Breit-Wigner-rho
    Complex BWr   = m_R.BreitWigner(s);
    Complex BWrr  = m_RR.BreitWigner(s);
    Complex BWrrr = m_RRR.BreitWigner(s);
    ret = ( BWr + m_beta*BWrr + m_gamma*BWrrr )/( 1.+m_beta+m_gamma );
    return ret;
  }
  if( m_ff == 2 ) {         // Resonance Chiral Theory
    double MG_R, MG_RR, MG_RRR;
    Complex AA = A( m_m2_pi/s, m_m2_pi/m_R.Mass2() ) + 0.5*A( m_m2_K/s, m_m2_K/m_R.Mass2() );
    double expon = -1.*s/(96.*sqr(M_PI*m_fpi))*AA.real();
    Complex BW_1, BW_2, BW_3;
    if (m_R.Running()) {
      MG_R   = -m_gammaR  *1.*m_R.Mass2()  *s/(96.*sqr(M_PI*m_fpi)) * AA.imag();
      MG_RR  = -m_gammaRR *1.*m_RR.Mass2() *s/(96.*sqr(M_PI*m_fpi)) * AA.imag();
      MG_RRR = -m_gammaRRR*1.*m_RRR.Mass2()*s/(96.*sqr(M_PI*m_fpi)) * AA.imag();
      BW_1 = Tools::BreitWigner( s, m_R.Mass2(), MG_R );
      BW_2 = Tools::BreitWigner( s, m_RR.Mass2(), MG_RR );
      BW_3 = Tools::BreitWigner( s, m_RRR.Mass2(), MG_RRR );
    }
    else {
      MG_R   = m_R.MassWidth();
      MG_RR  = m_RR.MassWidth();
      MG_RRR = m_RRR.MassWidth();
      BW_1 = Tools::BreitWignerFix( s, m_R.Mass2(), MG_R );
      BW_2 = Tools::BreitWignerFix( s, m_RR.Mass2(), MG_RR );
      BW_3 = Tools::BreitWignerFix( s, m_RRR.Mass2(), MG_RRR );
    }
    ret = (BW_1+m_beta*BW_2+m_gamma*BW_3)/(1.+m_beta+m_gamma) * exp(expon);
    return ret;
  }
  return Complex(1.,0.);
}

void Tau_Two_Pion::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
  double  q2 = (_p[m_pion_ch] + _p[m_pion0] ).Abs2();
  Complex FF = FormFactor(q2);
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
      _ampls_tensor->push_back( 
          ( F.X(m_nutau,m_pion_ch,0,h,m_cR,m_cL)
          - F.X(m_nutau,m_pion0,0,h,m_cR,m_cL) ) * m_global*FF
      );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pion does not have spin index
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   pion-kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Pion_Kaon,Tau_Pion_Kaon_Getter,"Tau_Pion_Kaon");

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
  double Vus  =_md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real() );
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
  m_global   = GF*Vus/2.;
  m_cR       = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL       = Complex(0.,_md("v",1.)+_md("a",1.));
  m_Delta_KP = m_ms[m_kaon] - m_ms[m_pion];
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : p_ff = new RChT(_md);
             break;               // use RChT on own risk: not tested sufficiently
    case 1 : p_ff = new KS(_md);
             break;
  }
  p_ff->SetMasses2( m_ms[m_pion], m_ms[m_kaon], sqr(Flavour(kf::eta).PSMass()) );
}

Tau_Pion_Kaon::FF_Base::FF_Base( GeneralModel _md )
{
  int    running = int( _md("RUNNING_WIDTH", 1 ) );
  double MR      = _md("Mass_K*(892)+", Flavour(kf::K_star_892_plus).PSMass());
  double GR      = _md("Width_K*(892)+", Flavour(kf::K_star_892_plus).Width());
  double MR0     = _md("Mass_K(0)*(1430)+", Flavour(kf::K_0_star_1430_plus).PSMass());
  double GR0     = _md("Width_K(0)*(1430)+", Flavour(kf::K_0_star_1430_plus).Width());
  m_R       = ResonanceFlavour(kf::K_star_892_plus, MR, GR, running);
  m_R0      = ResonanceFlavour(kf::K_0_star_1430_plus, MR0, GR0, running );
  m_fpi2     = sqr( _md("fpi", 0.0924) );
}
 
void Tau_Pion_Kaon::FF_Base::SetMasses2( double _mPi2, double _mK2, double _mEta2 )
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
 

// Resonance Chiral Theory

Tau_Pion_Kaon::RChT::RChT( GeneralModel _md )
  : FF_Base(_md)
{
  m_MK2      = m_R.Mass2();
  m_GK       = m_R.Width();
  m_MK02     = m_R0.Mass2();
  m_GK0      = m_R0.Width();
  m_renorm2  = sqr( _md("renorm",_md("Mass_rho(770)+",Flavour(kf::rho_770_plus).PSMass())));
  m_cd       = _md("const_cd", 0.014);
  m_cm       = m_fpi2/4./m_cd;
}

double Tau_Pion_Kaon::RChT::MassWidthVector( double s )
{
  double ret (0.);
  if (s>sqr(m_mK+m_mPi))  ret += pow( Tools::Lambda(s,m_mK2,m_mPi2), 1.5 );
  if (s>sqr(m_mK+m_mEta)) ret += pow( Tools::Lambda(s,m_mK2,m_mEta2), 1.5 );
  ret *= m_MK2 /( 128.*M_PI*m_fpi2*sqr(s) );
  return ret;
}

double Tau_Pion_Kaon::RChT::MassWidthScalar( double s )
{
  double ret (0.);
  if (s>sqr(m_mK+m_mPi))  
    ret += 3./(32.*M_PI*sqr(m_fpi2)*m_MK02*s)*sqr( m_cd*(s-m_Sigma_KP )+m_cm*m_Sigma_KP )*
      pow( Tools::Lambda(s,m_mK2,m_mPi2), 1.5 );
  if (s>sqr(m_mK+m_mEta)) 
    ret += 1./(864.*M_PI*sqr(m_fpi2)*m_MK02*s)*sqr( m_cd*(s-7.*m_mK2-m_mPi2)+m_cm*(5.*m_mK2-3.*m_mPi2) ) *
      pow( Tools::Lambda(s,m_mK2,m_mEta2), 1.5 );
  return ret;
}

Complex Tau_Pion_Kaon::RChT::JBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double nu = sqrt( Tools::Lambda(s,MP2,MQ2) );
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
  Complex M_part = Mr(s,m_mK2,m_mPi2) + Mr(s,m_mK2,m_mEta2);
  Complex L_part = L(s,m_mK2,m_mPi2) + L(s,m_mK2,m_mEta2);
  double   expon = 3./2./m_fpi2 *( s*M_part.real() - L_part.real() );
  ret = BW * exp(expon);
  return ret;
}

Complex Tau_Pion_Kaon::RChT::ScalarFormFactor( double s )
{
  Complex ret(1.,0.);
  double  MG_K0 = MassWidthScalar(s);
  Complex    BW = Tools::BreitWigner( s, m_MK02, MG_K0 );
  Complex    F4 = 1./(8.*m_fpi2)
    * ( 5.*s - 2.*m_Sigma_KP - 3.*sqr(m_Delta_KP)/s )
    * JBar(s,m_mK2,m_mPi2,m_mK2+m_mPi2,m_mK2-m_mPi2)
    + 1./(24.*m_fpi2)
    * ( 3.*s - 2.*m_Sigma_KP - sqr(m_Delta_KP)/s )
    * JBar(s,m_mK2,m_mEta2,m_mK2+m_mEta2,m_mK2-m_mEta2);
  double  inter = 1. - (1.-m_fpi2/4./sqr(m_cd))*m_Sigma_KP/m_MK02;
  Complex  expon = Complex( F4.real(), F4.imag()/(1.+sqr(F4.imag())) );
  ret = BW * inter * exp(expon);
  return ret;
}

// Kuehn Santamaria Model

Tau_Pion_Kaon::KS::KS( GeneralModel _md )
  : FF_Base(_md)
{ }
 

Complex Tau_Pion_Kaon::KS::VectorFormFactor( double s )
{
  return m_R.BreitWigner(s);
}

Complex Tau_Pion_Kaon::KS::ScalarFormFactor( double s )
{
  return m_R0.BreitWigner(s);
}

// general framework

void Tau_Pion_Kaon::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )

{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
  double  q2 = (_p[m_pion]+_p[m_kaon]).Abs2();
  Complex FS = p_ff->ScalarFormFactor(q2);
  Complex FV = p_ff->VectorFormFactor(q2);
  Complex termK = m_Delta_KP/q2*(FS-FV)+FV;
  Complex termP = m_Delta_KP/q2*(FS-FV)-FV;
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
      _ampls_tensor->push_back( 
          ( F.X(m_nutau,m_kaon,0,h,m_cR,m_cL) * termK
          + F.X(m_nutau,m_pion,0,h,m_cR,m_cL) * termP ) * m_global
      );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pion/kaon do not have spin index
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  3 pseudo mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Three_Pseudo,Tau_Three_Pseudo_Getter,"Tau_Three_Pseudo");
 
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
}
 
void Tau_Three_Pseudo::SetModelParameters( GeneralModel _md ) 
{ 
  m_Vud      = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real() );
  m_Vus      = _md("Vus", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 1).real() );
  m_cR       = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL       = Complex(0.,_md("v",1.)+_md("a",1.));
  double fpi = _md("fpi", 0.0924 );
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")) ); 
   
  m_B123 = 0.;
  double A123;
  switch( m_mode ) {
    case 3000 : /* pi- pi- pi+ mode */
                A123 = m_Vud;
                break;
    case 1200 : /* pi0 pi0 pi- mode */
                A123 = m_Vud;
                break;
    case 1020 : /* K- pi- K+ */
                m_B123 = 2.;
                A123 = -0.5*m_Vud;
                break;
    case 1002 : /* K0 pi- K0b */
                m_B123 = -2.;
                A123 = -0.5*m_Vud;
                break;
    case  111 : /* K- pi0 K0 */
                A123 = 1.5*m_Vud*SQRT_05;
                break;
    case  210 : /* pi0 pi0 K- mode */
                A123 = m_Vus/4.;
                break;
    case 2010 : /* K- pi- pi+ */
                m_B123 = -2.;
                A123 = -0.5*m_Vus;
                break;
    case 1101 : /* pi- K0b pi0 */
                m_B123 = 4./3.;
                A123 = 1.5*m_Vus*SQRT_05;
                break;
    case   30 : /* K- K- K+ */
                A123 = m_Vus;
                break;
    default   : msg.Error()<<"Warning in HADRONS::Tau_Decay_MEs.C in Tau_Three_Pseudo::GetA123() :\n"
                           <<"     Obviously this three pseudoscalar channel (code "<<m_mode<<")\n"
                           <<"     doesn't have a global A123. Maybe it is not implemented yet.\n"
                           <<"     Take A123=1., will continue and hope for the best."<<endl;
                A123 = 1.;
                break;
  }
  m_global   = 2.*GF/(3.*fpi)*A123;
  double ms[3];
  ms[0] = m_ms[m_pseudo_1];
  ms[1] = m_ms[m_pseudo_2];
  ms[2] = m_ms[m_pseudo_3];
   
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : p_ff = new RChT(m_mode,m_path,_md,ms);
             break;
    case 1 : p_ff = new KS(m_mode,m_path,_md,ms);
             break;
  }
}
 
Tau_Three_Pseudo::FF_Base::FF_Base(int mode, std::string path, GeneralModel _md, double * _ms)
  : m_mode (mode)
{
  // set resonances and parameters
  m_ms[0] = _ms[0];
  m_ms[1] = _ms[1];
  m_ms[2] = _ms[2];
  m_deltas = 0;
  kf::code resA, resAA, resV[2], resVV[2];
  bool anomaly (false);
  kf::code resAnoV (kf::start), resAnoVV (kf::start), resAnoVVV (kf::start);
  switch( m_mode ) {
    case 3000 : /* pi- pi- pi+ mode */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      resV[0] = kf::rho_770;
      resV[1] = kf::rho_770;
      break;
    case 1200 : /* pi0 pi0 pi- mode */
      m_X123  = Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      resV[0] = kf::rho_770_plus;
      resV[1] = kf::rho_770_plus;
      break;
    case 1020 : /* K- pi- K+ */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      resV[0] = kf::rho_770;
      resV[1] = kf::K_star_892; 
      anomaly = true;
      break;
    case 1002 : /* K0 pi- K0b */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      resV[0] = kf::rho_770;
      resV[1] = kf::K_star_892_plus; 
      anomaly = true;
      break;
    case  111 : /* K- pi0 K0 */
      m_X123  = 0.;
      m_ms123 = Mass2(1);
      m_G123  = 0;
      resV[0] = kf::rho_770_plus;
      resV[1] = kf::rho_770_plus;
      break;
    case  210 : /* pi0 pi0 K- mode */
      m_X123  = -2.*(Mass2(1)+Mass2(2));
      m_ms123 = Mass2(2);
      m_G123  = 1;
      m_deltas= 1;
      resV[0] = kf::K_star_892_plus;
      resV[1] = kf::K_star_892_plus;
      break;
    case 2010 : /* K- pi- pi+ */
      m_X123  = Mass2(1)+Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      m_deltas= 1;
      resV[0] = kf::K_star_892; 
      resV[1] = kf::rho_770;
      anomaly = true;
      break;
    case 1101 : /* pi- K0b pi0 */
      m_X123  = 0.;
      m_ms123 = Mass2(1);
      m_G123  = 0;
      m_deltas= 1;
      resV[0] = kf::rho_770_plus;
      resV[1] = kf::rho_770_plus;
      anomaly = true;
      break;
    case   30 : /* K- K- K+ */
      m_X123  = 2.*Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      m_deltas= 1;
      resV[0] = kf::rho_770;
      resV[1] = kf::rho_770;
      break;
    case   12 : /* K- K0b K0 */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      m_deltas= 1;
      resV[0] = kf::rho_770_plus;
      resV[1] = kf::rho_770_plus;
      break;
  }
  resA = (m_deltas)? kf::K_1_1400_plus : kf::a_1_1260_plus;

  // higher resonances
  resAA = (m_deltas)? kf::K_1_1270_plus : kf::a_1_1260_plus;
  for( int i=0; i<2; ++i ) {
    if( resV[i] == kf::rho_770 )          resVV[i] = kf::rho_1450;
    if( resV[i] == kf::rho_770_plus )     resVV[i] = kf::rho_1450_plus;
    if( resV[i] == kf::K_star_892 )       resVV[i] = kf::K_star_1410;
    if( resV[i] == kf::K_star_892_plus )  resVV[i] = kf::K_star_1410_plus;
  }

  // set correct parameters 
  int running    = int( _md("RUNNING_WIDTH", 3 ) );     // running width
  double mA      = _md("Mass_"+Flavour(resA).IDName(),     Flavour(resA).PSMass());          // mass of axial resonance
  double mAA     = _md("Mass_"+Flavour(resAA).IDName(),    Flavour(resAA).PSMass());         // mass of axial resonance'
  double wA      = _md("Width_"+Flavour(resA).IDName(),     Flavour(resA).Width());          // width of axial resonance
  double wAA     = _md("Width_"+Flavour(resAA).IDName(),    Flavour(resAA).Width());         // width of axial resonance'
  if( m_deltas ) running = running&2;
  m_A = ResonanceFlavour( resA, mA, wA, running&1, path );
  m_AA = ResonanceFlavour( resAA, mAA, wAA, running&1, path );
  ResonanceFlavour help_rho = ResonanceFlavour(
      kf::rho_770_plus, 
      _md("Mass_rho(770)+", _md("Mass_rho(770)", Flavour(kf::rho_770_plus))),
      _md("Width_rho(770)+", _md("Width_rho(770)", Flavour(kf::rho_770_plus))),
      running&2 
      );
  ResonanceFlavour help_rhoP = ResonanceFlavour(
      kf::rho_1450_plus, 
      _md("Mass_rho(1450)+", _md("Mass_rho(1450)", Flavour(kf::rho_770_plus))),
      _md("Width_rho(1450)+", _md("Width_rho(1450)", Flavour(kf::rho_770_plus))),
      running&2 
      );
  double help_beta = _md("beta_rho(1450)+", _md("beta_rho(1450)", 0.) );
  m_A.InitialiseThreeBodyResonance(help_rho, help_rhoP, help_beta);
  m_AA.InitialiseThreeBodyResonance(help_rho, help_rhoP, help_beta);
  m_alpha        = _md("alpha_"+Flavour(resAA).IDName(),  0. );                              // weight factor for A'

  double mV[2], mVV[2], wV[2], wVV[2];
  for( int i=0; i<2; ++i ) {
    mV[i]       = _md("Mass_"+Flavour(resV[i]).IDName(),  Flavour(resV[i]).PSMass());       // mass^2 of vector resonance ij
    wV[i]       = _md("Width_"+Flavour(resV[i]).IDName(), Flavour(resV[i]).Width());        // width^2 of vector resonance ij
    mVV[i]      = _md("Mass_"+Flavour(resVV[i]).IDName(), Flavour(resVV[i]).PSMass());      // mass^2 of vector resonance' ij
    wVV[i]      = _md("Width_"+Flavour(resVV[i]).IDName(),Flavour(resVV[i]).Width());       // width^2 of vector resonance' ij
    m_V[i]      = ResonanceFlavour( resV[i], mV[i], wV[i], running&2 );
    m_VV[i]     = ResonanceFlavour( resVV[i], mVV[i], wVV[i], running&2 );
    m_Beta[i]   = _md("beta_"+Flavour(resVV[i]).IDName(),  0. );                            // weight factor for Vij'
  }

  resAnoV   = (m_deltas)? kf::K_star_892_plus : kf::rho_770_plus;
  resAnoVV  = (m_deltas)? kf::K_star_1410_plus : kf::rho_1450_plus;
  resAnoVVV = (m_deltas)? kf::K_star_1680_plus : kf::rho_1700_plus;

  if( anomaly ) {
    double MV   = _md("Mass_anomaly_"+Flavour(resAnoV).IDName(),   Flavour(resAnoV).PSMass()   );   // mass V
    double MVV  = _md("Mass_anomaly_"+Flavour(resAnoVV).IDName(),  Flavour(resAnoVV).PSMass()   );  // mass V'
    double MVVV = _md("Mass_anomaly_"+Flavour(resAnoVVV).IDName(), Flavour(resAnoVVV).PSMass()   ); // mass V'
    double GV   = _md("Width_anomaly_"+Flavour(resAnoV).IDName(),   Flavour(resAnoV).Width()   );   // width V
    double GVV  = _md("Width_anomaly_"+Flavour(resAnoVV).IDName(),  Flavour(resAnoVV).Width()   );  // width V'
    double GVVV = _md("Width_anomaly_"+Flavour(resAnoVVV).IDName(), Flavour(resAnoVVV).Width()   ); // width V'
    m_AnoV      = ResonanceFlavour( resAnoV, MV, GV, running&2, path );
    m_AnoVV     = ResonanceFlavour( resAnoVV, MVV, GVV, running&2, path );
    m_AnoVVV    = ResonanceFlavour( resAnoVVV, MVVV, GVVV, running&2, path );
    m_AlphaV    = _md("alpha_anomaly_K*(892)", _md("alpha_anomaly_K*(892)+", 0. ));   // alpha for K*
    m_BetaV[0]  = _md("beta_anomaly_"+Flavour(resAnoVV).IDName(), 0. );        // beta for V'
    m_BetaV[1]  = _md("gamma_anomaly_"+Flavour(resAnoVVV).IDName(), 0. );      // gamma for V''
  }
  else {
    m_AnoV      = ResonanceFlavour( kf::rho_770_plus, 0., 0., 0 );
    m_AnoVV     = ResonanceFlavour( kf::rho_770_plus, 0., 0., 0 );
    m_AnoVVV    = ResonanceFlavour( kf::rho_770_plus, 0., 0., 0 );
    m_AlphaV    = 0.;
    m_BetaV[0]  = 0.;
    m_BetaV[1]  = 0.;
  }

  m_fpi2      = sqr( _md("fpi", 0.0924) );          // pion decay constant
}
 

// Parameterisation
// DUMM, PICH, PORTOLES hep-ph/0312183 
 
Tau_Three_Pseudo::RChT::RChT(int mode, string path, GeneralModel _md, double * _ms)
  : FF_Base(mode,path,_md,_ms)
{
  // set resonances and parameters
  if( m_mode != 1200 && 
      m_mode != 3000 && 
      m_mode != 1020 ) {
    msg.Error()<<"Error: The mode "<<m_mode<<endl
      <<"     hasn't been implemented yet (RChT). Please use KS model."
      <<"     Don't know what to do. Will abort"<<endl;
    abort();
  }
   
  m_MO     = _md("Mass_omega(782)", Flavour(kf::omega_782).PSMass());
  m_MO2    = sqr(m_MO);
  m_GO     = _md("Width_omega(782)", Flavour(kf::omega_782).Width());
      
  m_l0     = _md("lambda0", 1.);                    // fit parameter lambda0
  m_gammaR = _md("gamma_rho(770)", 1.);              // global factor for rho width
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
}

double Tau_Three_Pseudo::RChT::MassWidthVector( int a, double s )
{ 
  if(m_V[0].Running()) {
    switch(m_mode) {
      case 3000:
      case 1200: return MassWidthVector( s );
      case 1020: return (a==0)? MassWidthVector(s) : m_V[1].MassWidth();
    }
    // default:
    msg.Error()<<"Warning: this form factor (RChT) for the three-pseudoe mode "<<m_mode<<"\n"
      <<"     hasn't been implemented yet. Please use KS model."<<endl;
  }
  return( m_V[a].MassWidth() );
}

double Tau_Three_Pseudo::RChT::MassWidthVector( double s )
{ 
  double MVGV (0.);
    if( s>4.*m_m2 )  MVGV += pow( 1.-4.*m_m2/s, 1.5 );
    if( s>4.*m_mK2 ) MVGV += pow( 1.-4.*m_mK2/s, 1.5 ) / 2.;
    MVGV *= m_gammaR*m_V[0].Mass2()*s/(96.*M_PI*m_fpi2); 
  return MVGV;
}

double Tau_Three_Pseudo::RChT::MassWidthAxial( double Q2 )
{
  if( !m_deltas && m_A.Running() )
    return(  m_A.OffShellMassWidth(Q2)*pow(m_A.Mass2()/Q2,m_exp_alpha-2.) );
  return m_A.MassWidth();
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
                   double MV2 = m_V[0].Mass2();
                   Complex alpha = 1. - 3./2.*x/Complex(x-MV2, MVGV_x);
                   Complex beta =  -3./2.*x/Complex(x-MV2, MVGV_x)
                     + F_Q2_x*(2.*Q2+x-u)/Complex(x-MV2, MVGV_x)
                     + F_Q2_y*(u-x)/Complex(y-MV2, MVGV_y);
                   return alpha - Q2/Complex(Q2-m_A.Mass2(),MAGA)*beta;
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
                             Complex FAchi = 1;
                             Complex FA1r  = 0.5 * m_FV*m_GV/m_fpi2 * ( 
                                  1./Complex(m_V[a].Mass2()-x,-1.*MassWidthVector(a,x))*
                                    ( 3.*x + Mass2(2)-Mass2(a) + (1.-2.*m_GV/m_FV)*(2.*Q2-2.*x-u+Mass2(a)-Mass2(b)) )
                                 +1./Complex(m_V[b].Mass2()-y,-1.*MassWidthVector(b,y))*
                                    ( 2.*(Mass2(b)-Mass2(2)) + (1.-2.*m_GV/m_FV)*(u-x+Mass2(2)-Mass2(b)) )    
                                 );
                             Complex FA2r  = -sqrt(2.) * m_FA*m_GV/m_fpi2 * Q2/Complex(m_A.Mass2()-Q2,-1.*MassWidthAxial(Q2)) * (
                                  1./Complex(m_V[a].Mass2()-x,-1.*MassWidthVector(a,x))*
                                    ( m_lsum*(-3.*x+Mass2(a)-Mass2(2)) + FFunc(x,Mass2(b),Q2)*(2.*Q2+x-u+Mass2(2)-Mass2(b)) )
                                 +1./Complex(m_V[b].Mass2()-y,-1.*MassWidthVector(b,y))*
                                    ( 2.*m_lsum*(Mass2(2)+Mass2(b)) + FFunc(y,Mass2(a),Q2)*(u-x+Mass2(b)-Mass2(2)) )
                                 );
                             return FAchi + FA1r + FA2r;
                           }
                   case 3: {                // pseudoscalar contribution
                             Complex FSchi = 3./2.* Tools::BreitWigner(Q2,Mass2(1),0.) *
                                                ( 1. + (Mass2(2)-u)/Q2 );
                             Complex FS1r  = 3./4.*sqr(m_GV)/m_fpi2 * Tools::BreitWigner(Q2,Mass2(1),0.)/Q2 *
                                                (   s*(t-u)/Complex(m_V[0].Mass2()-s,-1.*MassWidthVector(0,s))
                                                  + (t*(s-u)+(Q2-Mass2(0))*(Mass2(1)-Mass2(2)))/Complex(m_V[1].Mass2()-t,-1.*MassWidthVector(1,t)) 
                                                );
                             return FSchi + FS1r;
                           }
                   case 4: {                // vector contribution
                             Complex FVchi = 3./(4.*sqr(M_PI)*m_fpi2);
                             Complex FV1r  = 9./(32.*sqr(M_PI)*m_fpi2) * (
                                 Tools::BreitWigner(s,m_MO2,m_MO*m_GO)
                                 +Tools::BreitWigner(t,m_V[1].Mass2(),MassWidthVector(1,t)) );
                             Complex FV2r  = 3./(4.*Complex(m_V[0].Mass2()-Q2,-1.*MassWidthVector(0,Q2))) * (
                                 1./Complex(m_MO2-s,-1.*m_MO*m_GO) * (
                                   (1.-3.*m_V[0].Mass2()/(8.*sqr(M_PI)*m_fpi2))*(Q2+s)+Mass2(1) )
                                 + 1./Complex(m_V[1].Mass2()-t,-1.*MassWidthVector(1,t)) * (
                                   (1.-3.*m_V[0].Mass2()/(8.*sqr(M_PI)*m_fpi2))*(Q2+t)+Mass2(0) )
                                 );
                             return FVchi + FV1r + FV2r;
                           }
                 }
               }
  }
  return  (j>3)? Complex(0.,0.) : Complex(1.,0.); 
}

double Tau_Three_Pseudo::RChT::FFunc( double a, double b, double c)
{
  return m_l2 + m_l1*a/c - m_l0*b/c;
}

// Parameterisation
// DECKER, FINKEMEIER, MIRKES hep-ph/9310270

Tau_Three_Pseudo::KS::KS(int mode, string path, GeneralModel _md, double * _ms)
  : FF_Base(mode,path,_md,_ms)
{
} 
 
// methods for axial and scalar FF
 
Complex Tau_Three_Pseudo::KS::BW_V( int a, double s )
{
  if( !m_G123 && a==1 ) return Complex(0.,0.);          // BW_V23 = 0 if G123=0
  return m_V[a].BreitWigner(s);
}
 
Complex Tau_Three_Pseudo::KS::BW_VV( int a, double s )
{
  if( !m_G123 && a==1 ) return Complex(0.,0.);          // BW_V23' = 0 if G123=0
  return m_VV[a].BreitWigner(s);
}

Complex Tau_Three_Pseudo::KS::BW_A( double s )
{
//  PRINT_INFO(m_A.Name()<<" "<<m_A.Mass()<<" "<<m_AA.Name()<<" "<<m_AA.Mass() );
//  PRINT_INFO(m_alpha); abort();
  if (!IsZero(m_alpha)) {
    return( ( m_A.BreitWigner(s)+m_alpha*m_AA.BreitWigner(s) ) / (1.+m_alpha) );
  }
  return m_A.BreitWigner(s);
}


Complex Tau_Three_Pseudo::KS::Tvector1( int a, int b, double x )
{
  Complex ret =           
                BW_V(a,x) *( 1.-(Mass2(a)-Mass2(b))/(3.*m_V[a].Mass2()) )
  + m_Beta[a]*( BW_VV(a,x)*( 1.-(Mass2(a)-Mass2(b))/(3.*m_VV[a].Mass2()) ) );
  return ret/(1.+m_Beta[a]);
}

Complex Tau_Three_Pseudo::KS::Tvector2( int a, int b, double x )
{
  Complex ret = BW_V(a,x)/m_V[a].Mass2() + m_Beta[a]*BW_VV(a,x)/m_V[a].Mass2();
  return ret*( 2.*(Mass2(a)-Mass2(b)) )/( 3.*(1.+m_Beta[a]) );    
}

Complex Tau_Three_Pseudo::KS::TSvector( int a, int b, int c, double Q2, double s, double t )
{
  Complex ret =
                BW_V(a,s)  * (   m_ms123*( Q2-2.*t-s+2.*Mass2(a)+Mass2(c) )
                               - (Mass2(a)-Mass2(b))/m_V[a].Mass2()*( m_ms123*(Q2+s-Mass2(c))-Q2*(s-m_V[a].Mass2()) ) )
    + m_Beta[a]*BW_VV(a,s) * (   m_ms123*( Q2-2.*t-s+2.*Mass2(a)+Mass2(c) ) 
                               - (Mass2(a)-Mass2(b))/m_VV[a].Mass2()*( m_ms123*(Q2+s-Mass2(c))-Q2*(s-m_VV[a].Mass2()) ) );
  return ret/(1.+m_Beta[a]);
}

Complex Tau_Three_Pseudo::KS::Tgen( int a, int b, int c, double s, double t)
{
  return( Tvector1(a-1,c-1,s) + Tvector2(b-1,c-1,t) );
}

// methods for vector FF
 
Complex Tau_Three_Pseudo::KS::T_V( double q2 )
{
  Complex ret(0.,0.);
  ret  =           m_AnoV.BreitWigner(q2)
    + m_BetaV[0] * m_AnoVV.BreitWigner(q2)
    + m_BetaV[1] * m_AnoVVV.BreitWigner(q2);
  ret /= 1.+m_BetaV[0]+m_BetaV[1];
  return ret;
}

Complex Tau_Three_Pseudo::KS::T_V13V23( double s, double t )
{
  Complex ret(0.,0.);
  int n_rho (1-1), n_k (2-1);                       // what is the rho and K* resonance
  double x (s), y (t);
  if( m_mode==2010 ) { 
    n_rho = 2-1; n_k = 1-1;                         // swap them in Kpipi channel
    x=t; y=s;
  }    
  ret  = ( BW_V(n_rho,x) + m_Beta[n_rho]*BW_VV(n_rho,x) )/( 1.+m_Beta[n_rho] );
  ret += m_AlphaV * BW_V(n_k,y);
  ret /= 1.+m_AlphaV;
  return ret;
}

// handling of form factors

Complex Tau_Three_Pseudo::KS::FormFactor( int j, double Q2, double s, double t )
{
  Complex FF(0.,0.);
  switch( j ) {
    case 1 : { FF = BW_A(Q2)*Tgen(1,2,3,s,t);                   // axial
             break; }
    case 2 : { FF = BW_A(Q2)*Tgen(2,1,3,t,s);                   // axial
             break; }
    case 3 : { FF = m_X123 + BW_A(Q2)*(Q2-m_A.Mass2())/(m_A.Mass2()*Q2)     // scalar
                         *( TSvector(1-1,3-1,2-1,Q2,s,t) + TSvector(2-1,3-1,1-1,Q2,t,s) );
             FF /= 2.*(Q2-m_ms123);          
             break; }
    case 4 : { FF = 3./(8.*sqr(M_PI)*m_fpi2)*T_V(Q2)*T_V13V23(s,t); // vector
               break; }
  }
  return FF;
}

// General framework

Complex Tau_Three_Pseudo::FormFactor( int j, double Q2, double s, double t )
{
  return( p_ff->FormFactor(j,Q2,s,t) );     
}

void Tau_Three_Pseudo::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
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
  Complex FV = m_B123*FormFactor( 4, Q2, s, t );
  Vec4D qq = Tools::Cross(p1,p2,p3);
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
    _ampls_tensor->push_back( 
          ( F.X( m_nutau, m_pseudo_1, 0, h, m_cR, m_cL ) * ( FS + F1*(1.-d1) - F2*d2 )
          + F.X( m_nutau, m_pseudo_2, 0, h, m_cR, m_cL ) * ( FS - F1*d1 + F2*(1.-d2) )
          + F.X( m_nutau, m_pseudo_3, 0, h, m_cR, m_cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) ) 
          + Complex(0.,1.)*F.X( m_nutau, cross(p1,p2,p3), 0, h, m_cR, m_cL ) * FV
          ) * m_global
        );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pseudos do not have spin index
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode (3prong)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Four_Pion_3,Tau_Four_Pion_3_Getter,"Tau_Four_Pion_3");

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
  m_cR     = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL     = Complex(0.,_md("v",1.)+_md("a",1.));
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real());
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")));
  m_global   = GF*SQRT_05*Vud;

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

  p_lorenz = new KS(m_path,_md);
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

Tau_Four_Pion_3::KS::KS( string path, GeneralModel _md )
  : LorenzBase()
{
  m_fpi2    = sqr(_md("fpi",0.0924))*2.;        // redefine fpi
  m_grop    = _md("grop", 12.924);
  m_Go3p    = _md("Go3p", 1476.);

  m_mpi2    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_mpi02   = sqr( Flavour(kf::pi).PSMass() );

  double MR      = _md("Mass_rho(770)+",  Flavour(kf::rho_770_plus).PSMass()  );
  double MRR     = _md("Mass_rho(1450)+", Flavour(kf::rho_1450_plus).PSMass() );
  double MRRR    = _md("Mass_rho(1700)+", Flavour(kf::rho_1700_plus).PSMass() );
  double MO      = _md("Mass_omega(782)", Flavour(kf::omega_782).PSMass() );
  double MF      = _md("Mass_f(0)(980)",    Flavour(kf::f_0_980).PSMass() );
  double MS      = _md("Mass_sigma",      Flavour(kf::f_0_980).PSMass()  );
  double MA      = _md("Mass_a(1)(1260)+",  Flavour(kf::a_1_1260_plus).PSMass());
  double GR      = _md("Width_rho(770)+",  Flavour(kf::rho_770_plus).Width()  );
  double GRR     = _md("Width_rho(1450)+", Flavour(kf::rho_1450_plus).Width() );
  double GRRR    = _md("Width_rho(1700)+", Flavour(kf::rho_1700_plus).Width() );
  double GO      = _md("Width_omega(782)", Flavour(kf::omega_782).Width() );
  double GF      = _md("Width_f(0)(980)",    Flavour(kf::f_0_980).Width() );
  double GS      = _md("Width_sigma",      Flavour(kf::f_0_980).Width()  );
  double GA      = _md("Width_a(1)(1260)+",  Flavour(kf::a_1_1260_plus).Width());
  m_Rho = ResonanceFlavour( kf::rho_770_plus, MR, GR, 1 );
  m_RR  = ResonanceFlavour( kf::rho_1450_plus, MRR, GRR, 1 );
  m_RRR = ResonanceFlavour( kf::rho_1700_plus, MRRR, GRRR, 1 );
  m_O   = ResonanceFlavour( kf::omega_782, MO, GO, 0 );
  m_F   = ResonanceFlavour( kf::f_0_980, MF, GF, 1 );
  m_S   = ResonanceFlavour( kf::f_0_980, MS, GS, 1 );
  m_A   = ResonanceFlavour( kf::a_1_1260_plus, MA, GA, 1, path );
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
  m_sigma   = _md("sigma", 0.);
  m_A.InitialiseThreeBodyResonance(m_Rho, m_RR, m_beta);
   
   
  m_Frho    = _md("frho", 0.266)*m_Rho.Mass2();

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
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (_beta[0] + _beta[1]*BW_R + _beta[2]*BW_RR + _beta[3]*BW_RRR)/(_beta[0]+_beta[1]+_beta[2]+_beta[3]) );
}

Complex Tau_Four_Pion_3::KS::Trho( double x )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Complex Tau_Four_Pion_3::KS::TTrho( double x )
{
  Complex BW_R   = m_Rho.BreitWignerAlt(x);
  Complex BW_RR  = m_RR.BreitWignerAlt(x);
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
    sum1 += m_O.BreitWignerAlt(m_r[k].Abs2())*sum2;
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
  A = m_Rho.BreitWigner(m_s[3]);
  B = m_Rho.BreitWigner(m_s[4]);
  C = m_A.BreitWigner(R.Abs2());
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
    A = m_Rho.BreitWigner(m_s[2]);
    B = m_Rho.BreitWigner(m_s[ind-3]);
    C = m_A.BreitWigner(R.Abs2());
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
    BW_S = m_S.BreitWigner(m_s[k]);
    BW_R = m_Rho.BreitWigner(m_s[ind-3]);
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
    BW_S = m_F.BreitWigner(m_s[k]);
    BW_R = m_Rho.BreitWigner(m_s[ind-3]);
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

void Tau_Four_Pion_3::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // internal numeration and convenient variables
  for (int i=1; i<=5; i++ ) m_p[i] = _p[m_inter[i]];
  // create amplitudes tensor
  _ampls_tensor->clear();
  Complex help(0.,0.);
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
    // pre-calculate X-funcs
    for (int k=1; k<=4; k++) 
      m_X[k] = F.X( m_nutau, m_inter[k], 0, h, m_cR, m_cL );
    m_X[0] = m_X[1] + m_X[2] + m_X[3] + m_X[4];
    // sum over all contributions
    p_lorenz->SetPrivates( m_X, m_p );
    help = Complex(0.,0.);
    for (int k=0; k<m_ncontrib; k++) {
      help += m_Alpha[k] * (*p_lorenz)(k) / m_SumAlpha;
    }
    _ampls_tensor->push_back( help*m_global );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pions do not have spin index
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode (1prong)  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEFINE_ME_GETTER(Tau_Four_Pion_1,Tau_Four_Pion_1_Getter,"Tau_Four_Pion_1");

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
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real());
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")));
  m_global   = GF*SQRT_05*Vud;
  m_cR     = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL     = Complex(0.,_md("v",1.)+_md("a",1.));

  p_lorenz = new KS(_md);
}
 
void Tau_Four_Pion_1::LorenzBase::SetPrivates( Complex * _x, ATOOLS::Vec4D * _p ) 
{
  p_X  = _x;
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

Tau_Four_Pion_1::KS::KS( GeneralModel _md )
  : LorenzBase()
{
  m_fpi2    = sqr(_md("fpi",0.0924))*2.;        // redefine fpi with factor sqrt(2)
  m_mpi2    = sqr( Flavour(kf::pi_plus).PSMass() );
  m_mpi02   = sqr( Flavour(kf::pi).PSMass() );

  double MR      = _md("Mass_rho(770)+",  Flavour(kf::rho_770_plus).PSMass()  );
  double MRR     = _md("Mass_rho(1450)+", Flavour(kf::rho_1450_plus).PSMass() );
  double MRRR    = _md("Mass_rho(1700)+", Flavour(kf::rho_1700_plus).PSMass() );
  double GR      = _md("Width_rho(770)+",  Flavour(kf::rho_770_plus).Width()  );
  double GRR     = _md("Width_rho(1450)+", Flavour(kf::rho_1450_plus).Width() );
  double GRRR    = _md("Width_rho(1700)+", Flavour(kf::rho_1700_plus).Width() );
  m_Rho = ResonanceFlavour( kf::rho_770_plus, MR, GR, 1 );
  m_RR  = ResonanceFlavour( kf::rho_1450_plus, MRR, GRR, 1 );
  m_RRR = ResonanceFlavour( kf::rho_1700_plus, MRRR, GRRR, 1 );
   
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
}

Complex Tau_Four_Pion_1::KS::Trho( double x )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
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

void Tau_Four_Pion_1::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // internal numeration and convenient variables
  for (int i=1; i<=5; i++ ) {
    m_p[i] = _p[m_inter[i]];
  }
  // create amplitudes tensor
  _ampls_tensor->clear();
  Complex help(0.,0.);
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
    // pre-calculate X-funcs
    for (int k=1; k<=4; k++) 
      m_X[k] = F.X( m_nutau, m_inter[k], 0, h, m_cR, m_cL );
    m_X[0] = m_X[1] + m_X[2] + m_X[3] + m_X[4];
    // get Lorentz structure
    p_lorenz->SetPrivates( m_X, m_p );
    help = Complex(0.,0.);
    help = (*p_lorenz)();
    _ampls_tensor->push_back( help*m_global );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pions do not have spin index
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
  m_cR       = Complex(0.,_md("v",1.)-_md("a",1.));
  m_cL       = Complex(0.,_md("v",1.)+_md("a",1.));
  double fpi = _md("fpi", 0.0924);
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")));
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real());
  m_global =  2./3.*GF*Vud/fpi;  
}

void Tau_Eta_Two_Pion::operator()( 
    const Vec4D         * _p, 
    vector<Complex>     * _ampls_tensor, 
    vector<pair<int,int> > * _indices,
    int k0_n )
{
  XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // create amplitudes tensor
  _ampls_tensor->clear();
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
  Complex F1 = Complex(1.,0.);
  Complex F2 = Complex(1.,0.);
  Complex FS = Complex(0.,0.);
  for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
      _ampls_tensor->push_back( 
          ( F.X( m_nutau, m_pion, 0, h, m_cR, m_cL ) * ( FS + F1*(1.-d1) - F2*d2 )
          + F.X( m_nutau, m_pion0, 0, h, m_cR, m_cL ) * ( FS - F1*d1 + F2*(1.-d2) )
          + F.X( m_nutau, m_eta, 0, h, m_cR, m_cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) ) ) * m_global
      );
  }
  F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
  _indices->clear();
  _indices->push_back( pair<int,int>(0,1) );
  _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pseudos do not have spin index
}
