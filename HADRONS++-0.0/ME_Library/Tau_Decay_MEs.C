#include "Tau_Decay_MEs.H"
#include "Message.H"
#include "XYZFuncs.H"
#include "Histogram.H"
#include "Traces.H"

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
		p_flavs[i].Kfcode() == kf::mu ) { m_lep = i; break; }			// that's the lepton
  }
  // find the corresponding neutrinos
  for( int i=1; i<4; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
	if( p_flavs[i].Kfcode() == p_flavs[m_lep].Kfcode()+1 ) m_nulep = i;
  }
}
 
void Tau_Lepton::SetModelParameters( GeneralModel _md ) 
{ 
  GF   = _md.pm.GF;
  cR1 = (0.,_md.pm.a-_md.pm.b);
  cL1 = (0.,_md.pm.a+_md.pm.b);
  cR2 = (0.,_md.pm.a2-_md.pm.b2);
  cL2 = (0.,_md.pm.a2+_md.pm.b2);
  a   = _md.pm.a;
  b   = _md.pm.b;
  a2  = _md.pm.a2;
  b2  = _md.pm.b2;
}
 
double Tau_Lepton::Using_Hels( const Vec4D * _p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  double MW2 = sqr(Flavour(kf::W).Mass());
  for( int h1=0; h1<4; h1++ ) for( int h2=0; h2<4; h2++ ) {
	int h = (h1<<2)+h2;		// helicity combination (nutau,tau,lep,nulep)
	ret += norm( F.Z( m_nutau, 0, m_lep, m_nulep, h, cR1, cL1, cR2, cL2 ) );
  }
  F.Delete();
  return ret*0.25;
} // its value is 100% identical to simple traces
 
double Tau_Lepton::Using_Traces_Simple( const Vec4D *_p )
{
  double M = p_masses[0];
  double m = p_masses[m_lep]; 
  Vec4D P = _p[0];
  Vec4D k1 = _p[m_nutau];
  Vec4D p = _p[m_lep];
  Vec4D k2 = _p[m_nulep]; 
  Vec4D q = P-k1;
  double q2 = q.Abs();
  double dot1=0., dot2=0.;
  double M2=0., M2_part=0.;
  double MW2 = sqr(Flavour(kf::W).Mass());
  Complex term1(0.,0.);
  Complex term11(0.,0.);
  Complex term12(0.,0.);
  Complex term2(0.,0.);
  
  // | VA structure |^2
  Complex sum(0.,0.);
  for( int lambda=0; lambda<4; lambda++ ) for( int rho=0; rho<4; rho++ ) {
	term11  = Trace(k1,lambda,P,rho) - Trace5(k1,lambda,P,rho);
	term11 *= sqr(a+b)/4.;
	term12  = Trace(k1,lambda,P,rho) + Trace5(k1,lambda,P,rho); 
	term12 *= sqr(a-b)/4.;
	term1 = term11 + term12;
	term2 = a2*2.*Trace(p,lambda,k2,rho) - b2*2.*Trace5(p,lambda,k2,rho);
	sum += term1*term2*G(lambda,lambda)*G(rho,rho);
  }
  M2_part = sum.real()/128.;
  M2 += M2_part;
  
  return M2*64.;
}

double Tau_Lepton::Using_Traces( const Vec4D *_p )
{
  double M = p_masses[0];
  double m = p_masses[m_lep]; 
  Vec4D P = _p[0];
  Vec4D k1 = _p[m_nutau];
  Vec4D p = _p[m_lep];
  Vec4D k2 = _p[m_nulep]; 
  Vec4D S =  Vec4D(0.,0.,0.,1.);
  Vec4D q = P-k1;
  double q2 = q.Abs();
  double dot1=0., dot2=0.;
  double M2=0., M2_part=0.;
  double MW2 = sqr(Flavour(kf::W).Mass());
  Complex term1(0.,0.);
  Complex term11(0.,0.);
  Complex term12(0.,0.);
  Complex term2(0.,0.);
  
  // | VA structure |^2
  Complex sum(0.,0.);
  for( int lambda=0; lambda<4; lambda++ ) for( int rho=0; rho<4; rho++ ) {
	term11  = Trace(k1,lambda,P-M*S,rho) - Trace5(k1,lambda,P-M*S,rho);
	term11 *= sqr(a+b)/4.;
	term12  = Trace(k1,lambda,P+M*S,rho) + Trace5(k1,lambda,P+M*S,rho); 
	term12 *= sqr(a-b)/4.;
	term1 = term11 + term12;
	term2 = a2*2.*Trace(p,lambda,k2,rho) - b2*2.*Trace5(p,lambda,k2,rho);
	sum += term1*term2*G(lambda,lambda)*G(rho,rho);
  }
  M2_part = sum.real()/128.;
  M2 += M2_part;
  
  return M2*64.;
}

double Tau_Lepton::Using_Traces_2( const Vec4D *_p )
{
  double M = p_masses[0];
  double m = p_masses[m_lep]; 
  Vec4D P = _p[0];
  Vec4D k1 = _p[m_nutau];
  Vec4D p = _p[m_lep];
  Vec4D k2 = _p[m_nulep]; 
  Vec4D S =  Vec4D(0.,0.,0.,1.);
  Vec4D q = P-k1;
  double q2 = q.Abs();
  double dot1=0., dot2=0.;
  double M2=0., M2_part=0.;
  Complex term1  (0.,0.);
  Complex term11 (0.,0.);
  Complex term12 (0.,0.);
  Complex term2  (0.,0.);
  // | VA structure |^2
  Complex sum    (0.,0.);
  for( int lambda=0; lambda<4; lambda++ ) for( int rho=0; rho<4; rho++ ) {
	term11  = Trace(k1,lambda,P-M*S,rho) - Trace5(k1,lambda,P-M*S,rho);
	term11 *= sqr(a+b)/4.;
	term12  = Trace(k1,lambda,P+M*S,rho) + Trace5(k1,lambda,P+M*S,rho); 
	term12 *= sqr(a-b)/4.;
	term1 = term11 + term12;
	term2 = a2*2.*Trace(p,lambda,k2,rho) - b2*2.*Trace5(p,lambda,k2,rho);
	sum += term1*term2*G(lambda,lambda)*G(rho,rho);
  }
  M2_part = sum.real()/128.;
  M2 += M2_part;
  return M2*64.;
}

double Tau_Lepton::operator()( const Vec4D *_p )
{
  double T (1.);
  T = Using_Hels(_p);
  return T*sqr(GF);
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
	else m_pion = i;							// that's the pion
  }
}
 
void Tau_Pseudo::SetModelParameters( GeneralModel _md ) 
{ 
  Vxx  = ( p_flavs[m_pion].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.Vud : _md.pm.Vus;
  fxx  = ( p_flavs[m_pion].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.fpi : _md.pm.fK;
  GF   = _md.pm.GF;
  cR   = (0.,_md.pm.a-_md.pm.b);
  cL   = (0.,_md.pm.a+_md.pm.b);
}
  
double Tau_Pseudo::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {
	ret += norm( F.X(m_nutau,m_pion,0,h,cR,cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Pseudo::operator()( const Vec4D *_p )
{
  double T(1.);
  double q2 = _p[m_pion].Abs2(); 
  T = Using_Hels(_p); 
  return T*0.5*sqr(GF*Vxx*fxx); 
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
  m_metype = string("Tau_TwoPion");
  for( int i=1; i<4; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
	else {
	  if( p_flavs[i] == kf::pi || p_flavs[i] == kf::K ) m_pion0 = i;
	  else m_pion_ch = i;
	}
  }
}

void Tau_Two_Pion::SetModelParameters( GeneralModel _md ) 
{ 
  m        = p_masses[m_pion_ch];
  m2       = m*m;
  m_md.run = _md.run;
  m_md.ff  = _md.ff; 
  Vud  	 = _md.pm.Vud; 
  fxx  	 = ( p_flavs[m_pion_ch].Kfcode() == ATOOLS::kf::pi_plus )? _md.pm.fpi : _md.pm.fK;
  GF   	 = _md.pm.GF;
  frho 	 = _md.pm.frho;
  grpp 	 = _md.pm.grpp;
  CG   	 = ( p_flavs[m_pion_ch].Kfcode() == ATOOLS::kf::pi_plus )? 1. : 0.5; 	// Clebsch-Gordon
  cR     = (0.,_md.pm.a - _md.pm.b);
  cL     = (0.,_md.pm.a + _md.pm.b);

  MR 	= _md.pm.MR;
  GR 	= _md.pm.GR;
  MRR 	= _md.pm.MRR;
  GRR 	= _md.pm.GRR;
  MRRR 	= _md.pm.MRRR;
  GRRR 	= _md.pm.GRRR;
  MO 	= _md.pm.MO;
  GO 	= _md.pm.GO;

  MR2 	= MR*MR;
  MRR2 	= MRR*MRR;
  MRRR2 = MRRR*MRRR;
  MO2 	= MO*MO;
   
  beta  	= _md.pm.beta;
  delta 	= _md.pm.delta;
  lambda 	= _md.pm.lambda;
  gamma 	= _md.pm.gamma;

  gammaR	= _md.pm.gammaR;
  gammaRR 	= _md.pm.gammaRR;
  gammaRRR	= _md.pm.gammaRRR;
}
 
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
  if( m_md.ff == 1 ) {			// Breit-Wigner-rho
	double MG_R   = MR*GR;
	double MG_RR  = MRR*GRR;
	double MG_RRR = MRRR*GRRR;
	if (m_md.run) {
	  MG_R   = Tools::OffShellMassWidth( s, MR2, GR, m2, 1. );
	  MG_RR  = Tools::OffShellMassWidth( s, MRR2, GRR, m2, 1. );
	  MG_RRR = Tools::OffShellMassWidth( s, MRRR2, GRRR, m2, 1. );
	}
	Complex BWr   = Tools::BreitWigner( s, MR2, MG_R );
	Complex BWrr  = Tools::BreitWigner( s, MRR2, MG_RR );
	Complex BWrrr = Tools::BreitWigner( s, MRRR2, MG_RRR );
	Complex BWo   = Tools::BreitWigner( s, MO2, MO*GO );
	Complex help = BWr * (1.+delta*BWo)/(1.+delta); 
	ret = ( help + beta*BWrr + gamma*BWrrr )/( 1.+beta+gamma );
  }
  if( m_md.ff == 2 ) {			// Resonance Chiral Theory
	double m2_pi = sqr( Flavour(kf::pi_plus).Mass() );
	double m2_K  = sqr( Flavour(kf::K_plus).Mass() );
	Complex AA = A( m2_pi/s, m2_pi/MR2 ) + 0.5*A( m2_K/s, m2_K/MR2 );
	double expon = -1.*s/(96.*sqr(M_PI*fxx))*AA.real();
	double MG_R, MG_RR, MG_RRR;
	if (m_md.run) {
	  MG_R 	 = -gammaR  *1.*MR2  *s/(96.*sqr(M_PI*fxx)) * AA.imag();
	  MG_RR	 = -gammaRR *1.*MRR2 *s/(96.*sqr(M_PI*fxx)) * AA.imag();
	  MG_RRR = -gammaRRR*1.*MRRR2*s/(96.*sqr(M_PI*fxx)) * AA.imag();
	}
	else  {
	  MG_R   = MR*GR;
	  MG_RR  = MRR*GRR;
	  MG_RRR = MRRR*GRRR;
	}
	Complex BW_1 = Tools::BreitWigner( s, MR2, MG_R );
	Complex BW_2 = Tools::BreitWigner( s, MRR2, MG_RR );
	Complex BW_3 = Tools::BreitWigner( s, MRRR2, MG_RRR );
	ret = (BW_1+beta*BW_2+gamma*BW_3)/(1.+beta+gamma) * exp(expon);
  }
  return ret;
}

double Tau_Two_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  for( int h=0; h<4; h++ ) {				// only nutau and tau have a helicity !
	ret += norm( F.X(m_nutau,m_pion_ch,0,h,cR,cL)
				 - F.X(m_nutau,m_pion0,0,h,cR,cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Two_Pion::operator()( const Vec4D *_p )
{
  double T = Using_Hels(_p);
  double q2 = (_p[m_pion_ch] + _p[m_pion0] ).Abs2();
  Complex FF = FormFactor(q2);
  return T*0.5*sqr(GF*Vud)*norm(FF)*CG;
}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//   pion-kaon mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tau_Pion_Kaon::Tau_Pion_Kaon( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion(-1),
  m_kaon(-1)
{
  m_metype = string("Tau_PionKaon");
  m_ms[0] = sqr(p_flavs[0].Mass());
  for( int i=1; i<4; i++ ) {
	m_ms[i] = sqr(p_flavs[i].Mass());
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;           // neutrino
	else {
	  if( p_flavs[i] == kf::pi || p_flavs[i] == kf::pi_plus ) m_pion = i;     // pion
	  else m_kaon = i;                                                        // kaon
	}
  }
}

void Tau_Pion_Kaon::SetModelParameters( GeneralModel _md ) 
{ 
  Vus2     = sqr(_md.pm.Vud); 
  fpi2 	   = sqr(_md.pm.fpi);
  GF2  	   = sqr(_md.pm.GF);
  cR       = (0.,_md.pm.a - _md.pm.b);
  cL       = (0.,_md.pm.a + _md.pm.b);
  Delta_KP = m_ms[m_kaon] - m_ms[m_pion];
  switch( _md.ff ) {
	case 2 : p_ff = new RChT();
			 break;
	case 1 : p_ff = new KS();
			 break;
  }
  p_ff->SetModelParameters( _md );
  p_ff->SetMasses2( m_ms[m_pion], m_ms[m_kaon], sqr(Flavour(kf::eta).Mass()) );
}

// Resonance Chiral Theory

void Tau_Pion_Kaon::RChT::SetModelParameters( GeneralModel _md ) 
{
  MK2    = sqr(_md.pm.MK);
  GK     = _md.pm.GK;
  MK02   = sqr(_md.pm.MK0);
  GK0    = _md.pm.GK0;
  fpi2   = sqr(_md.pm.fpi);
  renorm2  = sqr(_md.pm.dummy1);
}

void Tau_Pion_Kaon::RChT::SetMasses2( double _mPi2, double _mK2, double _mEta2 )
{
  mPi2  = _mPi2;
  mK2   = _mK2;
  mEta2 = _mEta2;
  mPi   = sqrt(mPi2);
  mK    = sqrt(mK2);
  mEta  = sqrt(mEta2);
  Sigma_KP = mK2+mPi2;
  Delta_KP = mK2-mPi2;
}
 
double Tau_Pion_Kaon::RChT::MassWidthVector( double s )
{
  double ret (0.);
  if (s>sqr(mK+mPi))  ret += pow( Tools::Lambda(s,mK2,mPi2), 1.5 );
  if (s>sqr(mK+mEta)) ret += pow( Tools::Lambda(s,mK2,mEta2), 1.5 );
  ret *= MK2 /( 128.*M_PI*fpi2*sqr(s) );
}

double Tau_Pion_Kaon::RChT::MassWidthScalar( double s )
{
  double G_s = GK0*MK02/s*pow( Tools::Lambda(s,mK2,mPi2)/Tools::Lambda(MK02,mK2,mPi2), 0.5 );
  return sqrt(MK02) * G_s;
}

Complex Tau_Pion_Kaon::RChT::JBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double nu = sqrt( sqr(s) + sqr(MP2) + sqr(MQ2) - 2.*s*Sigma - 2.*s*MP2*MQ2 );
  double  J = 2. + Delta/s*log(MQ2/MP2) - Sigma/Delta*log(MQ2/MP2) - nu/s*log( (sqr(s+nu)-sqr(Delta))/(sqr(s-nu)-sqr(Delta)) ); 
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
  double mu2   = renorm2;
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
  Complex     BW = Tools::BreitWigner( s, MK2, MG_K );
  Complex M_part = Mr(s,mK2,mPi2) - Mr(s,mK2,mEta2);
  Complex L_part = L(s,mK2,mPi2) - L(s,mK2,mEta2);
  double   expon = 3./2./fpi2 *( s*M_part.real() - L_part.real() );
  ret = BW * exp(expon);
  return ret;
}

Complex Tau_Pion_Kaon::RChT::ScalarFormFactor( double s )
{
  Complex ret(1.,0.);
  double  MG_K0 = MassWidthScalar(s);
  Complex    BW = Tools::BreitWigner( s, MK02, MG_K0 );
  double    cd2 = sqr(0.032);
  Complex    F4 = 1./(8.*fpi2)*( 5.*s - 2.*Sigma_KP - 3.*sqr(Delta_KP)/s )*JBar(s,mK2,mPi2,mK2+mPi2,mK2-mPi2)
	            + 1./(24.*fpi2)*( 3.*s - 2.*Sigma_KP - sqr(Delta_KP)/s )*JBar(s,mK2,mEta2,mK2+mEta2,mK2-mEta2);
  double  inter = 1. - (1.-fpi2/4./cd2)*Sigma_KP/MK02;
  Complex  expon = Complex( F4.real(), F4.imag()/(1.+sqr(F4.imag())) );
  ret = BW * inter * exp(expon);
  return ret;
}

// Kuehn Santamaria Model

void Tau_Pion_Kaon::KS::SetModelParameters( GeneralModel _md ) 
{
}
 

Complex Tau_Pion_Kaon::KS::VectorFormFactor( double s )
{
  Complex ret(1.,0.);
  return ret;
}

Complex Tau_Pion_Kaon::KS::ScalarFormFactor( double s )
{
  Complex ret(1.,0.);
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
  Complex termK = Delta_KP/q2*(FS-FV)+FV;
  Complex termP = Delta_KP/q2*(FS-FV)-FV;
  for( int h=0; h<4; h++ ) {				// only nutau and tau have a helicity !
	ret += norm( 
		  F.X(m_nutau,m_kaon,0,h,cR,cL) * termK
		+ F.X(m_nutau,m_pion,0,h,cR,cL) * termP
		);
  }
  F.Delete();
  return ret/4.;
} 

double Tau_Pion_Kaon::operator()( const Vec4D *_p )
{
  double T = Using_Hels(_p);
  return T*0.5*GF2*Vus2;
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
  m_ms[0] = sqr( p_flavs[0].Mass() );
  // count number of pions, kaons and calc. mass^2
  for( int i=1; i<5; i++ ) {
	if( p_flavs[i].Kfcode() == kf::pi_plus )   nPion_ch++;
	if( p_flavs[i].Kfcode() == kf::pi )        nPion_0++;
	if( p_flavs[i].Kfcode() == kf::K_plus )	   nKaon_ch++;
	if( p_flavs[i].Kfcode() == kf::K )         nKaon_0++;
	m_ms[i] = sqr( p_flavs[i].Mass() );
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
	case   12 : /* K- K0b K0 */
				for( int i=1; i<5; i++ ) {
				  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
				  else {
					if( p_flavs[i].Kfcode() == kf::K ) {
					  if( p_flavs[i].IsAnti() == p_flavs[0].IsAnti() ) m_pseudo_3 = i;
					  else m_pseudo_1 = i;
					}
					else m_pseudo_2 = i;
				  }
				}
				break;
	case 1002 : /* K0 pi- K0b */
				for( int i=1; i<5; i++ ) {
				  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
				  else {
					if( p_flavs[i].Kfcode() == kf::K ) {
					  if( p_flavs[i].IsAnti() != p_flavs[0].IsAnti() ) m_pseudo_3 = i;
					  else m_pseudo_1 = i;
					}
					else m_pseudo_2 = i;
				  }
				}
				break;
	case  111 : /* K- pi0 K0 */
				for( int i=1; i<5; i++ ) {
				  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
				  else {
					if( p_flavs[i].Kfcode() == kf::K ) m_pseudo_3 = i;
					if( p_flavs[i].Kfcode() == kf::pi ) m_pseudo_2 = i;
					if( p_flavs[i].Kfcode() == kf::K_plus ) m_pseudo_1 = i;
				  }
				}
				break;
	case 1101 : /* pi- K0b pi0 */
				for( int i=1; i<5; i++ ) {
				  if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
				  else {
					if( p_flavs[i].Kfcode() == kf::pi ) m_pseudo_3 = i;
					if( p_flavs[i].Kfcode() == kf::K ) m_pseudo_2 = i;
					if( p_flavs[i].Kfcode() == kf::pi_plus ) m_pseudo_1 = i;
				  }
				}
				break;
  }
  // get internal number of a particle (in KS order)
  m_part[m_pseudo_1-1] = 1;
  m_part[m_pseudo_2-1] = 2;
  m_part[m_pseudo_3-1] = 3;
  m_part[m_nutau-1] = 4;
}
 
void Tau_Three_Pseudo::SetModelParameters( GeneralModel _md ) 
{ 
  Vud      = _md.pm.Vud;
  Vus      = _md.pm.Vus;
  fpi      = _md.pm.fpi;
  GF       = _md.pm.GF;
  cR = (0.,_md.pm.a - _md.pm.b);
  cL = (0.,_md.pm.a + _md.pm.b);
  SetA123();
  switch( _md.ff ) {
	case 2 : p_ff = new RChT();
			 break;
	case 1 : p_ff = new KS();
			 break;
  }
  p_ff->SetPath( m_path );
  p_ff->SetMode(m_mode);
  p_ff->SetOutgoingMasses2( m_ms[m_pseudo_1], m_ms[m_pseudo_2], m_ms[m_pseudo_3]  );
  p_ff->SetInternalNumbers( m_part );
  p_ff->SetModelParameters( _md );
}
 
void Tau_Three_Pseudo::SetA123()
{
  switch( m_mode ) {
	case 1200 : /* pi0 pi0 pi- mode */
				m_A123 = Vud;
				break;
	case   30 : /* K- K- K+ */
				m_A123 = Vus;
				break;
	case 2010 : /* K- pi- pi+ */
				m_A123 = -0.5*Vus;
				break;
	case 1020 : /* K- pi- K+ */
				m_A123 = -0.5*Vud;
				break;
	case 3000 : /* pi- pi- pi+ mode */
				m_A123 = Vud;
				break;
	case 1002 : /* K0 pi- K0b */
				m_A123 = -0.5*Vud;
				break;
	case  111 : /* K- pi0 K0 */
				m_A123 = 1.5*Vud*SQRT_05;
				break;
	case  210 : /* pi0 pi0 K- mode */
				m_A123 = Vus/4.;
				break;
	case 1101 : /* pi- K0b pi0 */
				m_A123 = 1.5*Vus*SQRT_05;
				break;
	case   12 : /* K- K0b K0 */
				m_A123 = -0.5*Vus;
				break;
  }
}

// Parameterisation
// DUMM, PICH, PORTOLES hep-ph/0312183 

void Tau_Three_Pseudo::RChT::SetModelParameters( GeneralModel _md ) 
{ 
  fpi  = _md.pm.fpi;
  l0   = _md.pm.lambda0;		   	// fit parameter lambda0
  MA   = _md.pm.MA;					// mass of axial resonance
  MV   = _md.pm.MV[0];				// mass of vector resonance
  MV2  = ATOOLS::sqr(MV);
  MA2  = ATOOLS::sqr(MA);
  gammaR = _md.pm.gammaR;			// global factor in off-shell rho width
  m    = ATOOLS::Flavour( ATOOLS::kf::pi_plus ).Mass();
  m2   = ATOOLS::sqr(m);
  mK2  = ATOOLS::sqr( ATOOLS::Flavour( ATOOLS::kf::K_plus ).Mass() );
  exp_alpha = _md.pm.exp_alpha;		// exponent in off-shell GA
  GA_at_MA2 = _md.pm.GA;			// on-shell axial width
  CreatePhiHistogram();
  Phi_at_MA2 = Phi( MA2 );
}
 
void Tau_Three_Pseudo::RChT::CreatePhiHistogram()
{
  // create file name
  char fn[100];
  sprintf(fn, "%sPhi(Q2)_MV=%.3f_fpi=%.4f_gR=%.2f.dat",m_path.c_str(),MV,fpi,gammaR);

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {							// if file exists
	// read table and create histogram
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::RChT::CreatePhiHistogram : \n"
	         <<"     Read phi(q2) for MV="<<MV<<" GeV from "<<fn<<"."<<endl;
	myHist = new Histogram( fn );
  }
  else {							// if file does not exist
	// create histogram (i.e. table of values)
	msg.Out()<<"Create necessery phase space function for choosen parameters.\n"
	         <<"This may take some time. Please wait..."<<endl;
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::RChT::CreatePhiHistogram : \n"
	         <<"     Create phi(q2) for MV="<<MV<<" GeV in "<<fn<<"."<<endl;
	double low (0.),
	up  (3.2);
	int nbins  (50);
	double step = (up-low)/nbins;
	myHist = new Histogram( 0, low, up+step, nbins+1 );
	double q2  (low),
	phi (0.);
	while( q2<=up+step ) {
	  phi = IntegralPhi(q2);				// get phi value
	  myHist->Insert( q2, phi );			// insert into histogram
	  q2 += step;		
	}
	myHist->Output(string(fn));
  }
  p_phi = myHist;
}

double Tau_Three_Pseudo::RChT::IntegralPhi( double Q2 )
{
  int Ns=500, Nt=500;					// number of subintervals
  double sum (0.);
  double s_max = sqr( sqrt(Q2)-m );
  double s_min = 4.*m2;
  double ds = (s_max-s_min)/Ns;
  double t_max (0.);
  double t_min (0.);
  double dt = (0.);
  double s (s_min), t, u;
  double V12, V22, V1V2;
  
  while ( s<s_max ) {
	t_max = ( sqr(Q2-m2) - sqr( Sqrt_Lambda(Q2, s, m2) - Sqrt_Lambda(s, m2, m2) ) )/(4.*s);
	t_min = ( sqr(Q2-m2) - sqr( Sqrt_Lambda(Q2, s, m2) + Sqrt_Lambda(s, m2, m2) ) )/(4.*s);
	dt = (t_max-t_min)/Nt;
	t = t_min;
	while ( t<t_max ) {
	  u = Q2 - s - t + 3.*m2;
	  V12  = 4.*m2-s - sqr(u-t)/4./Q2;
	  V22  = 4.*m2-t - sqr(u-s)/4./Q2;
	  V1V2 = (u-s-t+4.*m2)/2. - (u-t)*(u-s)/4./Q2;
	  Complex BW_s = Tools::BreitWigner( s, MV2, MassWidthVector(s) );
	  Complex BW_t = Tools::BreitWigner( s, MV2, MassWidthVector(s) );
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

double Tau_Three_Pseudo::RChT::MassWidthVector( double s )
{
  double MVGV (0.);
	if( s>4.*m2 )  MVGV += pow( 1.-4.*m2/s, 1.5 );
	if( s>4.*mK2 ) MVGV += pow( 1.-4.*mK2/s, 1.5 ) / 2.;
	MVGV *= gammaR*MV2*s/(96.*M_PI*sqr(fpi)); 
  return MVGV;
}

double Tau_Three_Pseudo::RChT::MassWidthAxial( double Q2 )
{
  double MAGA (0.);
  if( Q2 > 9.*m2 ) MAGA = MA * GA_at_MA2 * Phi(Q2) / Phi_at_MA2 * pow( MA2/Q2, exp_alpha );
  return MAGA;
}
 
Complex Tau_Three_Pseudo::RChT::FormFactor( int j, double Q2, double s, double t )
{
  if (j!=3) {
	double u = Q2-s-t+3.*m2;
	double x = (j==1)? s : t;
	double y = (j==1)? t : s;
	double MVGV_x = MassWidthVector(x);
	double MVGV_y = MassWidthVector(y);
	double MAGA = MassWidthAxial(Q2);
	double F_Q2_x = x/2./Q2 - l0*m2/Q2;
	double F_Q2_y = y/2./Q2 - l0*m2/Q2;
	Complex alpha = 1. - 3./2.*x/Complex(x-MV2, MVGV_x);
	Complex beta =  -3./2.*x/Complex(x-MV2, MVGV_x)
	  + F_Q2_x*(2.*Q2+x-u)/Complex(x-MV2, MVGV_x)
	  + F_Q2_y*(u-x)/Complex(y-MV2, MVGV_y);
	return SQRT_05*( alpha - Q2/Complex(Q2-MA2,MAGA)*beta );
		// the factor of SQRT_05 was added in order to yield a good expression for BR
  }
  else {
	return Complex(0.,0.);
  }
}

// Parameterisation
// DECKER, FINKEMEIER, MIRKES hep-ph/9310270

void Tau_Three_Pseudo::KS::SetModelParameters( GeneralModel _md ) 
{
   
  int a[2], b[2], c[2];
  if (!m_twoident) {
	for( int i=0; i<2; i++ ) {
	  a[i] = m_part[ _md.vec[i].i-1 ];
	  b[i] = m_part[ _md.vec[i].j-1 ];
	  if ( a[i]!=3 && b[i]!=3 ) {
		ATOOLS::msg.Error()<<ATOOLS::om::red
		  <<"ERROR in Tau_Three_Pseudo::DFM::SetModelParameters\n"
		  <<"     Resonances aren't set correctly.\n"
		  <<"     Make sure you have the right settings under \"Resonances\" vector"<<i+1<<" -> _ _.\n"
		  <<"     mode number : "<<m_mode<<"\n"
		  <<"     Don't know what to, will abort."
		  <<ATOOLS::om::reset<<std::endl;
		abort();
	  }
	  else c[i] = (a[i]==3)? b[i] : a[i];
	}
  }
  else {
	c[0] = c[1] = 1;
  }
  MA2    = ATOOLS::sqr(_md.pm.MA);			// mass^2 of axial resonance
  MAA2   = ATOOLS::sqr(_md.pm.MAA);			// mass^2 of axial resonance'
  msV[0] = ATOOLS::sqr(_md.pm.MV[c[0]-1]);	// mass^2 of vector resonance 13
  msV[1] = ATOOLS::sqr(_md.pm.MV[c[1]-1]);	// mass^2 of vector resonance 23
  msv[0] = ATOOLS::sqr(_md.pm.Mv[c[0]-1]);	// mass^2 of vector resonance' 13
  msv[1] = ATOOLS::sqr(_md.pm.Mv[c[1]-1]);	// mass^2 of vector resonance' 23
  Beta[0] = _md.pm.Beta[c[0]-1];				// weight factor for vector resonance' 13
  Beta[1] = _md.pm.Beta[c[1]-1];				// weight factor for vector resonance' 23
  alpha   = _md.pm.alpha;						// weight factor for axial resonance'
   
  GA      = _md.pm.GA;						// on-shell axial width
  GAA     = _md.pm.GAA;						// on-shell axial' width
  widthV[0] = _md.pm.GV[c[0]-1];				// on-shell vector 13 width
  widthV[1] = _md.pm.GV[c[1]-1];				// on-shell vector 23 width
  widthv[0] = _md.pm.Gv[c[0]-1];				// on-shell vector' 13 width
  widthv[1] = _md.pm.Gv[c[1]-1];				// on-shell vector' 23 width

  running_width = _md.run;

  CreateGHistogram();
} 
void Tau_Three_Pseudo::KS::SetMode( int m ) { 
  m_mode = m;
  m_twoident = false;
  if (m_mode==1200 ||
	  m_mode==210 ||
	  m_mode==30 ||
	  m_mode==3000 ||
	  m_mode==12) m_twoident = true;
  switch( m_mode ) {
	case 1200 : /* pi0 pi0 pi- mode */
	  m_X123  = ms[2];
	  m_ms123 = ms[2];
	  m_G123  = 1;
	  break;
	case   30 : /* K- K- K+ */
	  m_X123  = 2.*ms[2];
	  m_ms123 = ms[2];
	  m_G123  = 1;
	  break;
	case 2010 : /* K- pi- pi+ */
	  m_X123  = 2.*ms[0];
	  m_ms123 = ms[0];
	  m_G123  = 1;
	  break;
	case 1020 : /* K- pi- K+ */
	  m_X123  = ms[1] + ms[0];
	  m_ms123 = ms[1];
	  m_G123  = 1;
	  break;
	case 3000 : /* pi- pi- pi+ mode */
	  m_X123  = 2.*ms[0];
	  m_ms123 = ms[0];
	  m_G123  = 1;
	  break;
	case 1002 : /* K0 pi- K0b */
	  m_X123  = ms[1] + ms[0];
	  m_ms123 = ms[1];
	  m_G123  = 1;
	  break;
	case  111 : /* K- pi0 K0 */
	  m_X123  = 0.;
	  m_ms123 = ms[1];
	  m_G123  = 0;
	  break;
	case  210 : /* pi0 pi0 K- mode */
	  m_X123  = -2.*(ms[1]+ms[2]);
	  m_ms123 = ms[2];
	  m_G123  = 1;
	  break;
	case 1101 : /* pi- K0b pi0 */
	  m_X123  = 0.;
	  m_ms123 = ms[1];
	  m_G123  = 0;
	  break;
	case   12 : /* K- K0b K0 */
	  m_X123  = 2.*ms[0];
	  m_ms123 = ms[0];
	  m_G123  = 1;
	  break;
  }
}
 
void Tau_Three_Pseudo::KS::CreateGHistogram()
{
  // create file name
  char fn[100];
  if (m_G123) {
	sprintf(fn, "%sG(Q2)_MV13=%.3f_Mv13=%.3f_beta13=%.3f_MV23=%.3f_Mv23=%.3f_beta23=%.3f.dat",
		m_path.c_str(),
		sqrt(msV[0]), sqrt(msv[0]), Beta[0],
		sqrt(msV[1]), sqrt(msv[1]), Beta[1] );
  }
  else {
	sprintf(fn, "%sG(Q2)_noV13_MV23=%.3f_Mv23=%.3f_beta23=%.3f.dat",
		m_path.c_str(),
		sqrt(msV[1]), sqrt(msv[1]), Beta[1] );
  }

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {							// if file exists
	// read table and create histogram
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
	         <<"     Read phi(q2) for MV13="<<sqrt(msV[0])<<endl
			 <<"                      Mv13="<<sqrt(msv[0])<<endl
			 <<"                    beta13="<<Beta[0]<<endl
			 <<"                      MV23="<<sqrt(msV[1])<<endl
			 <<"                      Mv23="<<sqrt(msv[1])<<endl
			 <<"                    beta23="<<Beta[1]<<endl
			 <<"      from "<<fn<<"."<<endl;
	myHist = new Histogram( fn );
  }
  else {							// if file does not exist
	// create histogram (i.e. table of values)
	msg.Out()<<"Create necessery phase space function for choosen parameters.\n"
	         <<"This may take some time. Please wait..."<<endl;
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
	         <<"     Create phi(q2) for MV13="<<sqrt(msV[0])<<endl
			 <<"                       Mv13="<<sqrt(msv[0])<<endl
			 <<"                     beta13="<<Beta[0]<<endl
			 <<"                       MV23="<<sqrt(msV[1])<<endl
			 <<"                       Mv23="<<sqrt(msv[1])<<endl
			 <<"                     beta23="<<Beta[1]<<endl
			 <<"      in "<<fn<<"."<<endl;
	double low (0.),
	up  (3.2);
	int nbins  (50);
	double step = (up-low)/nbins;
	myHist = new Histogram( 0, low, up+step, nbins+1 );
	double q2  (low),
	phi (0.);
	while( q2<=up+step ) {
	  phi = IntegralG(q2);				// get phi value
	  myHist->Insert( q2, phi );		// insert into histogram
	  q2 += step;		
	}
	myHist->Output(string(fn));
  }
  p_G = myHist;
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
  ms[2] = ms[0];
  int Ns=500, Nt=500;					// number of subintervals
  double sum (0.);
  double s_max = Q2 + ms[2-1] - 2.*sqrt(Q2*ms[2-1]);
  double s_min = ms[1-1]+ms[3-1]+2.*sqrt(ms[1-1]*ms[3-1]);
  double ds = (s_max-s_min)/Ns;
  double t_max (0.);
  double t_min (0.);
  double dt (0.);
  double s (s_min), t, u;
  double V12, V22, V1V2;
  
  while ( s<s_max ) {
	t_max = (   sqr( Q2-ms[1-1]-ms[2-1]+ms[3-1] ) 
	          - sqr( Sqrt_Lambda(Q2, s, ms[2-1]) - Sqrt_Lambda(s,ms[1-1],ms[3-1]) ) 
			)/(4.*s);
	t_min = (   sqr( Q2-ms[1-1]-ms[2-1]+ms[3-1] ) 
	          - sqr( Sqrt_Lambda(Q2, s, ms[2-1]) + Sqrt_Lambda(s,ms[1-1],ms[3-1]) ) 
			)/(4.*s);
	dt = (t_max-t_min)/Nt;
	t = t_min;
	while ( t<t_max ) {
	  u = Q2 - s - t + ms[1-1]+ms[2-1]+ms[3-1];
	  V12  = 2.*(ms[1-1]+ms[2-1]) - s - sqr(u-t+ms[1-1]-ms[3-1])/(4.*Q2);
	  V22  = 2.*(ms[1-1]+ms[2-1]) - t - sqr(u-s+ms[2-1]-ms[3-1])/(4.*Q2);
	  V1V2 = (u-s-t+4.*ms[3-1])/2. - (u-t+ms[1-1]-ms[3-1])*(u-s+ms[2-1]-ms[3-1])/(4.*Q2);
	  Complex BW_1 = Tgen(1,2,3,s,t);
	  Complex BW_2 = Tgen(2,1,3,t,s);
	  sum += ( V12 * norm(BW_1) + V22*norm(BW_2) + 2.*V1V2*real( BW_1*conj(BW_2)) )*ds*dt;
	  t += dt;
	}
	s += ds;
  }
  return -1.*sum/Q2;
} 
 
Complex Tau_Three_Pseudo::KS::BW_A( double s )
{
  double MG_A (1.);
  double MG_AA (1.);
  if (running_width & 1) {
	MG_A  = sqrt(s)*GA*G(s)/G(MA2);
	MG_AA = sqrt(s)*GAA*G(s)/G(MAA2);
  }
  else {
	MG_A = sqrt(MA2)*GA;
	MG_AA = sqrt(MAA2)*GAA;
  }
  if (alpha!=0.) return ( Tools::BreitWigner(s,MA2,MG_A) + alpha*Tools::BreitWigner(s,MAA2,MG_AA) )/(1.+alpha);
  return Tools::BreitWigner(s,MA2,MG_A);
}


Complex Tau_Three_Pseudo::KS::Tvector1( int a, int b, double x )
{
  Complex ret =			  
              BW_V(a,b,x)*( 1.-(ms[a]-ms[b])/(3.*msV[a]) )
  + Beta[a]*( BW_v(a,b,x)*( 1.-(ms[a]-ms[b])/(3.*msv[a]) ) );
  return ret/(1.+Beta[a]);
}

Complex Tau_Three_Pseudo::KS::Tvector2( int a, int b, double x )
{
  Complex ret =
  			  BW_V(a,b,x)/msV[a] 
    + Beta[a]*BW_v(a,b,x)/msv[a];
  return ret*( 2.*(ms[a]-ms[b]) )/( 3.*(1.+Beta[a]) );	
}

Complex Tau_Three_Pseudo::KS::TSvector( int a, int b, int c, double Q2, double s, double t )
{
  Complex ret =
              BW_V(a,b,s) * (   Q2-2.*t-s+2.*ms[a]+ms[c] 
		                      - (ms[a]-ms[b])/msV[a]*( Q2+t-ms[c]-Q2*(t-msV[a]) ) )
	+ Beta[a]*BW_v(a,b,s) * (   Q2-2.*t-s+2.*ms[a]+ms[c] 
		                      - (ms[a]-ms[b])/msv[a]*( Q2+t-ms[c]-Q2*(t-msv[a]) ) );
  return ret/(1.+Beta[a]);	
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


Complex Tau_Three_Pseudo::KS::FormFactor( int j, double Q2, double s, double t )
{
  Complex FF(0.,0.);
  switch( j ) {
	case 1 : FF = BW_A(Q2)*Tgen(1,2,3,s,t);
			 break;
	case 2 : FF = BW_A(Q2)*Tgen(2,1,3,t,s); 
			 break;
	case 3 : FF = m_X123 + BW_A(Q2)*(Q2-MA2)*m_ms123/(MA2*Q2) 
	                     *( TSvector(1-1,3-1,2-1,Q2,s,t) + TSvector(2-1,3-1,1-1,Q2,t,s) );
			 FF /= 2.*(Q2-m_ms123);			 
			 break;
  }
  return FF;
}
 
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
  for( int h=0; h<4; h++ ) {	// helicity combination (nu,tau)
	ampl += norm(
	         F.X( m_nutau, m_pseudo_1, 0, h, cR, cL ) * ( FS + F1*(1.-d1) - F2*d2 )
	       + F.X( m_nutau, m_pseudo_2, 0, h, cR, cL ) * ( FS - F1*d1 + F2*(1.-d2) )
	       + F.X( m_nutau, m_pseudo_3, 0, h, cR, cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) )
		     ); 
  }
  F.Delete();
  return ampl/2.;
} // ME is invariant under exchange of two identical pseudos !


double Tau_Three_Pseudo::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*sqr(2.*GF*m_A123/fpi/3.); 
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  for( int i=1; i<6; i++ ) {
	m_ms[i] = sqr( p_flavs[i].Mass() );
  }
  for( int i=1; i<6; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) 	m_nutau = i;	// neutrino
	else {
	  if( p_flavs[i].Kfcode() == kf::pi ) 				m_pion0 = i;	// pi0
	  else {
		if( p_flavs[0].IsAnti()==p_flavs[i].IsAnti() )	m_pion3 = i;	// pi+
		else {
		  if( m_pion1==-1 ) 							m_pion1 = i;	// pi-
		  else 											m_pion2 = i;	// pi-
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
  Vud2     = sqr(_md.pm.Vud);
  Vus2     = sqr(_md.pm.Vus);
  GF2      = sqr(_md.pm.GF);
  cR       = Complex( 0., _md.pm.a-_md.pm.b );
  cL       = Complex( 0., _md.pm.a+_md.pm.b );

  Complex sum (0.,0.);
  for (int i=0; i<4; i++) {
	m_Alpha[i]  = _md.pm.Alpha[i];
	sum += m_Alpha[i];
  }
  m_SumAlpha = sum;


  p_lorenz = new KS();
  p_lorenz->SetModelParameters( _md );
}
 
void Tau_Four_Pion_3::LorenzBase::SetPrivates( Complex * _X, ATOOLS::Vec4D * _p ) 
{
  X  = _X;
  p  = _p;
  p[0] = p[1]+p[2]+p[3]+p[4];			// = q
  q2 = p[0].Abs2();
  for (int i=2; i<=4; i++ ) {
	r[i] = p[0]-p[i];
	s[i] = (p[1]+p[i]).Abs2();
  }
  // redefinition of variables
  s[0] = (p[2]+p[3]).Abs2();			// = t3
  s[1] = (p[2]+p[4]).Abs2();			// = t4
  X[0] = X[1]+X[2]+X[3]+X[4];			// X(nu,q,0)
  // unused variables yet
  r[0] = r[1] = ATOOLS::Vec4D(0.,0.,0.,0.);
};
 
// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

void Tau_Four_Pion_3::KS::SetModelParameters( GeneralModel _md )
{
  fpi2	  = sqr(_md.pm.fpi)*2.;		// redefine fpi
  grop    = _md.pm.grop;
  Go3p    = _md.pm.Go3p;

  mpi2 		= sqr( Flavour(kf::pi_plus).Mass() );
  mpi02 	= sqr( Flavour(kf::pi).Mass() );
  MR		= _md.pm.MR;
  MRR		= _md.pm.MRR;
  MRRR		= _md.pm.MRRR;
  MO      	= _md.pm.MO;
  MF        = _md.pm.Mdummy1;
  MS        = _md.pm.Mdummy2;
  MA		= _md.pm.MA;
   
  GR 		= _md.pm.GR;
  GRR 		= _md.pm.GRR;
  GRRR		= _md.pm.GRRR;
  GO      	= _md.pm.GO;
  GF        = _md.pm.Gdummy1;
  GS        = _md.pm.Gdummy2;
  GA		= _md.pm.GA;
   
  beta		= _md.pm.beta;
  gamma		= _md.pm.gamma;
  sigma     = _md.pm.sigma;
   
  MR2		= MR*MR;
  MRR2		= MRR*MRR;
  MRRR2		= MRRR*MRRR;
  MO2     	= MO*MO;
  MS2  		= MS*MS;
  MF2		= MF*MF;
  MA2		= MA*MA;
  Frho      = 0.266*MR2;

  R[0]=R[1] = 0.;
  R[2]      = -2.;
  R[3]=R[4] = 1.;

  for (int i=0; i<4; i++) {
	Beta_opi[i] = _md.pm.Beta_opi[i];
	Beta_api[i] = _md.pm.Beta_api[i];
	Beta_srh[i] = _md.pm.Beta_srh[i];
	Beta_frh[i] = _md.pm.Beta_frh[i];
  }
}

Complex Tau_Four_Pion_3::KS::Fk( double x, Complex * _beta )
{
//  Complex BW_R   = Tools::BreitWigner( x, MR2, GR, mpi2, 1. );
//  Complex BW_RR  = Tools::BreitWigner( x, MRR2, GRR, mpi2, 1. );
//  Complex BW_RRR = Tools::BreitWigner( x, MRRR2, GRRR, mpi2, 1. );
  Complex BW_R   = Tools::BreitWigner( x, MR2, sqrt(x)*GR );
  Complex BW_RR  = Tools::BreitWigner( x, MRR2, sqrt(x)*GRR );
  Complex BW_RRR = Tools::BreitWigner( x, MRRR2, sqrt(x)*GRRR );
  return( (_beta[0] + _beta[1]*BW_R + _beta[2]*BW_RR + _beta[3]*BW_RRR)/(_beta[0]+_beta[1]+_beta[2]+_beta[3]) );
}

Complex Tau_Four_Pion_3::KS::Trho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, MR2, GR, mpi2, 1. );
  Complex BW_RR  = Tools::BreitWigner( x, MRR2, GRR, mpi2, 1. );
  Complex BW_RRR = Tools::BreitWigner( x, MRRR2, GRRR, mpi2, 1. );
  return( (BW_R + beta*BW_RR + gamma*BW_RRR)/(1.+beta+gamma) );
}

Complex Tau_Four_Pion_3::KS::TTrho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, MR2, MR*GR )/MR2;
  Complex BW_RR  = Tools::BreitWigner( x, MRR2, MRR*GRR )/MRR2;
  return( BW_R + sigma*BW_RR );
}

double Tau_Four_Pion_3::KS::Dots( int k, int l )		// k=3,4; l=1,2,3,4 (l!=k)
  // special dot products used for anomalous part of J_omega
{
  int pre = l-1;				// predecessor of l
  if (pre==0) pre  = 4;
  if (pre==k) pre -= 1;		
  int suc  = l+1;				// successor of l
  if (suc==k)  suc += 1;
  if (suc==5)  suc  = 1;
  return( (r[k]*p[pre])*(p[suc]*p[k]) - (r[k]*p[suc])*(p[pre]*p[k]) );
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
	  sum2 += (X[0]-2.*X[l]) * (r[l]*(p[k]-p[1]))/r[l].Abs2();
	}
	sum1 += R[k]*Trho(s[k]) * (X[k]-X[1]-sum2);
  }
  J_chi = 2.*sqrt(3.)/fpi2*Trho(q2) * sum1;

  // anomalous part
  Complex J_a;
  sum1 = Complex(0.,0.);
  for( int k=3; k<=4; k++ ) {
	sum2 = Complex(0.,0.);
	for (int l=1; l<=4; l++) if (l!=k) {
	  sum2 += X[l] * Dots(k,l);
	}
	sum1 += Tools::BreitWigner( r[k].Abs2(), MO2, MO*GO )/MO2*sum2;
  }
  J_a = Go3p*Frho*grop*TTrho(q2) * sum1;

  // total
  return( (J_chi + J_a)*Fk(q2,Beta_opi) );
}

Complex Tau_Four_Pion_3::KS::AonePi()
{
  Complex term1 (0.,0.), term2 (0.,0.);
  Complex A, B, C;
  Vec4D   P, Q, R;

  //  1st term
  P = p[1]-p[3];
  Q = p[1]-p[4];
  R = r[2];
  A = Tools::BreitWigner(s[3],MR2,GR,mpi2,1.);
  B = Tools::BreitWigner(s[4],MR2,GR,mpi2,1.);
  C = Tools::BreitWigner(R.Abs2(),MA2,GA,mpi2,1.);
  term1  = A*(X[1]-X[3]) 
	+ B*(X[1]-X[4]) 
	- (X[0]-X[2])*( A*(R*P) + B*(R*Q) )/R.Abs2() 
	- X[0]*( A*(p[0]*P) + B*(p[0]*Q) - (p[0]*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() )/q2;
  term1 *= C;

  // 2nd and 3rd term
  int ind;
  Complex help;
  for (int k=3; k<=4; k++) {
	ind = (k==3)? 4 : 3;
	P = p[1]-p[2];
	Q = p[ind]-p[2];
	R = r[k];
	A = Tools::BreitWigner(s[2],MR2,GR,mpi02,mpi2,1.);
	B = Tools::BreitWigner(s[ind-3],MR2,GR,mpi02,mpi2,1.);	// -> s[0,1] = t[3,4]
	C = Tools::BreitWigner(R.Abs2(),MA2,GA,mpi2,1.);
	help   = A*(X[1]-X[2]) 
	  + B*(X[2]-X[ind]) 
	  - (X[0]-X[k])*( A*(R*P) + B*(R*Q) )/R.Abs2() 
	  - X[0]*( A*(p[0]*P) + B*(p[0]*Q) - (p[0]*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() )/q2;
	term2 += C*help;
  }

  // total 
  return (term1-term2)*Fk(q2,Beta_api);
}

Complex Tau_Four_Pion_3::KS::SigmaRho()
{
  int ind;
  Complex term (0.,0.);
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
	ind = (k==3)? 4 : 3;
	BW_S = Tools::BreitWigner(s[k],MS2,GS,mpi2,1.);
	BW_R = Tools::BreitWigner(s[ind-3],MR2,GR,mpi02,mpi2,1.);
	term += BW_S*BW_R * ( X[2] - X[ind] + X[0]*(p[0]*(p[ind]-p[2])) );
  }
  return term*Fk(q2,Beta_srh);
}

Complex Tau_Four_Pion_3::KS::FzeroRho()
{
  int ind;
  Complex term (0.,0.);
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
	ind = (k==3)? 4 : 3;
	BW_S = Tools::BreitWigner(s[k],MF2,GR,mpi2,1.);
	BW_R = Tools::BreitWigner(s[ind-3],MR2,GR,mpi02,mpi2,1.);
	term += BW_S*BW_R * ( X[2] - X[ind] + X[0]*(p[0]*(p[ind]-p[2])) );
  }
  return term*Fk(q2,Beta_frh);
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
	for (int h=0; h<4; h++) {			// helicity combination (nutau, tau)
	  for (int k=1; k<=4; k++) m_X[k] = F.X( m_nutau, m_inter[k], 0, h, cR, cL );
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
  return GF2/2.*Vud2 * T;
}
 
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  4 pion mode  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  for( int i=1; i<6; i++ ) {
	m_ms[i] = sqr( p_flavs[i].Mass() );
  }
  for( int i=1; i<6; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) 	m_nutau = i;	// neutrino
	else {
	  if( p_flavs[i].Kfcode() == kf::pi ) 				m_pion0 = i;	// pi0
	  else {
		if( p_flavs[0].IsAnti()==p_flavs[i].IsAnti() )	m_pion3 = i;	// pi+
		else {
		  if( m_pion1==-1 ) 							m_pion1 = i;	// pi-
		  else 											m_pion2 = i;	// pi-
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
  Vud2     = sqr(_md.pm.Vud);
  Vus2     = sqr(_md.pm.Vus);
  GF2      = sqr(_md.pm.GF);
  cR       = Complex( 0., _md.pm.a-_md.pm.b );
  cL       = Complex( 0., _md.pm.a+_md.pm.b );

  p_lorenz = new KS();
  p_lorenz->SetModelParameters( _md );
}
 
void Tau_Four_Pion_1::LorenzBase::SetPrivates( Complex * _X, ATOOLS::Vec4D * _p ) 
{
  X  = _X;
  p  = _p;
  p[0] = p[1]+p[2]+p[3]+p[4];			// = q
  q2 = p[0].Abs2();
  for (int i=2; i<=4; i++ ) {
	r[i] = p[0]-p[i];
	s[i] = (p[1]+p[i]).Abs2();
  }
  // redefinition of variables
  s[0] = (p[2]+p[3]).Abs2();			// = t3
  s[1] = (p[2]+p[4]).Abs2();			// = t4
  X[0] = X[1]+X[2]+X[3]+X[4];			// X(nu,q,0)
  // unused variables yet
  r[0] = r[1] = ATOOLS::Vec4D(0.,0.,0.,0.);
};
 
// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

void Tau_Four_Pion_1::KS::SetModelParameters( GeneralModel _md )
{
  fpi2	  = sqr(_md.pm.fpi)*2.;		// redefine fpi with factor sqrt(2.)

  mpi2 		= sqr( Flavour(kf::pi_plus).Mass() );
  mpi02 	= sqr( Flavour(kf::pi).Mass() );
  MR		= _md.pm.MR;
  MRR		= _md.pm.MRR;
  MRRR		= _md.pm.MRRR;
   
  GR 		= _md.pm.GR;
  GRR 		= _md.pm.GRR;
  GRRR		= _md.pm.GRRR;
   
  beta		= _md.pm.beta;
  gamma		= _md.pm.gamma;
   
  MR2		= MR*MR;
  MRR2		= MRR*MRR;
  MRRR2		= MRRR*MRRR;
}

Complex Tau_Four_Pion_1::KS::Trho( double x )
{
  Complex BW_R   = Tools::BreitWigner( x, MR2, GR, mpi2, 1. );
  Complex BW_RR  = Tools::BreitWigner( x, MRR2, GRR, mpi2, 1. );
  Complex BW_RRR = Tools::BreitWigner( x, MRRR2, GRRR, mpi2, 1. );
  return( (BW_R + beta*BW_RR + gamma*BW_RRR)/(1.+beta+gamma) );
}

Complex Tau_Four_Pion_1::KS::operator()()
{
  Complex J_chi;
  Complex sum1 (0.,0.), sum2 (0.,0.);
  sum1 = Complex(0.,0.);
  for (int k=2; k<=4; k++) {
	sum2 = Complex(0.,0.);
	for (int l=2; l<=4; l++) if (l!=k) {
	  sum2 += (X[0]-2.*X[l]) * (r[l]*(p[k]-p[1]))/r[l].Abs2();
	}
	sum1 += Trho(s[k]) * (X[k]-X[1]-sum2);
  }
  J_chi = 2.*sqrt(3.)/fpi2*Trho(q2) * sum1;

  // total
  return( J_chi );
}

// General framework

double Tau_Four_Pion_1::Using_Hels( const Vec4D * _p )
{
  // internal numeration and convenient variables
  for (int i=1; i<=5; i++ ) m_p[i] = _p[m_inter[i]];

  // summation over helicities
  XYZFunc F(m_nout,_p,p_flavs);
  double ampl (0.);
  Complex help(0.,0.);
	for (int h=0; h<4; h++) {			// helicity combination (nutau, tau)
	  for (int k=1; k<=4; k++) m_X[k] = F.X( m_nutau, m_inter[k], 0, h, cR, cL );
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
  return GF2/2.*Vud2 * T;
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
  m_pion_1(-1),
  m_pion_2(-1)
{
  m_metype = string("Tau_EtaTwoPion");
  // calc. mass^2
  for( int i=0; i<5; i++ ) {
	m_ms[i] = sqr( p_flavs[i].Mass() );
  }
}
 
void Tau_Eta_Two_Pion::SetModelParameters( GeneralModel _md ) 
{ 
  a        = _md.pm.a;
  b        = _md.pm.b;
  Vud      = _md.pm.Vud;
  fpi      = _md.pm.fpi;
  GF       = _md.pm.GF;
  cR = (0.,a - b);
  cL = (0.,a + b);
}

double Tau_Eta_Two_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  Vec4D p1( _p[m_pion_1] ),
        p2( _p[m_pion_2] ),
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
  for( int h=0; h<4; h++ ) {	// helicity combination (nu,tau)
	ampl += norm(
	         F.X( m_nutau, m_pion_1, 0, h, cR, cL ) * ( FS + F1*(1.-d1) - F2*d2 )
	       + F.X( m_nutau, m_pion_2, 0, h, cR, cL ) * ( FS - F1*d1 + F2*(1.-d2) )
	       + F.X( m_nutau, m_eta, 0, h, cR, cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) )
		     ); 
  }
  F.Delete();
  return ampl/2.;
} 


double Tau_Eta_Two_Pion::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*sqr(2.*GF*Vud/fpi/3.); 
}
 
