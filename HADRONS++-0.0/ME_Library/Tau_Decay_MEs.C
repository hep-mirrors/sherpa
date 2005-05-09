#include "Tau_Decay_MEs.H"
#include "Message.H"
#include "Traces.H"
#include "XYZFuncs.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;
using namespace HADRONS;

/////////////////////////////////////////////////////////
///  leptonic decay  ////////////////////////////////////
/////////////////////////////////////////////////////////

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

double Tau_Lepton::Using_Traces( const Vec4D *_p )
{
  double M = p_masses[0];
  double m = p_masses[m_lep]; 
  double a=m_md.pm.a, b=m_md.pm.b;
  double a2=m_md.pm.a2, b2=m_md.pm.b2;
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

double Tau_Lepton::Using_Traces_Simple( const Vec4D *_p )
{
  double M = p_masses[0];
  double m = p_masses[m_lep]; 
  double a=m_md.pm.a, b=m_md.pm.b;
  double a2=m_md.pm.a2, b2=m_md.pm.b2;
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

double Tau_Lepton::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  Complex cR1 = (0.,m_md.pm.a-m_md.pm.b);
  Complex cL1 = (0.,m_md.pm.a+m_md.pm.b);
  Complex cR2 = (0.,m_md.pm.a2-m_md.pm.b2);
  Complex cL2 = (0.,m_md.pm.a2+m_md.pm.b2);
  double MW2 = sqr(Flavour(kf::W).Mass());
  for( int h1=0; h1<4; h1++ ) for( int h2=0; h2<4; h2++ ) {
	int h = (h1<<2)+h2;
	ret += norm(    F.Z( m_nutau, 0, m_lep, m_nulep, h, cR1, cL1, cR2, cL2 )
				 - (F.X(m_nutau,0,0,h1,cR1,cL1)-F.X(m_nutau,m_nutau,0,h1,cR1,cL1))*
				   (F.X(m_lep,0,m_nulep,h2,cR2,cL2)-F.X(m_lep,m_nutau,m_nulep,h2,cR2,cL2))/MW2 );
  }
  F.Delete();
  return ret*0.25;
} // its value is 100% identical to simple traces

double Tau_Lepton::operator()( const Vec4D *_p )
{
  double a(1.);
  double q2 = (_p[0]-_p[m_nutau]).Abs2();
  double MW2 = sqr(Flavour(kf::W).Mass());
  double GW2 = sqr(Flavour(kf::W).Width());
  double GF = m_md.pm.GF;
  switch( m_md.me ) {
	case MD_TRACES			: a = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: a = Using_Traces_Simple(_p); break;
	case MD_XYZ				: a = Using_Hels(_p); break;						  
  }
  return a*sqr(GF*MW2)/(sqr(q2-MW2)+MW2*GW2);
}

/////////////////////////////////////////////////////////
///  1 pion/kaon mode  //////////////////////////////////
/////////////////////////////////////////////////////////

Tau_Pion::Tau_Pion( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion(-1)
{
  m_metype = string("Tau_Pion");
  for( int i=1; i<3; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
	else m_pion = i;							// that's the pion
  }
}

double Tau_Pion::Using_Traces_Simple( const Vec4D *_p )
{
  double a=m_md.pm.a, b=m_md.pm.b;
  Vec4D P = _p[0];
  Vec4D k = _p[m_nutau];
  Vec4D p = _p[m_pion];
  double M2=0.;
  double term1(0.);
  
  // | VA structure |^2
  term1 = Trace(k,p,P,p);
  M2  = (sqr(a)+sqr(b))*term1; 

  return M2;
}

double Tau_Pion::Using_Traces( const Vec4D *_p )
{
  double M = p_masses[0];
  double a=m_md.pm.a, b=m_md.pm.b;
  Vec4D P = _p[0];
  Vec4D k = _p[m_nutau];
  Vec4D p = _p[m_pion];
  Vec4D S(0.,0.,0.,1.);
  double M2=0.;
  double term1(0.), term2(0.);

  // | VA structure |^2
  term1 = Trace(k,p,P,p);
  term2 = Trace(k,p,S,p);
  M2  = (sqr(a)+sqr(b))*term1 - 2.*a*b*M*term2; 

  return M2;
}

double Tau_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  Complex cR = (0.,m_md.pm.a-m_md.pm.b);
  Complex cL = (0.,m_md.pm.a+m_md.pm.b);
  for( int h=0; h<4; h++ ) {
	ret += norm( F.X(m_nutau,m_pion,0,h,cR,cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Pion::operator()( const Vec4D *_p )
{
  double a(1.);
  double fxx = (p_flavs[m_pion].Kfcode() == kf::pi_plus )? m_md.pm.fpi : m_md.pm.fK;
  double Vxx = (p_flavs[m_pion].Kfcode() == kf::pi_plus )? m_md.pm.Vud : m_md.pm.Vus;
  double MW2 = sqr(Flavour(kf::W).Mass());
  double GW2 = sqr(Flavour(kf::W).Width());
  double q2 = _p[m_pion].Abs2(); 
  double GF = m_md.pm.GF;
  switch( m_md.me ){
	case MD_TRACES			: a = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: a = Using_Traces_Simple(_p); break;
	case MD_XYZ				: a = Using_Hels(_p)*sqr(1.-p_masses2[m_pion]/MW2); break;
  }
  return a*0.5*sqr(GF*MW2*Vxx*fxx*m_md.pm.b2)/(sqr(q2-MW2)+MW2*GW2); 
}


/////////////////////////////////////////////////////////
///  2 pion mode  ///////////////////////////////////////
/////////////////////////////////////////////////////////

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
	  if( p_flavs[i] == kf::pi ) m_pion0 = i;
	  else m_pion_ch = i;
	}
  }
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
  double MR = Flavour(kf::rho_770).Mass();
  double GR = Flavour(kf::rho_770).Width();
  double MR2 = sqr(MR);
  double MRGR = 0.;
  double m = p_masses[m_pion_ch];
  double frho(m_md.pm.frho), grpp(m_md.pm.grpp), fpi(m_md.pm.fpi);
  Complex ret(1.,0.);
  if( m_md.ff == 1 ) {			// Breit-Wigner-rho
	if( m_md.run ) MRGR = sqrt(s)*GR*MR2/s*pow((s-4.*sqr(m))/(MR2-4.*sqr(m)), 1.5);
	else MRGR = MR*GR;
	Complex denom = Complex(s-MR2,MRGR); 
	ret = MR*frho*grpp/denom; 
  }
  if( m_md.ff == 2 ) {			// Effective Field Theory 
	double m2_pi = sqr( Flavour(kf::pi_plus).Mass() );
	double m2_K  = sqr( Flavour(kf::K_plus).Mass() );
	Complex AA = A( m2_pi/s, m2_pi/MR2 ) + 0.5*A( m2_K/s, m2_K/MR2 );
	double expon = -1.*s/(96.*sqr(M_PI*fpi))*AA.real();
	if(m_md.run) MRGR = -1.*MR2*s/(96.*sqr(M_PI*fpi)) * AA.imag();
	else  MRGR = MR*GR;
	Complex denom = Complex(s-MR2,MRGR);
	ret = MR2/denom * exp(expon);
  }
  return ret;
}

double Tau_Two_Pion::Using_Traces_Simple( const Vec4D *_p )
{
  double a=m_md.pm.a, b=m_md.pm.b;
  Vec4D P = _p[0];
  Vec4D k = _p[m_nutau];
  Vec4D p2 = _p[m_pion_ch];
  Vec4D p3 = _p[m_pion0];
  double dot = p_masses2[m_pion_ch] - p_masses2[m_pion0];
  double MW2 = sqr(Flavour(kf::W).Mass());

  double M2=0.;
  double term1(0.);

  // | VA structure |^2
  double M2_part(0.);
	term1 = Trace(k,p2,P,p2); 
  M2_part += (sqr(a)+sqr(b))*term1 * sqr( 1.-dot/MW2 ); 
	term1 = Trace(k,p3,P,p3);
  M2_part += (sqr(a)+sqr(b))*term1 * sqr( 1.+dot/MW2 ); 
	term1 = Trace(k,p2,P,p3) + Trace(k,p3,P,p2);
  M2_part -= (sqr(a)+sqr(b))*term1 * ( 1.-sqr(dot/MW2) ); 

  M2 = M2_part;
  
  return M2;
}

double Tau_Two_Pion::Using_Traces( const Vec4D *_p )
{
  double M = p_masses[0];
  double a=m_md.pm.a, b=m_md.pm.b;
  Vec4D P = _p[0];
  Vec4D k = _p[m_nutau];
  Vec4D p2 = _p[m_pion_ch];
  Vec4D p3 = _p[m_pion0];
  Vec4D S(0.,0.,0.,1.);
  double q2 = (P-k).Abs2();
  double dot = p_masses2[m_pion_ch] - p_masses2[m_pion0];
  double MW2 = sqr(Flavour(kf::W).Mass());

  double M2=0.;
  double term1(0.), term2(0.);

  // | VA structure |^2
  double M2_part(0.);
	term1 = Trace(k,p2,P,p2); 
	term2 = Trace(k,p2,S,p2);
  M2_part += ( (sqr(a)+sqr(b))*term1 - 2.*a*b*M*term2 ) * sqr( 1.-dot/MW2 ); 
	term1 = Trace(k,p3,P,p3);
	term2 = Trace(k,p3,S,p3);
  M2_part += ( (sqr(a)+sqr(b))*term1 - 2.*a*b*M*term2 ) * sqr( 1.+dot/MW2 ); 
	term1 = Trace(k,p2,P,p3) + Trace(k,p3,P,p2);
	term2 = Trace(k,p2,S,p3) + Trace(k,p3,S,p2);
  M2_part -= ( (sqr(a)+sqr(b))*term1 - 2.*a*b*M*term2 ) * ( 1.-sqr(dot/MW2) ); 

  M2 = M2_part;
  
  return M2; 
}

double Tau_Two_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  Complex cR = (0.,m_md.pm.a - m_md.pm.b);
  Complex cL = (0.,m_md.pm.a + m_md.pm.b);
  double MW2 = sqr(Flavour(kf::W).Mass());
  double dot = p_masses2[m_pion_ch] - p_masses2[m_pion0];
  for( int h=0; h<4; h++ ) {
	ret += norm( F.X(m_nutau,m_pion_ch,0,h,cR,cL)*(1.-dot/MW2)
				 - F.X(m_nutau,m_pion0,0,h,cR,cL)*(1.+dot/MW2) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Two_Pion::operator()( const Vec4D *_p )
{
  double a(1.);
  double MW2 = sqr(Flavour(kf::W).Mass());
  double GW2 = sqr(Flavour(kf::W).Width());
  double GF = m_md.pm.GF;
  double Vud = m_md.pm.Vud;
  switch( m_md.me ) {
	case MD_TRACES			: a = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: a = Using_Traces_Simple(_p); break;
	case MD_XYZ				: a = Using_Hels(_p); break; 
  }
  double q2 = (_p[m_pion_ch] + _p[m_pion0] ).Abs2();
  Complex FF = FormFactor(q2);
  return a*0.5*sqr(GF*MW2*Vud*m_md.pm.a2)/(sqr(q2-MW2)+MW2*GW2)*norm(FF);
}

