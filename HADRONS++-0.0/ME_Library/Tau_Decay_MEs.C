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
  double T(1.);
  double q2 = (_p[0]-_p[m_nutau]).Abs2();
  switch( m_md.me ) {
	case MD_TRACES			: T = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: T = Using_Traces_Simple(_p); break;
	case MD_XYZ				: T = Using_Hels(_p); break;						  
  }
  return T*sqr(GF);
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
  for( int h=0; h<4; h++ ) {
	ret += norm( F.X(m_nutau,m_pion,0,h,cR,cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Pion::operator()( const Vec4D *_p )
{
  double T(1.);
  double q2 = _p[m_pion].Abs2(); 
  switch( m_md.me ){
	case MD_TRACES			: T = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: T = Using_Traces_Simple(_p); break;
	case MD_XYZ				: T = Using_Hels(_p); break;
  }
  return T*0.5*sqr(GF*Vxx*fxx); 
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
	  if( p_flavs[i] == kf::pi || p_flavs[i] == kf::K ) m_pion0 = i;
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
	double expon = -1.*s/(96.*sqr(M_PI*fxx))*AA.real();
	if(m_md.run) MRGR = -1.*MR2*s/(96.*sqr(M_PI*fxx)) * AA.imag();
	else  MRGR = MR*GR;
	Complex denom = Complex(s-MR2,MRGR);
	ret = MR2/denom * exp(expon);
  }
  return ret;
}

double Tau_Two_Pion::Using_Traces_Simple( const Vec4D *_p )
{
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
  for( int h=0; h<4; h++ ) {
	ret += norm( F.X(m_nutau,m_pion_ch,0,h,cR,cL)
				 - F.X(m_nutau,m_pion0,0,h,cR,cL) );
  }
  F.Delete();
  return ret;
} // its value is 100% identical to simple traces !

double Tau_Two_Pion::operator()( const Vec4D *_p )
{
  double T(1.);
  switch( m_md.me ) {
	case MD_TRACES			: T = Using_Traces(_p); break;
	case MD_TRACES_SIMPLE	: T = Using_Traces_Simple(_p); break;
	case MD_XYZ				: T = Using_Hels(_p); break; 
  }
  double q2 = (_p[m_pion_ch] + _p[m_pion0] ).Abs2();
  Complex FF = FormFactor(q2);
  return T*0.5*sqr(GF*Vud)*norm(FF)*CG;
}


/////////////////////////////////////////////////////////
///  3 pion mode  ///////////////////////////////////////
/////////////////////////////////////////////////////////

Tau_Three_Pion::Tau_Three_Pion( int _nout, Flavour *_fl ) :
  HD_ME_Base(_nout,_fl),
  m_nutau(-1),
  m_pion_ch(-1),
  m_pion_1(-1),
  m_pion_2(-1)
{
  m_metype = string("Tau_ThreePion");
  int n_charged(0);
  for( int i=1; i<5; i++ ) {
	cout<<p_flavs[i]<<" ";
	if( p_flavs[i].Kfcode() == kf::pi_plus ) n_charged++;
  }
  for( int i=1; i<5; i++ ) {
	if( p_flavs[i].Kfcode() == p_flavs[0].Kfcode()+1 ) m_nutau = i;
	else {
	  switch( n_charged ) {
		case 1:	/* pi- pi0 pi0 mode ==> M(-) */
		  if( p_flavs[i].Kfcode() == kf::pi_plus ) m_pion_ch = i;
		  else if( m_pion_1<0 ) m_pion_1 = i;
		       else m_pion_2 = i;
		  break;
		case 3: /* pi+ pi- pi- mode ==> M(+) */  
		  if ( p_flavs[i].IsAnti() == p_flavs[0].IsAnti() ) m_pion_ch = i;
		  else if( m_pion_1<0 ) m_pion_1 = i;
		       else m_pion_2 = i;
		  break;
	  }
	}
  }
}

double Tau_Three_Pion::Sqrt_Lambda( double _a, double _b, double _c )
{
  return sqrt(sqr(_a+_b-_c)-4.*_a*_b);
}

double Tau_Three_Pion::MassWidthVector( double s )
{
  double MVGV (0.);
	if( s>4.*m2 )  MVGV += pow( 1.-4.*m2/s, 1.5 );
	if( s>4.*mK2 ) MVGV += pow( 1.-4.*mK2/s, 1.5 ) / 2.;
	MVGV *= MV2*s/(96.*M_PI*sqr(fpi)); 
  return MVGV;
}

Complex Tau_Three_Pion::BreitWignerRho( double s )
{
  return Complex( MV2, 0.) / Complex( MV2-s, -1.*MassWidthVector(s) );
}

double Tau_Three_Pion::FitOrder19oddPhi( double Q2 )
{
  double ca (7.98935e-6);
  double cb (-0.000393075);
  double cc (0.00823104);
  double cd (-0.0955909);
  double ce (0.672044);
  double cf (-2.92582);
  double cg (7.71015);
  double ch (-11.2143);
  double ci (3.52726);
  double cj (-0.208412);
  return ca*pow(Q2,19)
	+ cb*pow(Q2,17)
	+ cc*pow(Q2,15)
	+ cd*pow(Q2,13)
	+ ce*pow(Q2,11)
	+ cf*pow(Q2,9)
	+ cg*pow(Q2,7)
	+ ch*pow(Q2,5)
	+ ci*pow(Q2,3)
	+ cj*Q2;
}

double Tau_Three_Pion::FitOrder5Phi( double Q2 )
{
  double ca (-0.00097277);
  double cb (-0.195621);
  double cc (0.0982911);
  double cd (-0.0367806);
  double ce (0.00799563);
  double cf (-0.000377053); 
  return ca*pow(Q2,5)
	+ cb*pow(Q2,4)
	+ cc*pow(Q2,3)
	+ cd*pow(Q2,2)
	+ ce*Q2 + cf;
}

double Tau_Three_Pion::IntegralPhi( double Q2 )
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
	  V1V2 = (u-s-t+3.*m2)/2. - (u-t)*(u-s)/4./Q2;
	  Complex BW_s = BreitWignerRho(s);
	  Complex BW_t = BreitWignerRho(t);
	  sum += ( V12 * norm(BW_s) + V22*norm(BW_t) + 2.*V1V2*real( BW_s*conj(BW_t)) )*ds*dt;
	  t += dt;
	}
	s += ds;
  }
  return sum*Q2;
}

double Tau_Three_Pion::Phi( double x )
{ 
  return FitOrder19oddPhi(x);
//   return FitOrder5Phi(x);
//  return IntegralPhi(x);
}

double Tau_Three_Pion::MassWidthAxial( double Q2 )
{
  double MAGA (0.);
  if( Q2 > 9.*m2 ) MAGA = MA * GA_at_MA2 * Phi(Q2) / Phi_at_MA2 * pow( MA2/Q2, exp_alpha );
//   if( Q2 > 9.*m2 ) MAGA = MA * GA_at_MA2 * Phi(Q2) / Phi(MA2) * pow( MA2/Q2, exp_alpha );
  return MAGA;
}
Complex Tau_Three_Pion::FormFactor( int j, double Q2, double s, double t, double u )
{
  double MVGV = MassWidthVector(s);
  double MAGA = MassWidthAxial(Q2);
  ofstream f("width.dat",ios::out|ios::app);
  f<<Q2<<" "<<MAGA/MA<<endl;
  double x = (j==1)? s : t;
  double y = (j==1)? t : s;
  double F_Q2_x = x/2./Q2 - l0*m2/Q2;
  double F_Q2_y = y/2./Q2 - l0*m2/Q2;
  Complex alpha = 1. - 3./2.*x/Complex(x-MV2, MVGV);
  Complex beta =  -3./2.*x/Complex(x-MV2, MVGV)
	            + F_Q2_x*(2.*Q2+x-u)/Complex(x-MV2, MVGV)
	            + F_Q2_y*(u-x)/Complex(t-MV2, MVGV);
  return alpha - Q2/Complex(Q2-MA2,MAGA)*beta;
}

double Tau_Three_Pion::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
  Vec4D p1( _p[m_pion_1] ),
        p2( _p[m_pion_2] ),
        p3( _p[m_pion_ch] );
  Vec4D Q( p1+p2+p3 );		
  double s = (p1+p3).Abs2(),
         t = (p2+p3).Abs2(),
         u = (p1+p2).Abs2();
  double Q2 = Q.Abs2();
  double dot1 = Q*(p1-p3),
         dot2 = Q*(p2-p3);
  Complex F1 = FormFactor( 1, Q2, s, t, u );		 
  Complex F2 = FormFactor( 2, Q2, s, t, u );		 
  Complex ampl(0.,0.);
  for( int h=0; h<4; h++ ) {
	ampl +=  F.X( m_nutau, m_pion_1, 0, h, cR, cL ) * ( F1*(1.-dot1/Q2)-F2 )
	       + F.X( m_nutau, m_pion_2, 0, h, cR, cL ) * ( F2*(1.-dot2/Q2)-F1 )
	       - F.X( m_nutau, m_pion_ch,0, h, cR, cL ) * ( F1*(1.+dot1/Q2)+F2*(1.+dot2/Q2) ); 
  }
  F.Delete();
  return norm( ampl );
} 


double Tau_Three_Pion::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*sqr(2.*GF*Vud/fpi/3.); 
}
