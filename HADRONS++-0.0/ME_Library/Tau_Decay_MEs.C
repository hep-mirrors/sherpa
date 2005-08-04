#include "Tau_Decay_MEs.H"
#include "Message.H"
#include "Traces.H"
#include "XYZFuncs.H"
#include "Histogram.H"
#include "Channel_Basics.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

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

double Tau_Pseudo::Using_Traces_Simple( const Vec4D *_p )
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

double Tau_Pseudo::Using_Traces( const Vec4D *_p )
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
  double MRGR = 0.;

  double m = p_masses[m_pion_ch];
  Complex ret(1.,0.);
  if( m_md.ff == 1 ) {			// Breit-Wigner-rho
	Complex BWo = MO2/Complex(MO2-s, -1.*sqrt(s)*GO);
	double Gamma_r = GR;
	double Gamma_rr = GRR;
	double Gamma_rrr = GRRR;
	if (m_md.run) {
	  Gamma_r = (s<4.*MP2)? 0. : GR*pow( MR2/s, lambda )*pow( (s-4.*MP2)/(MR2-4.*MP2), 1.5 );
	  Gamma_rr = (s<4.*MP2)? 0. : GRR*pow( MRR2/s, lambda )*pow( (s-4.*MP2)/(MRR2-4.*MP2), 1.5 );
	  Gamma_rrr = (s<4.*MP2)? 0. : GRRR*pow( MRRR2/s, lambda )*pow( (s-4.*MP2)/(MRRR2-4.*MP2), 1.5 );
	}
	Complex BWr = MR2/Complex(MR2-s,-1.*sqrt(s)*Gamma_r);
	Complex BWrr = MRR2/Complex(MRR2-s,-1.*sqrt(s)*Gamma_rr);
	Complex BWrrr = MRRR2/Complex(MRRR2-s,-1.*sqrt(s)*Gamma_rrr);
	Complex help = BWr * (1.+delta*BWo)/(1.+delta); 
	ret = ( help + Complex(beta,0.)*BWrr + Complex(gamma,0.)*BWrrr )/( 1.+beta+gamma );
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
///  3 pseudo mode  /////////////////////////////////////
/////////////////////////////////////////////////////////

double Tau_Three_Pseudo::DPP::Sqrt_Lambda( double _a, double _b, double _c )
{
//  return PHASIC::Channel_Basics::SqLam(_a,_b,_c);
  double L = sqr(_a-_b-_c) - 4.*_b*_c;
  if (L>=0.) return sqrt(L);
  return 0.;
}
 
double Tau_Three_Pseudo::DFM::Sqrt_Lambda( double _a, double _b, double _c )
{
  double L = sqr(_a-_b-_c) - 4.*_b*_c;
  if (L>=0.) return sqrt(L);
  return 0.;
}
 
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
  // get internal number of a particle (in DFM order)
  m_part[m_pseudo_1-1] = 1;
  m_part[m_pseudo_2-1] = 2;
  m_part[m_pseudo_3-1] = 3;
  m_part[m_nutau-1] = 4;
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

void Tau_Three_Pseudo::DPP::CreatePhiHistogram()
{
  // create file name
  char fn[100];
  sprintf(fn, "%sPhi(Q2)_MV=%.3f_fpi=%.4f_gR=%.2f.dat",m_path.c_str(),MV,fpi,gammaR);

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {							// if file exists
	// read table and create histogram
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::DPP::CreatePhiHistogram : \n"
	         <<"     Read phi(q2) for MV="<<MV<<" GeV from "<<fn<<"."<<endl;
	myHist = new Histogram( fn );
  }
  else {							// if file does not exist
	// create histogram (i.e. table of values)
	msg.Out()<<"Create necessery phase space function for choosen parameters.\n"
	         <<"This may take some time. Please wait..."<<endl;
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::DPP::CreatePhiHistogram : \n"
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

double Tau_Three_Pseudo::DPP::MassWidthVector( double s )
{
  double MVGV (0.);
	if( s>4.*m2 )  MVGV += pow( 1.-4.*m2/s, 1.5 );
	if( s>4.*mK2 ) MVGV += pow( 1.-4.*mK2/s, 1.5 ) / 2.;
	MVGV *= gammaR*MV2*s/(96.*M_PI*sqr(fpi)); 
  return MVGV;
}

Complex Tau_Three_Pseudo::DPP::BreitWignerRho( double s )
{
  return Complex( MV2, 0.) / Complex( MV2-s, -1.*MassWidthVector(s) );
}

double Tau_Three_Pseudo::DPP::FitOrder7oddLowPhi( double x )
{
  double ca (0.00124087);
  double cb (-0.0244092);
  double cc (0.170489);
  double cd (-2.71711);
  return( ca*x + cb*pow(x,3) + cc*pow(x,5) + cd*pow(x,7) );
}

double Tau_Three_Pseudo::DPP::FitOrder4HighPhi( double x )
{
  double ca (-2.56323);
  double cb (11.5306);
  double cc (-11.827);
  double cd (0.425535);
  double ce (-0.200207);
  return( ca + cb*x + cc*pow(x,2) + cd*pow(x,3) + ce*pow(x,4) );
}
 
double Tau_Three_Pseudo::DPP::FitOrder4InterPhi( double x )
{
  double ca (27.3171);
  double cb (-152.776);
  double cc (313.397);
  double cd (-276.024);
  double ce (85.4819);
  return( ca + cb*x + cc*pow(x,2) + cd*pow(x,3) + ce*pow(x,4) );
}

double Tau_Three_Pseudo::DPP::OptimisedPhi( double x )
{
  double p (0.);
  if( x<0.6 ) p = FitOrder7oddLowPhi(x);	// low region fit
  else {
	if (x<1.) p = FitOrder4InterPhi(x);		// intermediate region fit
	else      p = FitOrder4HighPhi(x);		// high region fit
  }
  return p;
}

double Tau_Three_Pseudo::DPP::IntegralPhi( double Q2 )
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
	  Complex BW_s = BreitWignerRho(s);
	  Complex BW_t = BreitWignerRho(t);
	  sum += ( V12 * norm(BW_s) + V22*norm(BW_t) + 2.*V1V2*real( BW_s*conj(BW_t)) )*ds*dt;
	  t += dt;
	}
	s += ds;
  }
  return sum*Q2;
}

double Tau_Three_Pseudo::DPP::Phi( double x )
{ 
  double val;	
  p_phi->Extrapolate( x+p_phi->BinSize(), &val, 0 );
  							// shift +BinSize to get right value
  double p = val;
  return p;
}

double Tau_Three_Pseudo::DPP::MassWidthAxial( double Q2 )
{
  double MAGA (0.);
  if( Q2 > 9.*m2 ) MAGA = MA * GA_at_MA2 * Phi(Q2) / Phi_at_MA2 * pow( MA2/Q2, exp_alpha );
  return MAGA;
}
 
Complex Tau_Three_Pseudo::DPP::FormFactor( int j, double Q2, double s, double t )
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
	return alpha - Q2/Complex(Q2-MA2,MAGA)*beta;
  }
  else {
	return Complex(0.,0.);
  }
}

// Parameterisation
// DECKER, FINKEMEIER, MIRKES hep-ph/9310270

void Tau_Three_Pseudo::DFM::CreateGHistogram()
{
  // create file name
  char fn[100];
  sprintf(fn, "%sG(Q2)_MV13=%.3f_Mv13=%.3f_beta13=%.3f_MV23=%.3f_Mv23=%.3f_beta23=%.3f.dat",
	  m_path.c_str(),
	  sqrt(msV[0]), sqrt(msv[0]), Beta[0],
	  sqrt(msV[1]), sqrt(msv[1]), Beta[1] );

  // look if file already exists
  ifstream f( fn );
  Histogram * myHist;
  if (f) {							// if file exists
	// read table and create histogram
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::DFM::CreateGHistogram : \n"
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
	msg.Tracking()<<"HADRONS::Tau_Three_Pseudo::DFM::CreateGHistogram : \n"
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
	  myHist->Insert( q2, phi );			// insert into histogram
	  q2 += step;		
	}
	myHist->Output(string(fn));
  }
  p_G = myHist;
}


Complex Tau_Three_Pseudo::DFM::BW( double s, double M2, double MG )
{
  return M2/Complex(M2-s,-1.*MG );
}

Complex Tau_Three_Pseudo::DFM::BW_V( int a, int b, double s )
{
  double MassWidth  = sqrt(s)*widthV[a]*msV[a]/s
                    * pow( msV[a]/s*Sqrt_Lambda(s,ms[a],ms[b])/Sqrt_Lambda(msV[a],ms[a],ms[b]), 3 );
  return BW( s, msV[a],  MassWidth );
}

Complex Tau_Three_Pseudo::DFM::BW_v( int a, int b, double s )
{
  double MassWidth  = sqrt(s)*widthv[a]*msv[a]/s
                    * pow( msv[a]/s*Sqrt_Lambda(s,ms[a],ms[b])/Sqrt_Lambda(msv[a],ms[a],ms[b]), 3 );
  return BW( s, msv[a], MassWidth );
}

double Tau_Three_Pseudo::DFM::G( double s )
{
  double val;	
  p_G->Extrapolate( s+p_G->BinSize(), &val, 0 );
  							// shift +BinSize to get right value
  double p = val;
  return p;
}

double Tau_Three_Pseudo::DFM::IntegralG( double Q2 )
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
 
Complex Tau_Three_Pseudo::DFM::Trho( double s )
{
  double MV (sqrt(msV[0])), Mv (sqrt(msv[0])), GV (widthV[0]), Gv (widthv[0]);
  double MV2 = sqr(MV), Mv2 = sqr(Mv), m2 = ms[0], beta = Beta[0];
  double MG_rho = MV*GV ;
  double MG_rhop = Mv*Gv;
  MG_rho  = sqrt(s)*GV*MV2/s*pow((s-4.*m2)/(MV2-4.*m2),1.5);
  MG_rhop = sqrt(s)*Gv*Mv2/s*pow((s-4.*m2)/(Mv2-4.*m2),1.5);
  Complex BW_rho  = BW( s, MV2, MG_rho );
  Complex BW_rhop = BW( s, Mv2, MG_rhop );
  return ( BW_rho + beta*BW_rhop ) / (1.+beta);
}

Complex Tau_Three_Pseudo::DFM::BW_A( double s )
{
  double MAGA (1.);
  if (running_width & 1) {
	MAGA = sqrt(s)*GA*G(s)/G(MA2);
  }
  else MAGA = sqrt(MA2)*GA;
  return BW( s, MA2, MAGA );
}


Complex Tau_Three_Pseudo::DFM::Tvector1( int a, int b, double x )
{
  Complex ret =			  
              BW_V(a,b,x)*( 1.-(ms[a]-ms[b])/(3.*msV[a]) )
  + Beta[a]*( BW_v(a,b,x)*( 1.-(ms[a]-ms[b])/(3.*msv[a]) ) );
  return ret/(1.+Beta[a]);
}

Complex Tau_Three_Pseudo::DFM::Tvector2( int a, int b, double x )
{
  Complex ret =
  			  BW_V(a,b,x)/msV[a] 
    + Beta[a]*BW_v(a,b,x)/msv[a];
  return ret*( 2.*(ms[a]-ms[b]) )/( 3.*(1.+Beta[a]) );	
}

Complex Tau_Three_Pseudo::DFM::TSvector( int a, int b, int c, double Q2, double s, double t )
{
  Complex ret =
              BW_V(a,b,s) * (   Q2-2.*t-s+2.*ms[a]+ms[c] 
		                      - (ms[a]-ms[b])/msV[a]*( Q2+t-ms[c]-Q2*(t-msV[a]) ) )
	+ Beta[a]*BW_v(a,b,s) * (   Q2-2.*t-s+2.*ms[a]+ms[c] 
		                      - (ms[a]-ms[b])/msv[a]*( Q2+t-ms[c]-Q2*(t-msv[a]) ) );
  return ret/(1.+Beta[a]);	
}

Complex Tau_Three_Pseudo::DFM::Tgen( int a, int b, int c, double s, double t)
{
  return( Tvector1(a-1,c-1,s) + Tvector2(b-1,c-1,t) );
}


Complex Tau_Three_Pseudo::DFM::FormFactor( int j, double Q2, double s, double t )
{
  double X123 = 2*ms[0];
  double ms123 = ms[0];
  Complex FF(0.,0.);
  switch( j ) {
	case 1 : FF = BW_A(Q2)*Tgen(1,2,3,s,t);
			 break;
	case 2 : FF = BW_A(Q2)*Tgen(2,1,3,t,s); 
			 break;
	case 3 : FF = X123 + BW_A(Q2)*(Q2-MA2)*ms123/(MA2*Q2) 
	                     *( TSvector(1-1,3-1,2-1,Q2,s,t) + TSvector(2-1,3-1,1-1,Q2,t,s) );
			 FF /= 2.*(Q2-ms123);			 
			 break;
  }
  return FF;
}
 
// General framework

Complex Tau_Three_Pseudo::FormFactor( int j, double Q2, double s, double t )
{
  switch( m_md.ff ) {
	case 1 : return dfm.FormFactor(j,Q2,s,t);		
	case 2 : return dpp.FormFactor(j,Q2,s,t);	
  }
}

double Tau_Three_Pseudo::Using_Hels( const Vec4D *_p )
{
  XYZFunc F(m_nout,_p,p_flavs);
  double ret = 0.;
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
  Complex ampl(0.,0.);
  for( int h=0; h<4; h++ ) {
	ampl +=  F.X( m_nutau, m_pseudo_1, 0, h, cR, cL ) * ( FS + F1*(1.-d1) - F2*d2 )
	       + F.X( m_nutau, m_pseudo_2, 0, h, cR, cL ) * ( FS - F1*d1 + F2*(1.-d2) )
	       + F.X( m_nutau, m_pseudo_3, 0, h, cR, cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) ); 
  }
  F.Delete();
  return norm( ampl )/4.;
} // ME is invariant under exchange of two identical pseudos !


double Tau_Three_Pseudo::operator()( const Vec4D *_p )
{
  double T = Using_Hels( _p );
  return T*sqr(2.*GF*m_A123/fpi/3.); 
}
