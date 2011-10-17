#include "HADRONS++/Current_Library/VA_0_PP_strange.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void VA_0_PP_strange::SetModelParameters( struct GeneralModel _md )
{
  m_chpionmode = (m_flavs[p_i[0]].Kfcode() == kf_pi_plus) ? 1 : 0;
  
  for( int i=0; i<2; i++ ) {
    m_ms[i] = sqr(m_flavs[p_i[i]].HadMass());
  }
  double Vus  =_md("Vus", Tools::Vus);
  m_global   = Vus;
  if( !m_chpionmode ) m_global *= SQRT_05;
  m_Delta_KP = m_ms[1] - m_ms[0];
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : p_ff = new RChT(_md);
    break;               // use RChT on own risk: not tested sufficiently
    case 1 : p_ff = new KS(_md);
    break;
  }
  p_ff->SetMasses2( m_ms[0], m_ms[1], sqr(Flavour(kf_eta).HadMass()) );
  
}


void VA_0_PP_strange::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  // 0 is pion, 1 kaon
  double  q2 = (moms[p_i[0]]+moms[p_i[1]]).Abs2();
  Complex FS = p_ff->ScalarFormFactor(q2);
  Complex FV = p_ff->VectorFormFactor(q2);
  Complex termK = m_Delta_KP/q2*(FS-FV)+FV;
  Complex termP = m_Delta_KP/q2*(FS-FV)-FV;
  
  Insert((m_global*termK)*moms[p_i[1]]+(m_global*termP)*moms[p_i[0]],0);
}


VA_0_PP_strange::FF_Base::FF_Base( GeneralModel _md )
{
  int    running = int( _md("RUNNING_WIDTH", 1 ) );
  double MR      = _md("Mass_K*(892)+", 0.8921 );
  double GR      = _md("Width_K*(892)+", 0.0513 );
  double MR0     = _md("Mass_K(0)*(1430)+", 1.396 );
  double GR0     = _md("Width_K(0)*(1430)+", 0.294 );
  m_R       = ResonanceFlavour(kf_K_star_892_plus, MR, GR, running);
  m_R0      = ResonanceFlavour(kf_K_0_star_1430_plus, MR0, GR0, running );
  m_fpi2     = sqr( _md("fpi", 0.1307) );
}
 
void VA_0_PP_strange::FF_Base::SetMasses2( double _mPi2, double _mK2, double _mEta2 )
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

VA_0_PP_strange::RChT::RChT( GeneralModel _md )
  : FF_Base(_md)
{
  m_fpi2	 /= 2.;			// redefinition of fpi
  m_MK2      = m_R.Mass2();
  m_GK       = m_R.Width();
  m_MK02     = m_R0.Mass2();
  m_GK0      = m_R0.Width();
  m_renorm2  = sqr( _md("renorm",_md("Mass_rho(770)+",Flavour(kf_rho_770_plus).HadMass())));
  m_cd       = _md("const_cd", 0.014);
  m_cm       = m_fpi2/4./m_cd;
}

double VA_0_PP_strange::RChT::MassWidthVector( double s )
{
  double ret (0.);
  if (s>sqr(m_mK+m_mPi))  ret += pow( Tools::Lambda(s,m_mK2,m_mPi2), 1.5 );
  if (s>sqr(m_mK+m_mEta)) ret += pow( Tools::Lambda(s,m_mK2,m_mEta2), 1.5 );
  ret *= m_MK2 /( 128.*M_PI*m_fpi2*sqr(s) );
  return ret;
}

double VA_0_PP_strange::RChT::MassWidthScalar( double s )
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

Complex VA_0_PP_strange::RChT::JBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double nu = sqrt( Tools::Lambda(s,MP2,MQ2) );
  double  J = 2. + Delta/s*log(MQ2/MP2) 
        - Sigma/Delta*log(MQ2/MP2) 
        - nu/s*log( (sqr(s+nu)-sqr(Delta))/(sqr(s-nu)-sqr(Delta)) ); 
  return J/(32.*sqr(M_PI));
}

Complex VA_0_PP_strange::RChT::JBarBar( double s, double MP2, double MQ2, double Sigma, double Delta )
{
  double Jp_0 = ( Sigma/sqr(Delta) + 2.*MP2*MQ2/pow(Delta,3)*log(MQ2/MP2) )/( 32.*sqr(M_PI) );
  return JBar(s,MP2,MQ2,Sigma,Delta) - s*Jp_0;
}

Complex VA_0_PP_strange::RChT::Mr( double s, double MP2, double MQ2 )
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

Complex VA_0_PP_strange::RChT::L( double s, double MP2, double MQ2 )
{
  double Sigma = MP2 + MQ2;
  double Delta = MP2 - MQ2;
  Complex Jb   = JBar(s,MP2,MQ2,Sigma,Delta);
  return( sqr(Delta)/(4.*s)*Jb );
}

double VA_0_PP_strange::RChT::MuOf( double m2 )
{
  double mu2 = m_renorm2;
  return m2/(32.*sqr(M_PI)*m_fpi2) * log(m2/mu2);
}

Complex VA_0_PP_strange::RChT::VectorFormFactor( double s )
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

Complex VA_0_PP_strange::RChT::ScalarFormFactor( double s )
{
  Complex ret(1.,0.);
  double  MG_K0 = MassWidthScalar(s);
  Complex    BW = Tools::BreitWigner( s, m_MK02, MG_K0 );
  Complex    F4 = 1./(8.*m_fpi2)
      * ( 5.*s - 2.*m_Sigma_KP - 3.*sqr(m_Delta_KP)/s )
      * JBar(s,m_mK2,m_mPi2,m_mK2+m_mPi2,m_mK2-m_mPi2)
      + 1./(24.*m_fpi2)
      * ( 3.*s - 2.*m_Sigma_KP - sqr(m_Delta_KP)/s )
      * JBar(s,m_mK2,m_mEta2,m_mK2+m_mEta2,m_mK2-m_mEta2)
	  + s/(4.*m_Delta_KP)*(5.*MuOf(m_mPi2)-2.*MuOf(m_mK2)-3.*MuOf(m_mEta2));
  double  inter = 1. - (1.-m_fpi2/4./sqr(m_cd))*m_Sigma_KP/m_MK02;
  Complex  expon = Complex( F4.real(), F4.imag()/(1.+sqr(F4.imag())) );
  ret = BW * inter * exp(expon);
  return ret;
}

// Kuehn Santamaria Model

VA_0_PP_strange::KS::KS( GeneralModel _md )
  : FF_Base(_md)
{
  double MRR 		= _md("Mass_K*(1680)+", 1.700 );
  double GRR 		= _md("Width_K*(1680)+", 0.235 );
  m_RR       		= ResonanceFlavour(kf_K_star_1410_plus, MRR, GRR, int(_md("RUNNING_WIDTH",1)));
  m_beta 	 		= _md("beta", -0.038);
}
 

Complex VA_0_PP_strange::KS::VectorFormFactor( double s )
{
  return (m_R.BreitWigner(s)+m_beta*m_RR.BreitWigner(s))/(1.+m_beta);
}

Complex VA_0_PP_strange::KS::ScalarFormFactor( double s )
{
  return m_R0.BreitWigner(s);
}

DEFINE_CURRENT_GETTER(VA_0_PP_strange,VA_0_PP_strange_Getter,"VA_0_PP_strange")

void VA_0_PP_strange_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi K $ \n\n"
    <<"Order: 0 = $\\pi$, 1 = K \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: http://sherpa-mc.de/dokuwiki/\\_media/publications/theses/diplom\\_laubrich.pdf \n"
    <<std::endl;
}
