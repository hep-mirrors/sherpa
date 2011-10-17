#include "HADRONS++/Current_Library/VA_0_PP.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void VA_0_PP::SetModelParameters( struct GeneralModel _md )
{
  m_pionmode = (m_flavs[p_i[1]].Kfcode() == kf_pi_plus) ? 1 : 0;
  m_ff       = int( _md("FORM_FACTOR", 1 ) );
  m_fpi      = _md("fpi", 0.1307 );
  double Vud = _md("Vud", Tools::Vud);
  double CG  = m_pionmode ? 1. : SQRT_05 ;   // weight factor
  m_global   = CG * Vud / SQRT_05;           // GF * V_CKM * CG

  int    running  = int( _md("RUNNING_WIDTH", 1 ) );
  double MR       = _md("Mass_rho(770)+", 0.7769 );
  double MRR      = _md("Mass_rho(1450)+", 1.363 );
  double MRRR     = _md("Mass_rho(1700)+", 1.700 );
  double GR       = _md("Width_rho(770)+", 0.149 );
  double GRR      = _md("Width_rho(1450)+", 0.310 );
  double GRRR     = _md("Width_rho(1700)+", 0.235 );
  m_R   = ResonanceFlavour( kf_rho_770_plus, MR, GR, running );
  m_RR  = ResonanceFlavour( kf_rho_1450_plus, MRR, GRR, running );
  m_RRR = ResonanceFlavour( kf_rho_1700_plus, MRRR, GRRR, running );

  m_m2_pi    = sqr( Flavour(kf_pi_plus).HadMass() );
  m_m2_K     = sqr( Flavour(kf_K_plus).HadMass() );
   
  // coefficients for KS model
  m_beta     = _md("beta", -0.167 );
  m_gamma    = _md("gamma", 0.05 );

  // coefficients for RChT model
  m_gammaR   = _md("gamma_rho_770", 1. );
  m_gammaRR  = _md("gamma_rho_1450", 1. );
  m_gammaRRR = _md("gamma_rho_1700", 1. );

  // redefine fpi for RChT
  if( m_ff == 2 ) m_fpi *= SQRT_05;
}


// loop function
Complex VA_0_PP::A( double x, double y )
{
  Complex ret(0.,0.); 
  Complex sigma = csqrt(1.-4.*x);
  ret = log(y) + 8.*x - 5./3. + pow(sigma,3.)*log((sigma+1.)/(sigma-1.)); 
  return ret;
}


Complex VA_0_PP::FormFactor( double s )
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
    double MG_R;
    Complex AA = A( m_m2_pi/s, m_m2_pi/m_R.Mass2() ) + 0.5*A( m_m2_K/s, m_m2_K/m_R.Mass2() );
    double expon = -1.*s/(96.*sqr(M_PI*m_fpi))*AA.real();
    Complex BW_1;
    if (m_R.Running()) {
      MG_R   = -m_gammaR  *1.*m_R.Mass2()  *s/(96.*sqr(M_PI*m_fpi)) * AA.imag();
      BW_1 = Tools::BreitWigner( s, m_R.Mass2(), MG_R );
    }
    else {
      MG_R   = m_R.MassWidth();
      BW_1 = Tools::BreitWignerFix( s, m_R.Mass2(), MG_R );
    }
    ret = BW_1 * exp(expon);
    return ret;
  }
  return Complex(1.,0.);
}


void VA_0_PP::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  double  q2 = (moms[p_i[1]] + moms[p_i[0]] ).Abs2();
  Complex FF = FormFactor(q2);
  Insert(m_global*FF*(moms[p_i[1]]-moms[p_i[0]]), 0);
}

DEFINE_CURRENT_GETTER(VA_0_PP,VA_0_PP_Getter,"VA_0_PP")

void VA_0_PP_Getter::PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi $ \n\n"
    <<"Order: 0 = $\\pi^0$, 1 = $\\pi^\\pm$ \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: http://sherpa-mc.de/dokuwiki/\\_media/publications/theses/diplom\\_laubrich.pdf \n"
    <<std::endl;
}
