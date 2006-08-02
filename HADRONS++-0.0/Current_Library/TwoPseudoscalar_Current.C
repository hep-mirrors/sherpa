#include "TwoPseudoscalar_Current.H"
#include "Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void TwoPseudoscalar_Current::SetModelParameters( struct GeneralModel _md )
{
  m_pionmode = (p_flavs[1].Kfcode() == kf::pi_plus) ? 1 : 0;
  m_ff       = int( _md("FORM_FACTOR", 1 ) );
  m_fpi      = _md("fpi", 0.0924 );
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real() );
  double CG  = m_pionmode ? 1. : SQRT_05;   // Clebsch-Gordon
  m_global   = CG * Vud / SQRT_05;           // GF * V_CKM * CG

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
Complex TwoPseudoscalar_Current::A( double x, double y )
{
  Complex ret(0.,0.); 
  Complex sigma = csqrt(1.-4.*x);
  ret = log(y) + 8.*x - 5./3. + pow(sigma,3.)*log((sigma+1.)/(sigma-1.)); 
  return ret;
}


Complex TwoPseudoscalar_Current::FormFactor( double s )
{
  Complex ret(1.,0.);
  if( m_ff == 1 ) {         // Breit-Wigner-rho
//     std::cout<<"s="<<s<<std::endl;
    Complex BWr   = m_R.BreitWigner(s);
    Complex BWrr  = m_RR.BreitWigner(s);
    Complex BWrrr = m_RRR.BreitWigner(s);
    ret = ( BWr + m_beta*BWrr + m_gamma*BWrrr )/( 1.+m_beta+m_gamma );
//     std::cout<<"BWr="<<BWr<<std::endl;
//     std::cout<<"BWrr="<<BWrr<<std::endl;
//     std::cout<<"BWrrr="<<BWrrr<<std::endl;
//     std::cout<<"m_beta="<<m_beta<<std::endl;
//     std::cout<<"m_gamma="<<m_gamma<<std::endl;
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


void TwoPseudoscalar_Current::Calc()
{
//   std::cout<<"mom[0]="<<p_moms[0]<<std::endl;
//   std::cout<<"mom[1]="<<p_moms[1]<<std::endl;
  double  q2 = (p_moms[1] + p_moms[0] ).Abs2();
  Complex FF = FormFactor(q2);
//   std::cout<<om::red<<"FF="<<FF<<om::reset<<std::endl;
//   std::cout<<om::red<<"m_global="<<m_global<<om::reset<<std::endl;
  p_results[0] = m_global*FF*ComplexVec4D( (p_moms[1]-p_moms[0]), Vec4D(0.0,0.0,0.0,0.0) );
//   std::cout<<om::red<<"current="<<p_currents[0]<<om::reset<<std::endl;
}



DECLARE_GETTER(TwoPseudoscalar_Current_Getter, "TwoPseudoscalar_Current",
               Current_Base,Flavour_Info);

Current_Base* TwoPseudoscalar_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new TwoPseudoscalar_Current(parameters.flavs, parameters.nout, parameters.indices, "TwoPseudoscalar_Current");
}

void TwoPseudoscalar_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}
