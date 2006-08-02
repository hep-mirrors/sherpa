#include "FourPionThreeCharged_Current.H"
#include "Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

FourPionThreeCharged_Current::FourPionThreeCharged_Current(const ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
    Current_Base(flavs, n, indices, name) {}


void FourPionThreeCharged_Current::SetModelParameters( struct GeneralModel _md )
{
  double Vud = _md("Vud", rpa.gen.ComplexMatrixElement(string("CKM"), 0, 0).real());
  double GF  = _md("GF", rpa.gen.ScalarConstant(string("GF")));
  m_global   = Vud;

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


void FourPionThreeCharged_Current::LorenzBase::SetPrivates( ATOOLS::Vec4D * _p ) 
{
//   p_X  = _x;
  p_p  = _p;
  m_q = p_p[0]+p_p[1]+p_p[2]+p_p[3];         // = q
  m_q2 = m_q.Abs2();
  for (int i=2; i<=4; i++ ) {
    m_r[i] = m_q-p_p[i-1];     // this is labelled m_r[i] instead of m_r[i-1] only for legacy reasons
    m_s[i] = (p_p[0]+p_p[i-1]).Abs2();
  }
  // redefinition of variables
  m_s[0] = (p_p[1]+p_p[2]).Abs2();          // = t3
  m_s[1] = (p_p[1]+p_p[3]).Abs2();          // = t4
//   p_X[0] = p_X[1]+p_X[2]+p_X[3]+p_X[4];     // X(nu,q,0)
  // unused variables yet
  m_r[0] = m_r[1] = ATOOLS::Vec4D(0.,0.,0.,0.);
};
 
// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

FourPionThreeCharged_Current::KS::KS( string path, GeneralModel _md )
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

Complex FourPionThreeCharged_Current::KS::Fk( double x, Complex * _beta )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (_beta[0] + _beta[1]*BW_R + _beta[2]*BW_RR + _beta[3]*BW_RRR)/(_beta[0]+_beta[1]+_beta[2]+_beta[3]) );
}

Complex FourPionThreeCharged_Current::KS::Trho( double x )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Complex FourPionThreeCharged_Current::KS::TTrho( double x )
{
  Complex BW_R   = m_Rho.BreitWignerAlt(x);
  Complex BW_RR  = m_RR.BreitWignerAlt(x);
  return( BW_R + m_sigma*BW_RR );
}

double FourPionThreeCharged_Current::KS::Dots( int k, int l )        // k=3,4; l=1,2,3,4 (l!=k)
  // special dot products used for anomalous part of J_omega
{
  int pre = l-1;                // predecessor of l
  if (pre==0) pre  = 4;
  if (pre==k) pre -= 1;     
  int suc  = l+1;               // successor of l
  if (suc==k)  suc += 1;
  if (suc==5)  suc  = 1;
  return( (m_r[k]*p_p[pre-1])*(p_p[suc-1]*p_p[k-1]) 
      - (m_r[k]*p_p[suc-1])*(p_p[pre-1]*p_p[k-1]) );
}

ComplexVec4D FourPionThreeCharged_Current::KS::OmegaPi()
{
  // chiral part
  ComplexVec4D J_chi( 0.0,0.0,0.0,0.0 );
  ComplexVec4D sum1( 0.0,0.0,0.0,0.0 );
  ComplexVec4D sum2( 0.0,0.0,0.0,0.0 );
  for (int k=2; k<=4; k++) {
    sum2 = ComplexVec4D( 0.0,0.0,0.0,0.0 );
    for (int l=2; l<=4; l++) if (l!=k) {
      sum2 += ComplexVec4D(m_q-2.*p_p[l-1]) * ((m_r[l]*(p_p[k-1]-p_p[0]))/m_r[l].Abs2());
    }
    sum1 += m_R[k]*Trho(m_s[k]) * (ComplexVec4D(p_p[k-1]-p_p[0])-sum2);
  }
  J_chi = 2.*sqrt(3.)/m_fpi2*Trho(m_q2) * sum1;

  // anomalous part
  ComplexVec4D J_a( 0.0,0.0,0.0,0.0 );
  sum1=ComplexVec4D(Vec4D(0.,0.,0.,0.));
  for( int k=3; k<=4; k++ ) {
    sum2=ComplexVec4D(Vec4D(0.,0.,0.,0.));
    for (int l=1; l<=4; l++) if (l!=k) {
      sum2 += ComplexVec4D(p_p[l-1]) * Dots(k,l);
    }
    sum1 += m_O.BreitWignerAlt(m_r[k].Abs2())*sum2;
  }
  J_a = m_Go3p*m_Frho*m_grop*TTrho(m_q2) * sum1;

  // total
  return( (J_chi + J_a)*Fk(m_q2,m_Beta_opi) );
}

ComplexVec4D FourPionThreeCharged_Current::KS::AonePi()
{
  ComplexVec4D term1( 0.0,0.0,0.0,0.0 ), term2( 0.0,0.0,0.0,0.0 );
  Complex A, B, C;
  Vec4D   P, Q, R;

  //  1st term
  P = p_p[0]-p_p[2];
  Q = p_p[0]-p_p[3];
  R = m_r[2];
  A = m_Rho.BreitWigner(m_s[3]);
  B = m_Rho.BreitWigner(m_s[4]);
  C = m_A.BreitWigner(R.Abs2());
  term1  = A*ComplexVec4D(p_p[0]-p_p[2]) 
         + B*ComplexVec4D(p_p[0]-p_p[3]) 
         - ComplexVec4D(m_q-p_p[1]) * (( A*(R*P) + B*(R*Q) )/R.Abs2())
         - ComplexVec4D(m_q) * (( A*(m_q*P) 
                               + B*(m_q*Q) 
                               - (m_q*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() 
                               )/m_q2);
  term1 *= C;

  // 2nd and 3rd term
  int ind;
  ComplexVec4D help( 0.0,0.0,0.0,0.0 );
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    P = p_p[0]-p_p[1];
    Q = p_p[ind-1]-p_p[1];
    R = m_r[k];
    A = m_Rho.BreitWigner(m_s[2]);
    B = m_Rho.BreitWigner(m_s[ind-3]);
    C = m_A.BreitWigner(R.Abs2());
    help   = A*ComplexVec4D(p_p[0]-p_p[1]) 
        + B*ComplexVec4D(p_p[1]-p_p[ind-1]) 
        - ComplexVec4D(m_q-p_p[k-1]) * (( A*(R*P) + B*(R*Q) )/R.Abs2())
        - ComplexVec4D(m_q)*(( A*(m_q*P) 
                            + B*(m_q*Q) 
                            - (m_q*Q)*(A*(R*P)+B*(R*Q))/R.Abs2() 
                             )/m_q2);
    term2 += C*help;
  }

  // total 
  return (term1-term2)*Fk(m_q2,m_Beta_api);
}

ComplexVec4D FourPionThreeCharged_Current::KS::SigmaRho()
{
  int ind;
  ComplexVec4D term( 0.0,0.0,0.0,0.0 );
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    BW_S = m_S.BreitWigner(m_s[k]);
    BW_R = m_Rho.BreitWigner(m_s[ind-3]);
    term += BW_S*BW_R * ComplexVec4D( p_p[1] - p_p[ind-1] + (m_q*(p_p[ind-1]-p_p[1]))*m_q );
  }
  return term*Fk(m_q2,m_Beta_srh);
}

ComplexVec4D FourPionThreeCharged_Current::KS::FzeroRho()
{
  int ind;
  ComplexVec4D term( 0.0,0.0,0.0,0.0 );
  Complex BW_S, BW_R;
  for (int k=3; k<=4; k++) {
    ind = (k==3)? 4 : 3;
    BW_S = m_F.BreitWigner(m_s[k]);
    BW_R = m_Rho.BreitWigner(m_s[ind-3]);
    term += BW_S*BW_R * ComplexVec4D( p_p[1] - p_p[ind-1] + (m_q*(p_p[ind-1]-p_p[1]))*m_q );
  }
  return term*Fk(m_q2,m_Beta_frh);
}

ComplexVec4D FourPionThreeCharged_Current::KS::operator()( int number )
{
  PRINT_INFO("number="<<number);
  switch( number ) {
    case 0 : return OmegaPi();
    case 1 : return AonePi();
    case 2 : return SigmaRho();
    case 3 : return FzeroRho();
  }
}

// General framework

void FourPionThreeCharged_Current::Calc()
{
//   XYZFunc F(m_nout,_p,p_flavs,k0_n);
  // internal numeration and convenient variables
//   for (int i=1; i<=5; i++ ) m_p[i] = _p[m_inter[i]];
  // create amplitudes tensor
//   _ampls_tensor->clear();
  ComplexVec4D help( 0.0,0.0,0.0,0.0 );
//   for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
    // pre-calculate X-funcs
//     for (int k=1; k<=4; k++) 
//       m_X[k] = F.X( m_nutau, m_inter[k], 0, h, m_cR, m_cL );
//     m_X[0] = m_X[1] + m_X[2] + m_X[3] + m_X[4];
    // sum over all contributions
  p_lorenz->SetPrivates( p_moms );

  for (int k=0; k<m_ncontrib; k++) {
    help += (m_Alpha[k]/m_SumAlpha) * (*p_lorenz)(k) ;
  }
  p_results[0]=m_global*help;
//     _ampls_tensor->push_back( help*m_global );
//   }
//   F.Delete();
  // create index bookkeeping (using internal numbers 0 -> 1 2)
  // with pair (number, 2*spin); note: reversed order
//   _indices->clear();
//   _indices->push_back( pair<int,int>(0,1) );
//   _indices->push_back( pair<int,int>(m_nutau,1) );
  // note: pions do not have spin index
}



DECLARE_GETTER(FourPionThreeCharged_Current_Getter, "FourPionThreeCharged_Current",
               Current_Base,Flavour_Info);

Current_Base* FourPionThreeCharged_Current_Getter::operator()(const Flavour_Info &parameters) const
{
  return new FourPionThreeCharged_Current(parameters.flavs, parameters.nout, parameters.indices, "FourPionThreeCharged_Current");
}

void FourPionThreeCharged_Current_Getter::
    PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"implement me";
}
