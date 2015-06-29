#include "HADRONS++/Current_Library/VA_0_PPP.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

VA_0_PPP::FF_Base::~FF_Base()
{
}

VA_0_PPP::VA_0_PPP(const ATOOLS::Flavour_Vector& flavs,
                   const std::vector<int>& indices, const std::string& name) :
  Current_Base(flavs, indices, name)
{
  int nPion_0(0), nPion_ch(0), nKaon_0(0), nKaon_ch(0), nKaon_S (0), nKaon_L(0);
  // count number of pions, kaons and calc. mass^2
  for( int i=0; i<3; i++ ) {
    if( m_flavs[p_i[i]].Kfcode() == kf_pi_plus )   nPion_ch++;
    if( m_flavs[p_i[i]].Kfcode() == kf_pi )        nPion_0++;
    if( m_flavs[p_i[i]].Kfcode() == kf_K_plus )    nKaon_ch++;
    if( m_flavs[p_i[i]].Kfcode() == kf_K )         nKaon_0++;
    if( m_flavs[p_i[i]].Kfcode() == kf_K_L )	   nKaon_L++;
    if( m_flavs[p_i[i]].Kfcode() == kf_K_S )       nKaon_S++;
    m_ms[i] = sqr( m_flavs[p_i[i]].HadMass() );
  }
  // sanity check
  if (nPion_ch+nPion_0+nKaon_ch+nKaon_L+nKaon_S+nKaon_0 != 3) {
    msg_Error()<<"ERROR in HADRONS::VA_0_PPP constructor\n"
        <<"     number of three outgoing pseudoscalars != 3 ?!.\n"
        <<"     Don't know, what to do. Will abort."<<endl;
    abort();           
  }
  // define mode number
  m_mode = nPion_ch*1000 + nPion_0*100 + nKaon_ch*10 + (nKaon_L+nKaon_S+nKaon_0);
  m_kaon_mode = 100*nKaon_L + 10*nKaon_S + nKaon_0;
}

void VA_0_PPP::SetModelParameters( struct GeneralModel _md )
{
  m_Vud      = _md("Vud", Tools::Vud);
  m_Vus      = _md("Vus", Tools::Vus);
  double fpi = _md("fpi", 0.1307 );
  
  m_B123 = 0.;
  double A123;
  switch( m_mode ) {
    case 3000 : /* pi- pi- pi+ mode */
                A123 = m_Vud*SQRT_05;
				m_B123 = 0.;
                break;
    case 1200 : /* pi0 pi0 pi- mode */
                A123 = m_Vud*SQRT_05;
				m_B123 = 0.;
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
				m_B123 = -2./3.;
                break;
    case  210 : /* pi0 pi0 K- mode */
                A123 = m_Vus/4.;
				m_B123 = 1.;
                break;
    case 2010 : /* K- pi- pi+ */
                m_B123 = -1.;
                A123 = -0.5*m_Vus;
                break;
    case 1101 : /* pi- K0b pi0 */
                m_B123 = -2/3.;
                A123 = 1.5*m_Vus*SQRT_05;
                break;
    case   30 : /* K- K- K+ */
                A123 = m_Vus;
                break;
    case   12 : /* K- K0b K0 */
                A123 = -0.5*m_Vus;
                break;
    default   : msg_Error()<<"Warning in HADRONS::Tau_Decay_MEs.C in Tau_Three_Pseudo::GetA123() :\n"
                           <<"     Obviously this three pseudoscalar channel (code "<<m_mode<<")\n"
                           <<"     doesn't have a global A123. Maybe it is not implemented yet.\n"
                           <<"     Take A123=1., will continue and hope for the best."<<endl;
                A123 = 1.;
                break;
  }
  m_global   = 4./(3.*fpi)*A123;
  switch( int(_md("FORM_FACTOR", 1)) ) {
    case 2 : p_ff = new RChT(m_mode,m_kaon_mode, m_path,_md,m_ms);
    break;
    case 1 : p_ff = new KS95(m_mode,m_kaon_mode, m_path,_md,m_ms);
    break;
    case 3 : p_ff = new KS(m_mode,m_kaon_mode, m_path,_md,m_ms);
    break;
  }
}


void VA_0_PPP::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4D p1( moms[p_i[0]] ),
        p2( moms[p_i[1]] ),
        p3( moms[p_i[2]] );
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
//  F1 = Complex(0.,0.);
//  F2 = Complex(0.,0.);
//  FS = Complex(0.,0.);
  
  Insert( m_global*( FS+F1*(1.-d1)-F2*d2 )*Vec4C(p1)
	+ m_global*( FS-F1*d1+F2*(1.-d2) )*Vec4C(p2)
	+ m_global*( FS-F1*(1.+d1)-F2*(1.+d2) )*Vec4C(p3)
	+ m_global*Complex(0.,1.)*FV*Vec4C(cross(p1,p2,p3)) , 0);
//   for( int h=0; h<4; h++ ) {        // helicity comb. (nutau,tau)
//     _ampls_tensor->push_back( 
//         ( F.X( m_nutau, m_pseudo_1, 0, h, m_cR, m_cL ) * ( FS + F1*(1.-d1) - F2*d2 )
//         + F.X( m_nutau, m_pseudo_2, 0, h, m_cR, m_cL ) * ( FS - F1*d1 + F2*(1.-d2) )
//         + F.X( m_nutau, m_pseudo_3, 0, h, m_cR, m_cL ) * ( FS - F1*(1.+d1) - F2*(1.+d2) ) 
//         + Complex(0.,1.)*F.X( m_nutau, cross(p1,p2,p3), 0, h, m_cR, m_cL ) * FV
//         ) * m_global
//                             );
//   }
}


Complex VA_0_PPP::FormFactor( int j, double Q2, double s, double t )
{
  return( p_ff->FormFactor(j,Q2,s,t) );     
}


VA_0_PPP::FF_Base::FF_Base(int mode, int kaon_mode, std::string path, GeneralModel _md, double * _ms)
  : m_mode (mode), m_kaon_mode (kaon_mode)
{
  // set resonances and parameters
  m_ms[0] = _ms[0];
  m_ms[1] = _ms[1];
  m_ms[2] = _ms[2];
  m_deltas = 0;
  kf_code resA, resAA, resV[2], resVV[2];
  switch( m_mode ) {
    case 3000 : /* pi- pi- pi+ mode */
      resV[0] = kf_rho_770;
      resV[1] = kf_rho_770;
      break;
    case 1200 : /* pi0 pi0 pi- mode */
      resV[0] = kf_rho_770_plus;
      resV[1] = kf_rho_770_plus;
      break;
    case 1020 : /* K- pi- K+ */
      resV[0] = kf_rho_770;
      resV[1] = kf_K_star_892; 
      break;
    case 1002 : /* K0 pi- K0b */
      resV[0] = kf_rho_770;
      resV[1] = kf_K_star_892_plus; 
      break;
    case  111 : /* K- pi0 K0 */
      resV[0] = kf_rho_770_plus;
      resV[1] = kf_K_star_892;
      break;
    case  210 : /* pi0 pi0 K- mode */
      m_deltas= 1;
      resV[0] = kf_K_star_892_plus;
      resV[1] = kf_K_star_892_plus;
      break;
    case 2010 : /* K- pi- pi+ */
      m_deltas= 1;
      resV[0] = kf_K_star_892; 
      resV[1] = kf_rho_770;
      break;
    case 1101 : /* pi- K0b pi0 */
      m_deltas= 1;
      resV[0] = kf_rho_770_plus;
      resV[1] = kf_K_star_892;
      break;
    case   30 : /* K- K- K+ */
      m_deltas= 1;
      resV[0] = kf_rho_770;
      resV[1] = kf_rho_770;
      break;
    case   12 : /* K- K0b K0 */
      m_deltas= 1;
      resV[0] = kf_rho_770_plus;
      resV[1] = kf_rho_770_plus;
      break;
    default:
      msg_Error()<<METHOD<<" Warning: Three pseudoscalar mode not recognized. "
        <<"mode="<<m_mode<<endl;
      break;
  }

  // higher resonances
  for( int i=0; i<2; ++i ) {
    if( resV[i] == kf_rho_770 )          resVV[i] = kf_rho_1450;
    if( resV[i] == kf_rho_770_plus )     resVV[i] = kf_rho_1450_plus;
    if( resV[i] == kf_K_star_892 )       resVV[i] = kf_K_star_1410;
    if( resV[i] == kf_K_star_892_plus )  resVV[i] = kf_K_star_1410_plus;
  }

  // axial resonance
  int running    = int( _md("RUNNING_WIDTH", 3 ) );     // running width
  double mA, mAA, wA, wAA;
  if( m_deltas ) {
	running = 0;
	resA  = kf_K_1_1400_plus;
    resAA = kf_K_1_1270_plus;
	mA      = _md("Mass_"+Flavour(resA).IDName(),     1.402 );          // mass of axial resonance
	mAA     = _md("Mass_"+Flavour(resAA).IDName(),    1.270 );         // mass of axial resonance'
	wA      = _md("Width_"+Flavour(resA).IDName(),     0.174 );          // width of axial resonance
	wAA     = _md("Width_"+Flavour(resAA).IDName(),    0.090 );         // width of axial resonance'
	m_alpha        = _md("alpha_"+Flavour(resAA).IDName(),  1.3 );                              // weight factor for A'
  }
  else {
	resA  = kf_a_1_1260_plus;
	resAA = kf_a_1_1260_plus;
	mA      = _md("Mass_"+Flavour(resA).IDName(),     1.251 );          // mass of axial resonance
	mAA     = _md("Mass_"+Flavour(resAA).IDName(),    1.251 );         // mass of axial resonance'
	wA      = _md("Width_"+Flavour(resA).IDName(),     0.599 );          // width of axial resonance
	wAA     = _md("Width_"+Flavour(resAA).IDName(),    0.599 );         // width of axial resonance'
	m_alpha        = _md("alpha_"+Flavour(resAA).IDName(),  0. );                              // weight factor for A'
  }
  m_A = ResonanceFlavour( resA, mA, wA, running&1, path );
  m_AA = ResonanceFlavour( resAA, mAA, wAA, running&1, path );
  ResonanceFlavour help_rho = ResonanceFlavour(
      kf_rho_770_plus, 
      _md("Mass_rho(770)+", _md("Mass_rho(770)", 0.773 )),
      _md("Width_rho(770)+", _md("Width_rho(770)", 0.145 )),
      running&2 
      );
  ResonanceFlavour help_rhoP = ResonanceFlavour(
      kf_rho_1450_plus, 
      _md("Mass_rho(1450)+", _md("Mass_rho(1450)", 1.370 )),
      _md("Width_rho(1450)+", _md("Width_rho(1450)", 0.510 )),
      running&2 
      );
  double help_beta = _md("beta_rho(1450)+", _md("beta_rho(1450)", -0.145) );
  m_A.InitialiseThreeBodyResonance(help_rho, help_rhoP, help_beta);
  m_AA.InitialiseThreeBodyResonance(help_rho, help_rhoP, help_beta);

  double mV[2], mVV[2], wV[2], wVV[2];
  for( int i=0; i<2; ++i ) {
	if( resV[i] == kf_rho_770 || resV[i] == kf_rho_770_plus ) {
	  mV[i]       = _md("Mass_"+Flavour(resV[i]).IDName(),  0.773 );       // mass^2 of vector resonance ij
	  wV[i]       = _md("Width_"+Flavour(resV[i]).IDName(), 0.145 );        // width^2 of vector resonance ij
	  mVV[i]      = _md("Mass_"+Flavour(resVV[i]).IDName(), 1.370 );      // mass^2 of vector resonance' ij
	  wVV[i]      = _md("Width_"+Flavour(resVV[i]).IDName(),0.510 );       // width^2 of vector resonance' ij
	  m_Beta[i]   = _md("beta_"+Flavour(resVV[i]).IDName(),  -0.145 );                            // weight factor for Vij'
	}
	else if( resV[i] == kf_K_star_892 || resV[i] == kf_K_star_892_plus ) {
	  mV[i]       = _md("Mass_"+Flavour(resV[i]).IDName(),  0.8921 );       // mass^2 of vector resonance ij
	  wV[i]       = _md("Width_"+Flavour(resV[i]).IDName(), 0.0513 );        // width^2 of vector resonance ij
	  mVV[i]      = _md("Mass_"+Flavour(resVV[i]).IDName(), 1.412 );      // mass^2 of vector resonance' ij
	  wVV[i]      = _md("Width_"+Flavour(resVV[i]).IDName(),0.227 );       // width^2 of vector resonance' ij
	  m_Beta[i]   = _md("beta_"+Flavour(resVV[i]).IDName(),  -0.135 );                            // weight factor for Vij'
	}
	else {
	  mV[i]       = _md("Mass_"+Flavour(resV[i]).IDName(),  Flavour(resV[i]).HadMass());       // mass^2 of vector resonance ij
	  wV[i]       = _md("Width_"+Flavour(resV[i]).IDName(), Flavour(resV[i]).Width());        // width^2 of vector resonance ij
	  mVV[i]      = _md("Mass_"+Flavour(resVV[i]).IDName(), Flavour(resVV[i]).HadMass());      // mass^2 of vector resonance' ij
	  wVV[i]      = _md("Width_"+Flavour(resVV[i]).IDName(),Flavour(resVV[i]).Width());       // width^2 of vector resonance' ij
	  m_Beta[i]   = _md("beta_"+Flavour(resVV[i]).IDName(),  0. );                            // weight factor for Vij'
	}
	m_V[i]      = ResonanceFlavour( resV[i], mV[i], wV[i], running&2 );
	m_VV[i]     = ResonanceFlavour( resVV[i], mVV[i], wVV[i], running&2 );
  }

  m_fpi2      = sqr( _md("fpi", 0.1307) );          // pion decay constant
}


VA_0_PPP::RChT::RChT(int mode, int kaon_mode, string path, GeneralModel _md, double * _ms)
  : FF_Base(mode,kaon_mode,path,_md,_ms)
{
  // set resonances and parameters
  if( m_mode != 1200 && 
      m_mode != 3000 && 
      m_mode != 1020 ) {
    msg_Error()<<"Error: The mode "<<m_mode<<endl
        <<"     hasn't been implemented yet (RChT). Please use KS model."
        <<"     Don't know what to do. Will abort"<<endl;
    abort();
      }
   
      m_MO     = _md("Mass_omega(782)", Flavour(kf_omega_782).HadMass());
      m_MO2    = sqr(m_MO);
      m_GO     = _md("Width_omega(782)", Flavour(kf_omega_782).Width());
      
      m_l0     = _md("lambda0", 1.);                    // fit parameter lambda0
      m_gammaR = _md("gamma_rho(770)", 1.);              // global factor for rho width
      m_m      = Flavour( kf_pi_plus ).HadMass();         // pion mass
      m_m2     = sqr(m_m);                              // pion mass^2
      m_mK2    = sqr( Flavour( kf_K_plus ).HadMass() );   // Kaon mass^2
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

double VA_0_PPP::RChT::MassWidthVector( int a, double s )
{ 
  if(m_V[0].Running()) {
    switch(m_mode) {
      case 3000:
        case 1200: return MassWidthVector( s );
        case 1020: return (a==0)? MassWidthVector(s) : m_V[1].MassWidth();
    }
    // default:
    msg_Error()<<"Warning: this form factor (RChT) for the three-pseudoe mode "<<m_mode<<"\n"
        <<"     hasn't been implemented yet. Please use KS model."<<endl;
  }
  return( m_V[a].MassWidth() );
}

double VA_0_PPP::RChT::MassWidthVector( double s )
{ 
  double MVGV (0.);
  if( s>4.*m_m2 )  MVGV += pow( 1.-4.*m_m2/s, 1.5 );
  if( s>4.*m_mK2 ) MVGV += pow( 1.-4.*m_mK2/s, 1.5 ) / 2.;
  MVGV *= m_gammaR*m_V[0].Mass2()*s/(96.*M_PI*m_fpi2); 
  return MVGV;
}

double VA_0_PPP::RChT::MassWidthAxial( double Q2 )
{
  if( !m_deltas && m_A.Running() )
    return(  m_A.OffShellMassWidth(Q2)*pow(m_A.Mass2()/Q2,m_exp_alpha-2.) );
  return m_A.MassWidth();
}
 
Complex VA_0_PPP::RChT::FormFactor( int j, double Q2, double s, double t )
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
              Complex FAchi = Complex(1.,0.);
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
              Complex FSchi = 3./2.* Mass2(1)/(Q2-Mass2(1)) *
                  ( 1. + (Mass2(2)-u)/Q2 );
              Complex FS1r  = 3./2.*sqr(m_GV)/m_fpi2 * Mass2(1)/(Q2*(Q2-Mass2(1))) *
                  (   s*(t-u)/Complex(m_V[0].Mass2()-s,-1.*MassWidthVector(0,s))
                      + (t*(s-u)+(Q2-Mass2(0))*(Mass2(1)-Mass2(2)))/Complex(m_V[1].Mass2()-t,-1.*MassWidthVector(1,t)) 
                  );
              return FSchi + FS1r;
            }
            case 4: {                // vector contribution
              double c1c2c5 = 0.;
//              double c1c2c52c6 = -3./(32.*sqr(M_PI)) * m_V[0].Mass()/(sqrt(2.)*m_FV);
              double c1c2c52c6 = -3./(96.*sqr(M_PI)) * m_V[0].Mass()/(sqrt(2.)*m_GV);
              double c1c28c3c5 = 0.;
              double c4 = 0.;
//              double d3 = -3./(64.*sqr(M_PI)) * m_V[0].Mass2()/m_FV2 + m_fpi2/(8.*m_FV2);
              double d3 = -3./(192.*sqr(M_PI)) * m_V[0].Mass2()/(m_FV*m_GV);
              double d18d2d3 = m_fpi2/(8.*m_FV2);
              double g12g2g3 = 0.;
              double g2 = 3./(192.*sqr(M_PI)) * m_V[0].Mass()/(sqrt(2.)*m_FV);
             
              Complex FVchi  = 3./(2.*sqr(M_PI)*m_fpi2); 
              Complex FV1r1  = 6.*m_GV/(sqrt(2.)*m_fpi2*m_V[0].Mass()) * (
                  1./Complex(m_MO2-s,-1.*m_MO*m_GO) * (
                    c1c2c5*Q2 - c1c2c52c6*s + c1c28c3c5*Mass2(1) ) +
                  1./Complex(m_V[1].Mass2()-t,-1.*MassWidthVector(1,t)) * (
                    c1c2c5*Q2 - c1c2c52c6*t + c1c28c3c5*Mass2(0) + 8.*c4*(Mass2(0)-Mass2(1)) )
                  );
              Complex FV1r2  = -12.*m_FV/(sqrt(2.)*m_V[0].Mass()*m_fpi2) * 
                1./Complex(m_V[0].Mass2()-Q2,-1*MassWidthVector(0,Q2)) * (
                  g12g2g3*(s+t)-2.*g2*Q2 );
              Complex FV2r   = -6.*m_FV*m_GV/m_fpi2 * 
                1./Complex(m_V[0].Mass2()-Q2,-1.*MassWidthVector(0,Q2)) * (
                    1./Complex(m_MO2-s,-1.*m_MO*m_GO) * (
                      d3*(Q2+s) + d18d2d3*Mass2(1) ) +
                    1./Complex(m_V[1].Mass2()-t,-1.*MassWidthVector(1,t)) * (
                      d3*(Q2+t) + d18d2d3*Mass2(0) )
                    );
              return FVchi + FV1r1 + FV1r2 + FV2r;
            }
        }
      }
  }
  return  (j>3)? Complex(0.,0.) : Complex(1.,0.); 
}

double VA_0_PPP::RChT::FFunc( double a, double b, double c)
{
  return m_l2 + m_l1*a/c - m_l0*b/c;
}


// Parameterisation
// DECKER, FINKEMEIER, MIRKES hep-ph/9310270

VA_0_PPP::KS::KS(int mode, int kaon_mode, string path, GeneralModel _md, double * _ms)
  : FF_Base(mode,kaon_mode,path,_md,_ms)
{
  bool anomaly (false);
   
  // special constants for this particular parameterisation
  switch( m_mode ) {
    case 3000 : /* pi- pi- pi+ mode */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      break;
    case 1200 : /* pi0 pi0 pi- mode */
      m_X123  = Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      break;
    case 1020 : /* K- pi- K+ */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      anomaly = true;
      break;
    case 1002 : /* K0 pi- K0b */
      m_X123  = Mass2(1) + Mass2(0);
      m_ms123 = Mass2(1);
      m_G123  = 1;
      anomaly = true;
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
    case 2010 : /* K- pi- pi+ */
      m_X123  = Mass2(1)+Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      anomaly = true;
      break;
    case 1101 : /* pi- K0b pi0 */
      m_X123  = 0.;
      m_ms123 = Mass2(1);
      m_G123  = 0;
      anomaly = true;
      break;
    case   30 : /* K- K- K+ */
      m_X123  = 2.*Mass2(2);
      m_ms123 = Mass2(2);
      m_G123  = 1;
      break;
    case   12 : /* K- K0b K0 */
      m_X123  = 2.*Mass2(0);
      m_ms123 = Mass2(0);
      m_G123  = 1;
      break;
    default:
      msg_Error()<<METHOD<<" Warning: Three pseudoscalar mode not recognized. "
        <<"mode="<<m_mode<<endl;
      break;
  }
   
  kf_code resAnoV   = (m_deltas)? kf_K_star_892_plus : kf_rho_770_plus;
  kf_code resAnoVV  = (m_deltas)? kf_K_star_1410_plus : kf_rho_1450_plus;
  kf_code resAnoVVV = (m_deltas)? kf_K_star_1680_plus : kf_rho_1700_plus;

  if( anomaly ) {
	int running    = int( _md("RUNNING_WIDTH", 3 ) );     // running width
    double MV   = _md("Mass_anomaly_"+Flavour(resAnoV).IDName(),   Flavour(resAnoV).HadMass()   );   // mass V
    double MVV  = _md("Mass_anomaly_"+Flavour(resAnoVV).IDName(),  Flavour(resAnoVV).HadMass()   );  // mass V'
    double MVVV = _md("Mass_anomaly_"+Flavour(resAnoVVV).IDName(), Flavour(resAnoVVV).HadMass()   ); // mass V'
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
    m_AnoV      = ResonanceFlavour( kf_rho_770_plus, 0., 0., 0 );
    m_AnoVV     = ResonanceFlavour( kf_rho_770_plus, 0., 0., 0 );
    m_AnoVVV    = ResonanceFlavour( kf_rho_770_plus, 0., 0., 0 );
    m_AlphaV    = 0.;
    m_BetaV[0]  = 0.;
    m_BetaV[1]  = 0.;
  }

} 
 
// methods for axial and scalar FF
 
Complex VA_0_PPP::KS::BW_V( int a, double s )
{
  if( !m_G123 && a==1 ) return Complex(0.,0.);          // BW_V23 = 0 if G123=0
  return m_V[a].BreitWigner(s);
}
 
Complex VA_0_PPP::KS::BW_VV( int a, double s )
{
  if( !m_G123 && a==1 ) return Complex(0.,0.);          // BW_V23' = 0 if G123=0
  return m_VV[a].BreitWigner(s);
}

Complex VA_0_PPP::KS::BW_A( double s )
{
  if (!IsZero(m_alpha)) {
    return( ( m_A.BreitWigner(s)+m_alpha*m_AA.BreitWigner(s) ) / (1.+m_alpha) );
  }
  return m_A.BreitWigner(s);
}


Complex VA_0_PPP::KS::Tvector1( int a, int b, double x )
{
  Complex ret =           
                BW_V(a,x) *( 1.-(Mass2(a)-Mass2(b))/(3.*m_V[a].Mass2()) )
  + m_Beta[a]*( BW_VV(a,x)*( 1.-(Mass2(a)-Mass2(b))/(3.*m_VV[a].Mass2()) ) );
  return ret/(1.+m_Beta[a]);
}

Complex VA_0_PPP::KS::Tvector2( int a, int b, double x )
{
  Complex ret = BW_V(a,x)/m_V[a].Mass2() + m_Beta[a]*BW_VV(a,x)/m_V[a].Mass2();
  return ret*( 2.*(Mass2(a)-Mass2(b)) )/( 3.*(1.+m_Beta[a]) );    
}

Complex VA_0_PPP::KS::TSvector( int a, int b, int c, double Q2, double s, double t )
{
  Complex ret =
                BW_V(a,s)  * (   m_ms123*( Q2-2.*t-s+2.*Mass2(a)+Mass2(c) )
                               - (Mass2(a)-Mass2(b))/m_V[a].Mass2()*( m_ms123*(Q2+s-Mass2(c))-Q2*(s-m_V[a].Mass2()) ) )
    + m_Beta[a]*BW_VV(a,s) * (   m_ms123*( Q2-2.*t-s+2.*Mass2(a)+Mass2(c) ) 
                               - (Mass2(a)-Mass2(b))/m_VV[a].Mass2()*( m_ms123*(Q2+s-Mass2(c))-Q2*(s-m_VV[a].Mass2()) ) );
  return ret/(1.+m_Beta[a]);
}

Complex VA_0_PPP::KS::Tgen( int a, int b, int c, double s, double t)
{
  return( Tvector1(a-1,c-1,s) + Tvector2(b-1,c-1,t) );
}

// methods for vector FF
 
Complex VA_0_PPP::KS::T_V( double q2 )
{
  Complex ret(0.,0.);
  ret  =           m_AnoV.BreitWigner(q2)
    + m_BetaV[0] * m_AnoVV.BreitWigner(q2)
    + m_BetaV[1] * m_AnoVVV.BreitWigner(q2);
  ret /= 1.+m_BetaV[0]+m_BetaV[1];
  return ret;
}

Complex VA_0_PPP::KS::T_V13V23( double s, double t )
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

Complex VA_0_PPP::KS::FormFactor( int j, double Q2, double s, double t )
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

// Parameterisation
// FINKEMEIER, MIRKES hep-ph/9503474

VA_0_PPP::KS95::KS95(int mode, int kaon_mode, string path, GeneralModel _md, double * _ms)
  : FF_Base(mode,kaon_mode,path,_md,_ms)
{
  kf_code resV, resVV, resVVV;
  int running    = int( _md("RUNNING_WIDTH", 3 ) );     // running width
  double MV, MVV, MVVV, GV, GVV, GVVV;
  if( m_deltas==0 ) {
	resV   = kf_rho_770_plus;
	resVV  = kf_rho_1450_plus;
	resVVV = kf_rho_1700_plus;
	MV   = _md("q2_Mass_"+Flavour(resV).IDName(),   0.7769 );   // mass V
	MVV  = _md("q2_Mass_"+Flavour(resVV).IDName(),  1.363 );  // mass V'
	MVVV = _md("q2_Mass_"+Flavour(resVVV).IDName(), 1.700 ); // mass V"
	GV   = _md("q2_Width_"+Flavour(resV).IDName(),   0.1490 );   // width V
	GVV  = _md("q2_Width_"+Flavour(resVV).IDName(),  0.310 );  // width V'
	GVVV = _md("q2_Width_"+Flavour(resVVV).IDName(), 0.235 ); // width V"
	m_beta_VectorFF   = _md("q2_beta_"+Flavour(resVV).IDName(), -0.108 );        // beta for V'
	m_gamma_VectorFF  = _md("q2_gamma_"+Flavour(resVVV).IDName(), 0.020 );      // gamma for V"
  } 
  else {
	resV   = kf_K_star_892_plus;
	resVV  = kf_K_star_1410_plus;
	resVVV = kf_K_star_1680_plus;
	MV   = _md("q2_Mass_"+Flavour(resV).IDName(),   0.892    );   // mass V
	MVV  = _md("q2_Mass_"+Flavour(resVV).IDName(),  1.412    );  // mass V'
	MVVV = _md("q2_Mass_"+Flavour(resVVV).IDName(), 1.714    ); // mass V"
	GV   = _md("q2_Width_"+Flavour(resV).IDName(),   0.051   );   // width V
	GVV  = _md("q2_Width_"+Flavour(resVV).IDName(),  0.227   );  // width V'
	GVVV = _md("q2_Width_"+Flavour(resVVV).IDName(), 0.323   ); // width V"
	m_beta_VectorFF   = _md("q2_beta_"+Flavour(resVV).IDName(), -0.25 );        // beta for V'
	m_gamma_VectorFF  = _md("q2_gamma_"+Flavour(resVVV).IDName(), 0.038 );      // gamma for V"
  }

  m_V_VectorFF      = ResonanceFlavour( resV, MV, GV, running&2 );
  m_VV_VectorFF     = ResonanceFlavour( resVV, MVV, GVV, running&2 );
  m_VVV_VectorFF    = ResonanceFlavour( resVVV, MVVV, GVVV, running&2 );

  // omega - phi - resonances
  kf_code resO = kf_omega_782;
  kf_code resP = kf_phi_1020;
  double MO   = _md("G3_Mass_"+Flavour(resO).IDName(),   0.782   );    // mass omega
  double MP   = _md("G3_Mass_"+Flavour(resP).IDName(),   1.020   );    // mass phi
  double GO   = _md("G3_Width_"+Flavour(resO).IDName(),   0.00843   );   // mass omega
  double GP   = _md("G3_Width_"+Flavour(resP).IDName(),   0.00443   );   // mass phi
  m_Omega_VectorFF	= ResonanceFlavour( resO, MO, GO, 0 );
  m_Phi_VectorFF	= ResonanceFlavour( resP, MP, GP, 0 );
  m_eps_VectorFF    = _md("G3_epsilon_"+Flavour(resP).IDName(), 0.05 );        			// epsilon for phi

} 
 
Complex VA_0_PPP::KS95::Axial(double s)
{
  if (!IsZero(m_alpha)) {
    return( ( m_A.BreitWigner(s)+m_alpha*m_AA.BreitWigner(s) ) / (1.+m_alpha) );
  }
  return m_A.BreitWigner(s);
}
 
Complex VA_0_PPP::KS95::TVector(int i, double s)
{
  return (m_V[i].BreitWigner(s)+m_Beta[i]*m_VV[i].BreitWigner(s))/(1.+m_Beta[i]);
}
 
Complex VA_0_PPP::KS95::Vector(int i, double s, double u)
{
  switch(m_mode) {
	case 1101 : { 	// pi- K0b pi0
				  if(m_kaon_mode==3) {			// K0 K0
					switch(i) {
					  case 0 : return TVector(0,s);
					  case 1 : return TVector(1,s);
					}
				  }
				  else if(m_kaon_mode==20) { 	// KS KS
					switch(i) {
					  case 0 : return 0.5*TVector(1,u);
					  case 1 : return -0.5*(TVector(1,s)+TVector(1,u));
					}
				  }
				  else if(m_kaon_mode==200) {	// KL KL
					switch(i) {
					  case 0 : return 0.5*TVector(1,u);
					  case 1 : return -0.5*(TVector(1,s)+TVector(1,u));
					}
					switch(i) {
					  case 0 : return -0.5*TVector(1,u);
					  case 1 : return 0.5*(TVector(1,s)+TVector(1,u));
					}
				  }
				  else if(m_kaon_mode==110) {	// KS KL
					switch(i) {
					  case 0 : return 0.5*(2.*TVector(0,s)+TVector(1,u));
					  case 1 : return 0.5*(TVector(1,s)-TVector(1,u));
					}
				  }
				}
	case  111 : {	// K- pi0 K0
				  switch(i) {
					case 0 : return 2./3.*TVector(0,s) + 1./3.*TVector(1,u);
					case 1 : return 1./3.*TVector(1,s) - 1./3.*TVector(1,u);
				  }
				}
  }
  return TVector(i,s);
}
 
Complex VA_0_PPP::KS95::TOmega(double s)
{
  return (m_Omega_VectorFF.BreitWigner(s) + m_eps_VectorFF*m_Phi_VectorFF.BreitWigner(s)) / (1.+m_eps_VectorFF);
}
 
Complex VA_0_PPP::KS95::G3(double q2, double s, double t, double u )
{
  Complex vector1 
	=                     m_V_VectorFF.BreitWigner(q2)
	   +  m_beta_VectorFF*m_VV_VectorFF.BreitWigner(q2)
	   + m_gamma_VectorFF*m_VVV_VectorFF.BreitWigner(q2);
  vector1 /= 1.+m_beta_VectorFF+m_gamma_VectorFF;	   
	    
  Complex vector(0.,0.);
  switch(m_mode) {
	case 3000 :
	case 1200 : return Complex(0.,0.);
	case  210 : // pi0 pi0 K-
				{
				  vector  = TVector(2-1,t) - TVector(1-1,s);
				}
	case 2010 : // K- pi- pi+
				{
				  vector  = TVector(2-1,t) + TVector(1-1,s);
				}
	case 1101 : // pi- K0b pi0
				{
				  vector  = 2.*TVector(1-1,s) + TVector(2-1,t) + TVector(2-1,u);
				  vector /= 3.;
				}
	case 1020 : // K- pi- K+
				{
				  vector  = sqrt(2.)*TOmega(s) + TVector(2-1,t);
				  vector *= (sqrt(2.)-1.);
				}
	case 1002 : // K0 pi- K0b
				{
				  if(m_kaon_mode==2) {	// K0 K0
					vector  = sqrt(2.)*TOmega(s) + TVector(2-1,t);
					vector *= (sqrt(2.)-1.);
				  }
				  else if(m_kaon_mode==20) {	// KS KS
					vector  = TVector(2-1,t)-TVector(2-1,u);
					vector *= -0.5*(sqrt(2.)-1.);
				  }
				  else if(m_kaon_mode==200) {	// KL KL
					vector  = TVector(2-1,t)-TVector(2-1,u);
					vector *= -0.5*(sqrt(2.)-1.);
				  }
				  else if(m_kaon_mode==110)	{	// KS KL
					vector  = 2.*sqrt(2.)*TOmega(s)+TVector(2-1,t)+TVector(2-1,u);
					vector *= 0.5*(sqrt(2.)-1.);
				  }
				}
	case  111 : // K- pi0 K0
				{
				  vector  = TVector(1-1,s) - TVector(2-1,t);
				  vector *= (sqrt(2.)-1.);
				}
  }
  return vector1*vector;
}

// handling of form factors

Complex VA_0_PPP::KS95::FormFactor( int j, double Q2, double s, double t )
{
  double u = Q2 + Mass2(1-1) + Mass2(2-1) + Mass2(3-1) - s - t;
  switch( j ) {
    case 1 : 	// axial
		return Axial(Q2)*Vector(1-1,s,u);
    case 2 : 	// axial
		return Axial(Q2)*Vector(2-1,t,u);
    case 3 : 	// scalar
		break;
    case 4 :  	// vector
		return 6./(8.*sqr(M_PI)*m_fpi2)*G3(Q2,s,t,u);
		break; 
  }
  Complex FF(0.,0.);
  return FF;
}

DEFINE_CURRENT_GETTER(VA_0_PPP,"VA_0_PPP")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PPP>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi \\pi K $ \n\n"
    <<"Order: special, see tau decays. \n\n"
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 1 :} Kuehn-Santamaria \n"
    <<"    \\item {\\tt FORM\\_FACTOR = 2 :} Resonance Chiral Theory \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
