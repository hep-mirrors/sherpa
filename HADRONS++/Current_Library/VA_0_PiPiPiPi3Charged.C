#include "HADRONS++/Current_Library/VA_0_PiPiPiPi3Charged.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_File.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void VA_0_PiPiPiPi3Charged::SetModelParameters( GeneralModel _md )
{
  double Vud = _md("Vud", Tools::Vud);
  m_global   = Vud;

  p_lorenz = new Novo(m_path,_md);
}

void VA_0_PiPiPiPi3Charged::LorenzBase::SetPrivates(
  const ATOOLS::Vec4D * _p,int * i)
{
  p_i = i;
  p_p  = _p;
  m_q = p_p[p_i[0]]+p_p[p_i[1]]+p_p[p_i[2]]+p_p[p_i[3]];         // = q
  m_q2 = m_q.Abs2();
}

// Novosibirsk parameterisation
// see hep-ph/0201149 for details and
// hep-ph/0312240 for errata

VA_0_PiPiPiPi3Charged::Novo::Novo( string path, GeneralModel _md )
  : LorenzBase()
{
  m_mpi2 = sqr( Flavour(kf_pi_plus).HadMass() );
  m_rho = SimpleResonanceFlavour(	  
	  string("kf_rho_770_plus"),
      _md("Mass_rho(770)+", _md("Mass_rho(770)", 0.7761 )),
      _md("Width_rho(770)+", _md("Width_rho(770)", 0.1445 ))
	  );
  m_sigma = SimpleResonanceFlavour(
      string("kf_f_0_600"),
      _md("Mass_sigma", 0.800),
	  _md("Width_sigma", 0.800)
	  );
  m_omega = SimpleResonanceFlavour(
      string("kf_omega_782"),
      _md("Mass_omega(782)", 0.782),
	  _md("Width_omega(782)", 0.00841)
	  );
  m_a1 = SimpleResonanceFlavour(
      string("a(1)(1260)+"),
	  _md("Mass_a(1)(1260)+", 1.230 ),
	  _md("Width_a(1)(1260)+", 0.450 )
  );

  m_Lambda2 = _md("Lambda2", 1.2);
  double z_abs = _md("z_abs", 1.3998721 );
  double z_phase = _md("z_phase", 0.43585036 );
  m_z = Complex( z_abs*cos(z_phase), z_abs*sin(z_phase) );

  // G(q2) function (read as histogram)
  My_In_File fG("",path+"PhaseSpaceFunctions/G_pi-pi+pi-pi0.dat");
  if( !fG.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"PhaseSpaceFunctions/G_pi-pi+pi-pi0.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi3Charged::Novo::Novo(...) \n"
                  <<"     Read G_{pi-pi+pi-pi0}(q2)."<<endl;
    std::string found_file_name = fG.Path()+fG.File();
    fG.Close();
    p_G = new Histogram( found_file_name );
  }

  // Gomega(q2) function (read as histogram)
  My_In_File fGo("",path+"PhaseSpaceFunctions/Gomega_pi-pi+pi-pi0.dat");
  if( !fGo.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"/PhaseSpaceFunctions/Gomega_pi+pi0pi-pi-.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi3Charged::Novo::Novo(...) \n"
                  <<"     Read Gomega_{pi+pi0pi-pi-}(q2)."<<endl;
    std::string found_file_name = fGo.Path()+fGo.File();
    fGo.Close();
    p_Go = new Histogram( found_file_name );
  }

  // RunningWidth_a1(q2) function (read as histogram)
  My_In_File fW("",path+"PhaseSpaceFunctions/RunningWidth_a1_4pi-channel.dat");
  if( !fW.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"/PhaseSpaceFunctions/RunningWidth_a1_4pi-channel.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi3Charged::Novo::Novo(...) \n"
                  <<"     Read a1's running width (q2)."<<endl;
    std::string found_file_name = fW.Path()+fW.File();
    fW.Close();
    p_a1width = new Histogram( found_file_name );
  }
}

// global form factor
double VA_0_PiPiPiPi3Charged::Novo::G( double q2 )
{
  if(sqrt(q2)<0.6) return 0.;
  double val(0.);
  p_G->Extrapolate( sqrt(q2)+p_G->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  return val*76.565*sqrt(0.7179*sqrt(q2)-0.27505)/(sqr(m_rho.Mass2())*sqrt(q2));
  		// note: magic correction from TAUOLA (hep-ph/0312240)
}

double VA_0_PiPiPiPi3Charged::Novo::Go( double q2 )
{
  if(sqrt(q2)<0.6) return 0.;
  double val(0.);
  p_Go->Extrapolate( sqrt(q2)+p_Go->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  return val*886.84*sqrt(0.70983*sqrt(q2)-0.26689)/(sqr(m_rho.Mass2())*sqrt(q2));
  		// note: magic correction from TAUOLA (hep-ph/0312240)
}

// a1 form factor
double VA_0_PiPiPiPi3Charged::Novo::F2_a1( double s )
{
  return sqr( (1.+m_a1.Mass2()/m_Lambda2)/(1.+s/m_Lambda2) );
}

// a1 propagator
Complex VA_0_PiPiPiPi3Charged::Novo::Da1( double s )
{
  double val(0);
  p_a1width->Extrapolate( s+p_a1width->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  return Complex( s/m_a1.Mass2()-1., sqrt(s)*val/m_a1.Mass2() );
}

// rho h function
double VA_0_PiPiPiPi3Charged::Novo::hrho( double s )
{
  double epsilon = 1.e-8;
  if( s>4*m_mpi2 ) {
	double root = sqrt(1.-4.*m_mpi2/s);
	return root*log((1.+root)/(1.-root))*(s-4.*m_mpi2)/M_PI;
  }
  else if( s>epsilon ) return 0.;
  else return -8.*m_mpi2/M_PI;
}

// rho propagator
Complex VA_0_PiPiPiPi3Charged::Novo::Drho( double s )
{
  double real, imag;
  double g2 = pow( m_rho.Mass2()-4.*m_mpi2, 1.5 )/sqrt(m_rho.Mass2());
  double root = sqrt(1.-4.*m_mpi2/m_rho.Mass2());
  double dhrhods = root/M_PI  * ( root + (1.+2.*m_mpi2/m_rho.Mass2())*log((1.+root)/(1.-root)) );
  double dm = ( hrho(s) - hrho(m_rho.Mass2()) - (s-m_rho.Mass2())*dhrhods ) / g2;
  real = s-m_rho.Mass2() - m_rho.Mass()*m_rho.Width()*dm;
  if( s-4*m_mpi2 <= 0. ) imag = 0.;
  else {
	double g1 = pow( s-4.*m_mpi2, 1.5 )/sqrt(s);
	imag = m_rho.Mass()*m_rho.Width()*g1/g2;
  }
  double dm0 = ( hrho(0.) - hrho(m_rho.Mass2()) + m_rho.Mass2()*dhrhods ) / g2;
  double normalisation = m_rho.Mass2() + m_rho.Mass()*m_rho.Width()*dm0;
  return Complex(real/normalisation,imag/normalisation);
}
 
// sigma propagator
Complex VA_0_PiPiPiPi3Charged::Novo::Dsigma( double s )
{
  double real, imag;
  real = s-m_sigma.Mass2();
  if( s-4*m_mpi2 <= 0. ) imag = 0.;
  else {
	double g1 = sqrt( 1.-4.*m_mpi2/s );
	double g2 = sqrt( 1.-4.*m_mpi2/m_sigma.Mass2() );
	imag = m_sigma.Mass()*m_sigma.Width()*g1/g2;
  }
  double normalisation = m_sigma.Mass2();
  return Complex(real/normalisation,imag/normalisation);
}
 
// omega form factor
double VA_0_PiPiPiPi3Charged::Novo::F2_o( double s )
{
  return 1.;
}

// omega propagator
Complex VA_0_PiPiPiPi3Charged::Novo::Do( double s )
{
  double real, imag;
  real = s-m_omega.Mass2();
  double q (sqrt(s)), g;
  if( q<1.0 ) {
	q -= m_omega.Mass();
	g = 1+q*(17.560+q*(141.110+q*(894.884+q*(4977.35+q*(7610.66-42524.4*q)))));
	if( g<0. ) g = 0.;
  }
  else {
	g = -1333.26+q*(4860.19+q*(-6000.81+2504.97*q));
  }
  imag = m_omega.Mass()*m_omega.Width()*g;
  double normalisation = m_omega.Mass2();
  return Complex(real/normalisation,imag/normalisation);
}
 
// t1 function
Vec4C VA_0_PiPiPiPi3Charged::Novo::t1( int a, int b, int c, int d )
{
  // subtract 1 for C++ convention
  a--; b--; c--; d--;

  Vec4D qa1 = m_q - p_p[p_i[a]]; 			// virtuality of a1
  double qa12 = qa1.Abs2();
  Complex global = Complex(-1.*F2_a1(qa12), 0.) / Da1(qa12) / Drho((p_p[p_i[c]]+p_p[p_i[d]]).Abs2());
  double dot0 = qa1*m_q;
  double dot1 = qa1*p_p[p_i[c]];
  double dot2 = qa1*p_p[p_i[d]];
  double dot3 = (m_q*p_p[p_i[d]]) * (p_p[p_i[a]]*p_p[p_i[c]]);
  double dot4 = (m_q*p_p[p_i[c]]) * (p_p[p_i[d]]*p_p[p_i[a]]);
  return global * ( dot0*(p_p[p_i[d]]*dot1 - p_p[p_i[c]]*dot2) + qa1*(dot3-dot4) );
}

// t2 function
Vec4C VA_0_PiPiPiPi3Charged::Novo::t2( int a, int b, int c, int d )
{
  // subtract 1 for C++ convention
  a--; b--; c--; d--;

  Vec4D qa1 = m_q - p_p[p_i[a]];			 // virtuality of a1
  double qa12 = qa1.Abs2();
  Complex global = m_z * F2_a1(qa12) / Da1(qa12) / Dsigma((p_p[p_i[c]]+p_p[p_i[d]]).Abs2());
  double dot0 = m_q*qa1;
  double dot1 = m_q*p_p[p_i[b]];
  return global * ( p_p[p_i[b]]*dot0 - qa1*dot1 ) * qa12;
}

// t3 function
Vec4C VA_0_PiPiPiPi3Charged::Novo::t3( int a, int b, int c, int d )
{
  // subtract 1 for C++ convention
  a--; b--; c--; d--;

  Vec4D qo = m_q - p_p[p_i[a]];			 // virtuality of omega
  double qo2 = qo.Abs2();
  Complex global = F2_o(qo2) / Do(qo2) / Drho((p_p[p_i[c]]+p_p[p_i[d]]).Abs2());
  double dotq2 = m_q*p_p[p_i[b]];
  double dotq3 = m_q*p_p[p_i[c]];
  double dotq4 = m_q*p_p[p_i[d]];
  double dot12 = p_p[p_i[a]]*p_p[p_i[b]];
  double dot13 = p_p[p_i[a]]*p_p[p_i[c]];
  double dot14 = p_p[p_i[a]]*p_p[p_i[d]];

  return global * ( 
	    p_p[p_i[b]]*(dotq3*dot14-dotq4*dot13) 
	  - dotq2*(p_p[p_i[c]]*dot14-p_p[p_i[d]]*dot13) 
	  + dot12*(p_p[p_i[c]]*dotq4-p_p[p_i[d]]*dotq3) );
}

Vec4C VA_0_PiPiPiPi3Charged::Novo::operator()()
{
  // a1 -> rho pi current
  Vec4C J_a1rpi = 
	G(m_q2) * ( 
		  t1(1,2,3,4)
		+ t1(3,2,1,4)
		+ t1(1,3,2,4)
		+ t1(3,1,2,4)
		+ t1(4,3,1,2)
		+ t1(4,1,3,2) );
  // a1 -> sigma pi current
  Vec4C J_a1spi =
	G(m_q2) * (
		  t2(4,3,1,2)
		+ t2(4,1,3,2)
		- t2(1,4,3,2)
		- t2(3,4,1,2) );
  // omega -> rho pi current
  Vec4C J_omrpi =
	Go(m_q2) * (
		  t3(1,2,3,4)
		+ t3(3,2,1,4)
		- t3(1,3,2,4)
		- t3(3,1,2,4)
		- t3(1,4,3,2)
		- t3(3,4,1,2) );

  // sum (factor 2! due to 2 identical particles)
  return (J_a1rpi+J_a1spi+J_omrpi)/sqrt(2.);
}

// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

VA_0_PiPiPiPi3Charged::KS::KS( string path, GeneralModel _md )
  : LorenzBase()
{
  for (int i=1; i<=3; i++ ) {
    m_r[i] = m_q-p_p[p_i[i]];
    m_s[i] = (p_p[p_i[0]]+p_p[p_i[i]]).Abs2();
  }
  // redefinition of variables
  m_s[0] = (p_p[p_i[1]]+p_p[p_i[3]]).Abs2();          // = t4
  // unused variables yet
  m_r[0] = ATOOLS::Vec4D(0.,0.,0.,0.);
   
  m_fpi2    = sqr(_md("fpi",0.0924))*2.;        // redefine fpi
  m_grop    = _md("grop", 12.924);
  m_Go3p    = _md("Go3p", 1476.);

  m_mpi2    = sqr( Flavour(kf_pi_plus).HadMass() );
  m_mpi02   = sqr( Flavour(kf_pi).HadMass() );

  double MR      = _md("Mass_rho(770)+",  Flavour(kf_rho_770_plus).HadMass()  );
  double MRR     = _md("Mass_rho(1450)+", Flavour(kf_rho_1450_plus).HadMass() );
  double MRRR    = _md("Mass_rho(1700)+", Flavour(kf_rho_1700_plus).HadMass() );
  double MO      = _md("Mass_omega(782)", Flavour(kf_omega_782).HadMass() );
  double MF      = _md("Mass_f(0)(980)",    Flavour(kf_f_0_980).HadMass() );
  double MS      = _md("Mass_sigma",      Flavour(kf_f_0_980).HadMass()  );
  double MA      = _md("Mass_a(1)(1260)+",  Flavour(kf_a_1_1260_plus).HadMass());
  double GR      = _md("Width_rho(770)+",  Flavour(kf_rho_770_plus).Width()  );
  double GRR     = _md("Width_rho(1450)+", Flavour(kf_rho_1450_plus).Width() );
  double GRRR    = _md("Width_rho(1700)+", Flavour(kf_rho_1700_plus).Width() );
  double GO      = _md("Width_omega(782)", Flavour(kf_omega_782).Width() );
  double GF      = _md("Width_f(0)(980)",    Flavour(kf_f_0_980).Width() );
  double GS      = _md("Width_sigma",      Flavour(kf_f_0_980).Width()  );
  double GA      = _md("Width_a(1)(1260)+",  Flavour(kf_a_1_1260_plus).Width());
  m_Rho = ResonanceFlavour( kf_rho_770_plus, MR, GR, 1 );
  m_RR  = ResonanceFlavour( kf_rho_1450_plus, MRR, GRR, 1 );
  m_RRR = ResonanceFlavour( kf_rho_1700_plus, MRRR, GRRR, 1 );
  m_O   = ResonanceFlavour( kf_omega_782, MO, GO, 0 );
  m_F   = ResonanceFlavour( kf_f_0_980, MF, GF, 1 );
  m_S   = ResonanceFlavour( kf_f_0_980, MS, GS, 1 );
  m_A   = ResonanceFlavour( kf_a_1_1260_plus, MA, GA, 1, path );
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
  m_sigma   = _md("sigma", 0.);
  m_A.InitialiseThreeBodyResonance(m_Rho, m_RR, m_beta);
   
   
  m_Frho    = _md("frho", 0.266)*m_Rho.Mass2();

//   m_R[0]=m_R[1] = 0.;
//   m_R[2]        = -2.;
//   m_R[3]=m_R[4] = 1.;
  m_R[0] = 0.;
  m_R[1] = -2.;
  m_R[2]=m_R[3] = 1.;
  
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

Complex VA_0_PiPiPiPi3Charged::KS::Fk( double x, Complex * _beta )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (_beta[0] + _beta[1]*BW_R + _beta[2]*BW_RR + _beta[3]*BW_RRR)/(_beta[0]+_beta[1]+_beta[2]+_beta[3]) );
}

Complex VA_0_PiPiPiPi3Charged::KS::Trho( double x )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Complex VA_0_PiPiPiPi3Charged::KS::TTrho( double x )
{
  Complex BW_R   = m_Rho.BreitWignerAlt(x);
  Complex BW_RR  = m_RR.BreitWignerAlt(x);
  return( BW_R + m_sigma*BW_RR );
}

double VA_0_PiPiPiPi3Charged::KS::Dots( int k, int l )        // k=2,3; l=0,1,2,3 (l!=k)
  // special dot products used for anomalous part of J_omega
{
  int pre = l-1;                // predecessor of l
  if (pre==-1) pre  = 3;
  if (pre==k) pre -= 1;     
  int suc  = l+1;               // successor of l
  if (suc==k)  suc += 1;
  if (suc==4)  suc  = 0;
  return( (m_r[k]*p_p[p_i[pre]])*(p_p[p_i[suc]]*p_p[p_i[k]]) 
        - (m_r[k]*p_p[p_i[suc]])*(p_p[p_i[pre]]*p_p[p_i[k]]) );
}

Vec4C VA_0_PiPiPiPi3Charged::KS::OmegaPi()
{
  // chiral part
  Vec4C J_chi( 0.0,0.0,0.0,0.0 );
  Vec4C sum1( 0.0,0.0,0.0,0.0 );
  Vec4C sum2( 0.0,0.0,0.0,0.0 );
  for (int k=1; k<=3; k++) {
    sum2 = Vec4C( 0.0,0.0,0.0,0.0 );
    for (int l=1; l<=3; l++) if (l!=k) {
      sum2 += Vec4C(m_q-2.*p_p[p_i[l]]) * ((m_r[l]*(p_p[p_i[k]]-p_p[p_i[0]]))/m_r[l].Abs2());
    }
    sum1 += m_R[k]*Trho(m_s[k]) * (Vec4C(p_p[p_i[k]])-Vec4C(p_p[p_i[0]])-sum2);
  }
  J_chi = 2.*sqrt(3.)/m_fpi2*Trho(m_q2) * sum1;

  // anomalous part
  sum1 = Vec4C( 0.0,0.0,0.0,0.0 );
  for( int k=2; k<=3; k++ ) {
    sum2 = Vec4C( 0.0,0.0,0.0,0.0 );
    for (int l=0; l<=3; l++) if (l!=k) {
      sum2 += Vec4C(p_p[p_i[l]]) * Dots(k,l);
    }
    sum1 += m_O.BreitWignerAlt(m_r[k].Abs2())*sum2;
  }
  Vec4C J_a = m_Go3p*m_Frho*m_grop*TTrho(m_q2) * sum1;

  // total
  return( (J_chi + J_a)*Fk(m_q2,m_Beta_opi) );
}

Vec4C VA_0_PiPiPiPi3Charged::KS::AonePi()
{
  Vec4C term1(0.,0.,0.,0.), term2(0.,0.,0.,0.);
  Complex A, B, C;
  Vec4D   P, Q, R;

  //  1st term
  P = p_p[p_i[0]]-p_p[p_i[2]];
  Q = p_p[p_i[0]]-p_p[p_i[3]];
  R = m_r[1];
  A = m_Rho.BreitWigner(m_s[2]);
  B = m_Rho.BreitWigner(m_s[3]);
  C = m_A.BreitWigner(R.Abs2());
  term1  = A*Vec4C(p_p[p_i[0]]-p_p[p_i[2]])
         + B*Vec4C(p_p[p_i[0]]-p_p[p_i[3]])
	 - Vec4C(m_q-p_p[p_i[1]])*(( A*(R*P) + B*(R*Q) )/R.Abs2())
	 - Vec4C(m_q)*(( A*(m_q*P)
			      + B*(m_q*Q)
			      - (m_q*Q)*(A*(R*P)+B*(R*Q))/R.Abs2()
	                      )/m_q2);
  term1 *= C;

  // 2nd and 3rd term
  int ind;
  Vec4C help;
  for (int k=2; k<=3; k++) {
    ind = (k==2)? 3 : 2;
    P = p_p[p_i[0]]-p_p[p_i[1]];
    Q = p_p[p_i[ind]]-p_p[p_i[1]];
    R = m_r[k];
    A = m_Rho.BreitWigner(m_s[1]);
    double t = (ind==2)? (p_p[p_i[1]]+p_p[p_i[2]]).Abs2() : m_s[0]; // s[0] = t[3]
    B = m_Rho.BreitWigner(t);
    C = m_A.BreitWigner(R.Abs2());
    help   = A*Vec4C(p_p[p_i[0]]-p_p[p_i[1]])
           + B*Vec4C(p_p[p_i[1]]-p_p[p_i[ind]])
	   - Vec4C(m_q-p_p[p_i[k]])*(( A*(R*P) + B*(R*Q) )/R.Abs2())
	   - Vec4C(m_q)*(( A*(m_q*P)
                               + B*(m_q*Q) 
                               - (m_q*Q)*(A*(R*P)+B*(R*Q))/R.Abs2()
	                        )/m_q2);
    term2 += C*help;
  }

  // total 
  return (term1-term2)*Fk(m_q2,m_Beta_api);
}

Vec4C VA_0_PiPiPiPi3Charged::KS::SigmaRho()
{
  int ind;
  Vec4C term( 0.0,0.0,0.0,0.0 );
  Complex BW_S, BW_R;
  for (int k=2; k<=3; k++) {
    ind = (k==2)? 3 : 2;
    BW_S = m_S.BreitWigner(m_s[k]);
    double t = (ind==2)? (p_p[p_i[1]]+p_p[p_i[2]]).Abs2() : m_s[0]; // s[0] = t[3]
    BW_R = m_Rho.BreitWigner(t);
    term += BW_S*BW_R * Vec4C( p_p[p_i[1]] - p_p[p_i[ind]] + (m_q*(p_p[p_i[ind]]-p_p[p_i[1]]))*m_q );
  }
  return term*Fk(m_q2,m_Beta_srh);
}

Vec4C VA_0_PiPiPiPi3Charged::KS::FzeroRho()
{
  int ind;
  Vec4C term( 0.0,0.0,0.0,0.0 );
  Complex BW_S, BW_R;
  for (int k=2; k<=3; k++) {
    ind = (k==2)? 3 : 2;
    BW_S = m_F.BreitWigner(m_s[k]);
    double t = (ind==2)? (p_p[p_i[1]]+p_p[p_i[2]]).Abs2() : m_s[0]; // s[0] = t[3]
    BW_R = m_Rho.BreitWigner(t);
    term += BW_S*BW_R * Vec4C( p_p[p_i[1]] - p_p[p_i[ind]] + (m_q*(p_p[p_i[ind]]-p_p[p_i[1]]))*m_q );
  }
  return term*Fk(m_q2,m_Beta_frh);
}

Vec4C VA_0_PiPiPiPi3Charged::KS::operator()()
{
  Vec4C help(0.,0.,0.,0.);
  help += m_Alpha[0] * OmegaPi();
  help += m_Alpha[1] * AonePi();
  help += m_Alpha[2] * SigmaRho();
  help += m_Alpha[3] * FzeroRho();
  return help/m_SumAlpha;
}

// General framework

void VA_0_PiPiPiPi3Charged::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4C help(0.,0.,0.,0.);
  // Lorentz structure
  p_lorenz->SetPrivates(&moms.front(), &p_i.front());
  help = (*p_lorenz)();
  Insert( help*m_global , 0);
}

DEFINE_CURRENT_GETTER(VA_0_PiPiPiPi3Charged,"VA_0_PiPiPiPi3Charged")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PiPiPiPi3Charged>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi^+ \\pi^+ \\pi^- \\pi^0 $ \n\n"
    <<"Order: 0 = $\\pi^\\mp$, 1 = $\\pi^0$, 2 = $\\pi^\\pm$, 3 = $\\pi^\\pm$ "
    <<"Available form factors: \n"
    <<"  \\begin{itemize} \n"
    <<"    \\item Kuehn-Santamaria (default) \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
