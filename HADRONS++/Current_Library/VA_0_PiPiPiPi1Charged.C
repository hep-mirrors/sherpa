#include "HADRONS++/Current_Library/VA_0_PiPiPiPi1Charged.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

void VA_0_PiPiPiPi1Charged::SetModelParameters( GeneralModel _md )
{ 
  double Vud = _md("Vud", Tools::Vud);
  m_global   = Vud;		

//  p_lorenz = new KS(_md);
  p_lorenz = new Novo(_md,m_path);
}
 
void VA_0_PiPiPiPi1Charged::LorenzBase::SetPrivates(
  const ATOOLS::Vec4D * _p, int* i)
{
  p_p  = _p;
  p_i  = i;
  m_q = p_p[p_i[0]]+p_p[p_i[1]]+p_p[p_i[2]]+p_p[p_i[3]];
  m_q2 = m_q.Abs2();
  // m_s[0] and m_r[0] unused
  for (int i=1; i<=3; i++ ) {
    m_r[i] = m_q-p_p[p_i[i]];
    m_s[i] = (p_p[p_i[0]]+p_p[p_i[i]]).Abs2();
  }
}
 
// Novosibirsk parameterisation
// see hep-ph/0201149 for details and
// hep-ph/0312240 for errata

VA_0_PiPiPiPi1Charged::Novo::Novo( GeneralModel _md, string path )
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
  My_In_File fG("",path+"PhaseSpaceFunctions/G_pi+pi0pi0pi0.dat");
  if( !fG.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"/PhaseSpaceFunctions/G_pi+pi0pi0pi0.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	Abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi1Charged::Novo::Novo(...) \n"
                  <<"     Read G_{pi+pi0pi0pi0}(q2)."<<endl;
    std::string found_file_name = fG.Path()+fG.File();
    fG.Close();
    p_G = new Histogram( found_file_name );
  }

  // RunningWidth_a1(q2) function (read as histogram)
  My_In_File fW("",path+"PhaseSpaceFunctions/RunningWidth_a1_4pi-channel.dat");
  if( !fW.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"/PhaseSpaceFunctions/RunningWidth_a1_4pi-channel.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	Abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi1Charged::Novo::Novo(...) \n"
                  <<"     Read a1's running width (q2)."<<endl;
    std::string found_file_name = fW.Path()+fW.File();
    fW.Close();
    p_a1width = new Histogram( found_file_name );
  }

  // magic correction function from TAUOLA (see hep-ph/0312240)
  My_In_File fF("",path+"PhaseSpaceFunctions/z_forma1_q2.dat");
  if( !fF.Open() ) {                            // if file does not exist
    msg_Error()<<"The file "<<path<<"/PhaseSpaceFunctions/z_forma1_q2.dat does"
               <<"not exist. Don't know what to do. Will abort."<<endl;
	Abort();		   
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::VA_0_PiPiPiPi1Charged::Novo::Novo(...) \n"
                  <<"     Read TAUOLA's magic correction function (q2)."<<endl;
    std::string found_file_name = fF.Path()+fF.File();
    fF.Close();
    p_zforma1 = new Histogram( found_file_name );
  }
}

// global form factor
double VA_0_PiPiPiPi1Charged::Novo::G( double q2 )
{
  if(sqrt(q2)<0.6) return 0.;
  double val(0.);
  p_G->Extrapolate( sqrt(q2)+p_G->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  double val2(0.);
  p_zforma1->Extrapolate( q2+p_zforma1->BinSize(), &val2, 0 );
                            // shift +BinSize to get right value
  return val*96.867*val2*sqrt(0.70907*sqrt(q2)-0.26413)/(sqr(m_rho.Mass2())*sqrt(q2));
  		// note: magic correction from TAUOLA (hep-ph/0312240)
}

// a1 form factor
double VA_0_PiPiPiPi1Charged::Novo::F2_a1( double s )
{
  return sqr( (1.+m_a1.Mass2()/m_Lambda2)/(1.+s/m_Lambda2) );
}

// a1 propagator
Complex VA_0_PiPiPiPi1Charged::Novo::Da1( double s )
{
  double val(0);
  p_a1width->Extrapolate( s+p_a1width->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  return Complex( s/m_a1.Mass2()-1., sqrt(s)*val/m_a1.Mass2() );
}

// rho h function
double VA_0_PiPiPiPi1Charged::Novo::hrho( double s )
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
Complex VA_0_PiPiPiPi1Charged::Novo::Drho( double s )
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
Complex VA_0_PiPiPiPi1Charged::Novo::Dsigma( double s )
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
 
// t1 function
Vec4C VA_0_PiPiPiPi1Charged::Novo::t1( int a, int b, int c, int d )
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
Vec4C VA_0_PiPiPiPi1Charged::Novo::t2( int a, int b, int c, int d )
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

Vec4C VA_0_PiPiPiPi1Charged::Novo::operator()()
{
  // a1 -> rho pi current
  Vec4C J_a1rpi = 
	G(m_q2) * ( 
		  t1(2,3,1,4)
		+ t1(2,4,1,3)
		+ t1(3,2,1,4)
		+ t1(3,4,1,2)
		+ t1(4,2,1,3)
		+ t1(4,3,1,2) );

  // a1 -> sigma pi current
  Vec4C J_a1spi =
	G(m_q2) * (
		  t2(2,1,3,4)
		+ t2(3,1,2,4)
		+ t2(4,1,3,2) 
		- t2(1,2,3,4)
		- t2(1,3,2,4)
		- t2(1,4,3,2) );

  // sum (factor 3! due to 3 identical particles)
  return (J_a1rpi+J_a1spi)/sqrt(6.);
}

// CLEO parameterisation
// see hep-ex/9908024 and CERN-TH.6793/93 for details

VA_0_PiPiPiPi1Charged::KS::KS( GeneralModel _md )
  : LorenzBase()
{
  m_fpi2    = sqr(_md("fpi",0.1307));        
  m_mpi2    = sqr( Flavour(kf_pi_plus).HadMass() );
  m_mpi02   = sqr( Flavour(kf_pi).HadMass() );

  double MR      = _md("Mass_rho(770)+",  Flavour(kf_rho_770_plus).HadMass()  );
  double MRR     = _md("Mass_rho(1450)+", Flavour(kf_rho_1450_plus).HadMass() );
  double MRRR    = _md("Mass_rho(1700)+", Flavour(kf_rho_1700_plus).HadMass() );
  double GR      = _md("Width_rho(770)+",  Flavour(kf_rho_770_plus).Width()  );
  double GRR     = _md("Width_rho(1450)+", Flavour(kf_rho_1450_plus).Width() );
  double GRRR    = _md("Width_rho(1700)+", Flavour(kf_rho_1700_plus).Width() );
  m_Rho = ResonanceFlavour( kf_rho_770_plus, MR, GR, 1 );
  m_RR  = ResonanceFlavour( kf_rho_1450_plus, MRR, GRR, 1 );
  m_RRR = ResonanceFlavour( kf_rho_1700_plus, MRRR, GRRR, 1 );
   
  m_beta    = _md("beta", 0.);
  m_gamma   = _md("gamma", 0.);
}

Complex VA_0_PiPiPiPi1Charged::KS::Trho( double x )
{
  Complex BW_R   = m_Rho.BreitWigner(x);
  Complex BW_RR  = m_RR.BreitWigner(x);
  Complex BW_RRR = m_RRR.BreitWigner(x);
  return( (BW_R + m_beta*BW_RR + m_gamma*BW_RRR)/(1.+m_beta+m_gamma) );
}

Vec4C VA_0_PiPiPiPi1Charged::KS::operator()()
{
  Vec4C J_chi;
  Vec4C sum1(0.,0.,0.,0.), sum2 (0.,0.,0.,0.);
  for (int k=1; k<=3; k++) {
    sum2 = Vec4C(0.,0.,0.,0.);
    for (int l=1; l<=3; l++) if (l!=k) {
      sum2 += Vec4C(m_q-2.*p_p[p_i[l]]) * ((m_r[l]*(p_p[p_i[k]]-p_p[p_i[0]]))/m_r[l].Abs2());
    }
    sum1 += Trho(m_s[k]) * (Vec4C(p_p[p_i[k]])-Vec4C(p_p[p_i[0]])-sum2);
  }
  J_chi = 2.*sqrt(3.)/m_fpi2*Trho(m_q2) * sum1;

  // total
  return( J_chi );
}

// General framework
void VA_0_PiPiPiPi1Charged::Calc(const ATOOLS::Vec4D_Vector& moms, bool m_anti)
{
  Vec4C help(0.,0.,0.,0.);
  // get Lorentz structure
  p_lorenz->SetPrivates(&moms.front(), &p_i.front());
  help = (*p_lorenz)();
  Insert( help*m_global , 0);
}

DEFINE_CURRENT_GETTER(VA_0_PiPiPiPi1Charged,"VA_0_PiPiPiPi1Charged")

void ATOOLS::Getter<Current_Base,ME_Parameters,VA_0_PiPiPiPi1Charged>::
PrintInfo(std::ostream &st,const size_t width) const {
  st<<"Example: $ 0 \\rightarrow \\pi^0 \\pi^0 \\pi^0 \\pi^+ $ \n\n"
    <<"Order: 0 = $\\pi^\\pm$, 1 = $\\pi^0$, 2 = $\\pi^0$, 3 = $\\pi^0$ "
    <<"Available form factors: \n "
    <<"  \\begin{itemize} \n"
    <<"    \\item Kuehn-Santamaria (default) \n"
    <<"  \\end{itemize} \n"
    <<"Reference: https://sherpa.hepforge.org/olddokuwiki/data/media/publications/theses/diplom\\__laubrich.pdf \n"
    <<std::endl;
}
