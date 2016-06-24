#include "ATOOLS/Phys/Flavour.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_File.H"

using namespace HADRONS;
using namespace ATOOLS;
using namespace std;

// 3 particle phase space function lambda
double ResonanceFlavour::Lambda( double a, double b, double c )
{
  double L = sqr(a-b-c)-4.*b*c;
  if (L>0.) return L;
  return 0.;
}

double ResonanceFlavour::Sqrt_Lambda( double a, double b, double c)
{
  return sqrt(Lambda(a,b,c));
}

SimpleResonanceFlavour::SimpleResonanceFlavour( std::string name, double _mass,
                                                double _width ) 
  : m_name(name), m_mass(_mass), m_width(_width), m_mass2(sqr(_mass))
{
}
  
ResonanceFlavour::ResonanceFlavour( kf_code _kfc, double _mass, double _width, int _run, string path )
  : SimpleResonanceFlavour( Flavour(_kfc).IDName(), _mass, _width)
{
  m_kfc = _kfc;
  m_running = _run;
  m_path = path;
  m_body = 2;
  p_hist = NULL;
  m_G_at_m2 = 1.;
}

void ResonanceFlavour::InitialiseThreeBodyResonance( ResonanceFlavour &res1 )
{
  InitialiseThreeBodyResonance(res1,res1,0.);
}

void ResonanceFlavour::InitialiseThreeBodyResonance( ResonanceFlavour &res1, ResonanceFlavour &res2, double beta )
{
  m_body=3;
  if( m_running ) {
    p_hist = CreateGHistogram( res1, res2, beta, kf_pi_plus );
    m_G_at_m2 = GetValueOfG(m_mass2);
  }
}

Histogram * ResonanceFlavour::CreateGHistogram( ResonanceFlavour res1, ResonanceFlavour res2, double beta, kf_code out )
{
  // create file name
  char fn[512];
  sprintf(fn, "GQ2_Mres=%.3f_Gres=%.3f_MresP=%.3f_GresP=%.3f_beta=%.3f_Mout=%.3f.dat",
      res1.Mass(), res1.Width(), res2.Mass(), res2.Width(), beta, Flavour(out).HadMass() );

  // look if file already exists
  My_In_File f("",m_path+"PhaseSpaceFunctions/"+fn);
  if( !f.Open() ) {                            // if file does not exist
    // create histogram (i.e. table of values)
    msg_Out()<<"Create necessary phase space function for chosen parameters.\n"
             <<"This may take some time. Please wait..."<<endl;
    msg_Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
             <<"     Create G(q2) in "<<fn<<"."<<endl;
    double low (0.), up (3.2);
    int nbins  (50);
    double step = (up-low)/nbins;
    Histogram* myHist = new Histogram( 0, low, up+step, nbins+1 );
    double q2  (low), phi (0.);
    while( q2<=up+step ) {
      phi = IntegralG(q2,res1,res2,beta,out);     // get G value
      myHist->Insert( q2, phi );        // insert into histogram
      q2 += step;       
    }
    myHist->Output(m_path+"PhaseSpaceFunctions/"+fn);         // write into file
    return myHist;
  }
  else {
    // read table and create histogram
    msg_Tracking()<<"HADRONS::Tau_Three_Pseudo::KS::CreateGHistogram : \n"
             <<"     Read G(q2) from "<<fn<<"."<<endl;
    std::string found_file_name = f.Path()+f.File();
    f.Close();
    return new Histogram( found_file_name );
  }
}

double ResonanceFlavour::GetValueOfG( double s )
{
  double val(0);
  p_hist->Extrapolate( s+p_hist->BinSize(), &val, 0 );
                            // shift +BinSize to get right value
  return val;
}

double ResonanceFlavour::IntegralG( double Q2, ResonanceFlavour res1, ResonanceFlavour res2, double beta, kf_code out )
{
  if (Q2==0.0) return 0.0;
  int Ns=500, Nt=500;                   // number of subintervals
  double sum (0.);
  double mpi2 = sqr( Flavour(out).HadMass() );
  double s_max = Q2 + mpi2 - 2.*sqrt(Q2*mpi2);
  double s_min = 4.*mpi2;
  double ds = (s_max-s_min)/Ns;
  double t_max (0.);
  double t_min (0.);
  double dt (0.);
  double s (s_min), t, u;
  double V12, V22, V1V2;
  
  while ( s<s_max ) {
    t_max = ( sqr(Q2-mpi2) - sqr( Sqrt_Lambda(Q2,s,mpi2)-Sqrt_Lambda(s,mpi2,mpi2) ) )/(4.*s);
    t_min = ( sqr(Q2-mpi2) - sqr( Sqrt_Lambda(Q2,s,mpi2)+Sqrt_Lambda(s,mpi2,mpi2) ) )/(4.*s);
    dt = (t_max-t_min)/Nt;
    t = t_min;
    while ( t<t_max ) {
      u    = Q2 - s - t + 3.*mpi2;
      V12  = 4.*mpi2 - s - sqr(u-t)/(4.*Q2);
      V22  = 4.*mpi2 - t - sqr(u-s)/(4.*Q2);
      V1V2 = (u-s-t+4.*mpi2)/2. - (u-t)*(u-s)/(4.*Q2);
      Complex BW_s = (res1.BreitWigner(s)+beta*res2.BreitWigner(s))/(1.+beta);
      Complex BW_t = (res1.BreitWigner(t)+beta*res2.BreitWigner(t))/(1.+beta);
      sum += ( V12 * norm(BW_s) + V22*norm(BW_t) + 2.*V1V2*real( BW_s*conj(BW_t)) )*ds*dt;
      t += dt;
    }
    s += ds;
  }
  return sum/Q2;            // global factors do not matter as they cancel anyway
} 
 
 
double ResonanceFlavour::OffShellMassWidth( double s )
{
  if( !m_running ) return m_mass*m_width;
  if( m_body == 3 ) return ThreeBodyResonanceMassWidth(s);
  // 2-body resonances
  switch( m_kfc ) {
    case kf_f_0_600:
      return TwoBodyResonanceMassWidth_12(s,Flavour(kf_pi_plus).HadMass());
    case kf_rho_770:
    case kf_rho_1450:
    case kf_rho_1700:
    case kf_f_0_980:
      return TwoBodyResonanceMassWidth(s,Flavour(kf_pi_plus).HadMass());
    case kf_rho_770_plus:
    case kf_rho_1450_plus:
    case kf_rho_1700_plus:
      return TwoBodyResonanceMassWidth(s,Flavour(kf_pi).HadMass(), Flavour(kf_pi_plus).HadMass());
    case kf_K_star_892:
    case kf_K_star_1410:
    case kf_K_star_1680:
    case kf_K_0_star_1430:
      return TwoBodyResonanceMassWidth(s,Flavour(kf_pi_plus).HadMass(), Flavour(kf_K_plus).HadMass());
    case kf_K_star_892_plus:
    case kf_K_star_1410_plus:
    case kf_K_star_1680_plus:
    case kf_K_0_star_1430_plus:
      return TwoBodyResonanceMassWidth(s,Flavour(kf_pi).HadMass(), Flavour(kf_K_plus).HadMass());
  }
  msg_Error()<<"WARNING in ResonanceFlavour::OffShellMassWidth(double s) : "<<endl
    <<"     OffShellWidth of "<<m_name<<" hasn't been implemented yet."<<endl
    <<"     Will ignore and return on-shell mass width"<<endl;
  cout<<om::red<<"I changed my mind. Will abort."<<om::reset<<endl; abort();
  return m_mass*m_width;
}

double ResonanceFlavour::ThreeBodyResonanceMassWidth( double s )
{
  if( p_hist ) {
    return m_mass*m_width*GetValueOfG(s)/m_G_at_m2;
  }
  else {
    msg_Error()<<"ERROR in ResonanceFlavour::ThreeBodyResonanceMassWidth() : "<<endl
      <<"     No histogram for "<<m_name<<" has been constructed. Use Method InitialiseThreeBodyResonance() first."<<endl
      <<"     Don't know what to do. Will abort."<<endl;
    abort();
  }
  return 0.;
}

double ResonanceFlavour::TwoBodyResonanceMassWidth_12( double s, double m )
{
  double ms=sqr(m);
  if (s>4.*ms && m_mass2>4.*ms)
    return( m_width*m_mass2/sqrt(s) * pow( (s-4.*ms)/(m_mass2-4.*ms), 0.5 ) );
  return 0.;	
}
 
double ResonanceFlavour::TwoBodyResonanceMassWidth( double s, double m )
{
  double ms=sqr(m);
  if (s>4.*ms && m_mass2>4.*ms)
    return( sqrt(s)*m_width*m_mass2/s * pow( (s-4.*ms)/(m_mass2-4.*ms), 1.5 ) );
  return 0.;	
}
 
double ResonanceFlavour::TwoBodyResonanceMassWidth( double s, double m1, double m2 )
{
  double threshold = sqr(m1+m2);
  double ms1 (m1*m1), ms2 (m2*m2);
  if (m_mass2>threshold && s>threshold)
    return( sqrt(s)*m_width*m_mass2/s * pow( m_mass2/s*Lambda(s,ms1,ms2)/Lambda(m_mass2,ms1,ms2), 1.5 ) );
  return 0;
}

Complex ResonanceFlavour::BreitWigner( double s )
{
  double MG;
  if( m_running ) {
    MG = OffShellMassWidth(s);
    return Complex(m_mass2,0.)/Complex( m_mass2-s, -1.*MG );
  }
  else {
    MG = m_mass * m_width;
    return Complex(m_mass2, -1.*MG)/Complex( m_mass2-s, -1.*MG );
  }
}
 
Complex ResonanceFlavour::BreitWignerAlt( double s )
{
  double MG = m_mass * m_width;
  return Complex(1.,0.)/Complex( m_mass2-s, -1.*MG );
}
