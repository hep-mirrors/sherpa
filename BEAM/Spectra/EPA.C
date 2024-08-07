#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/EPA_FF.H"
#include "BEAM/Spectra/EPA_Spectra_Plotter.H"
#include "ATOOLS/Math/Special_Functions.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

EPA::EPA(const Flavour beam, const double energy, const double pol, const int dir) :
  Beam_Base(beamspectrum::EPA, beam, energy, pol, dir),
  m_type(EPA_ff_type::point), p_ff(nullptr), 
  m_mass(m_beam.Mass(true)), m_charge(m_beam.Charge()), m_gamma(m_energy/m_mass)
{
  Initialise();
  m_Nbunches   = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_photon);
  m_bunches[1] = m_beam;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0.,0.,0.,0.);
  m_on         = true;
  exit(1);
}

void EPA::FixPosition() {
  /*
  // This is a bit of a poor-man's choice for a point-like source,
  // with a minmimal distance m_minR ... we would need some notion of
  // off'shellness here ...
  double ratio = m_maxR/m_minR, logratio = log(ratio), R, phi;
  if (ran->Get()< logratio/(0.5+logratio)) {
    R = m_minR * pow(ratio,ran->Get());
  }
  else {
    R = m_minR * sqrt(ran->Get());
  }
  phi = 2.*M_PI*ran->Get();
  m_position = R * Vec4D(0., cos(phi), sin(phi), 0.);
  */
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i) {
  if (i==0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab-out;
  }
}


double EPA::phi(double x, double qq) const {
  if (m_beam.Kfcode() == kf_p_plus) {
    const double a = 7.16;
    const double b = -3.96;
    const double c = .028;
    double y, qq1, f;
    qq1 = 1 + qq;
    y = x * x / (1 - x);
    f = (1 + a * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) +
                       1 / (3 * qq1 * qq1 * qq1));
    f += (1 - b) * y / (4 * qq * qq1 * qq1 * qq1);
    f += c * (1 + y / 4) *
         (log((qq1 - b) / qq1) + b / qq1 + b * b / (2 * qq1 * qq1) +
          b * b * b / (3 * qq1 * qq1 * qq1));
    return f;
  }
  if (m_beam.IsIon()) {
    // x := omega / omega0 is assumed in the following code!
    // ensure whether calls of phi for ions are done correctly
    // x_omega=x*E/omega0=x*E*R/gamma
    double f = 0.;
    // needed for gaussian shaped nucleus
    const double q0 = 0.06;
    const int atomicNumber = m_beam.GetAtomicNumber();
    const double radius = 1.2 / .197 * pow(atomicNumber, 1. / 3.);
    //CosInt Ci;
    // do form factor dependent calculation
    switch (m_type) {
    case EPA_ff_type::point:
      f = log(1. + (1. / (x * x))) / 2. + 1. / (1. + (1. / (x * x))) / 2. -
          1. / 2.;
      break;
    case EPA_ff_type::Gauss: 
      f = (1. + x * x / (q0 * q0 * radius * radius));
      f *= ExpIntegral(1, x * x / (q0 * q0 * radius * radius));
      f -= exp(-x * x / (q0 * q0 * radius * radius));
      f /= 2.;
      break;
    case EPA_ff_type::hcs:
      f += 3. / (16. * pow(x, 6.));
      f += 3. / (8. * pow(x, 4.));
      f -= cos(2. * x) * 3. / (16 * pow(x, 6.)) +
           cos(2. * x) * 7. / (40. * x * x);
      f -= cos(2. * x) * 1. / 20.;
      f -= sin(2. * x) * 3. / (8 * pow(x, 5.)) +
           sin(2. * x) * 1. / (10. * x * x * x);
      f += sin(2. * x) * 9. / (20. * x) - sin(2. * x) * x / 10.;
      //f -= Ci.GetCosInt(2. * x) * (1. + pow(x, 5.) / 5.); // integral-cosine
      break;
    case EPA_ff_type::smooth_hcs:
      if (x < 0.003) { // make n(x) smooth at low x
        f = 1.83698 * pow(x, -0.00652101) * M_PI * m_energy;
        // f=1.36549*pow(x,-0.059967)*M_PI*m_energy*atomicNumber;
        //  prefactor*c*x^a with c and a from a fit to x_omega*n(x_omega)
        f /= (2 * m_aqed * m_charge * m_charge * radius * m_beam.Mass());
      } else if (x > 1.33086) { // cut off oscillating parts at high x
        f = 0.;
      } else { // normal homogenously charged sphere
        f += 3. / (16. * pow(x, 6.));
        f += 3. / (8. * pow(x, 4.));
        f -= cos(2. * x) * 3. / (16 * pow(x, 6.)) +
             cos(2. * x) * 7. / (40. * x * x);
        f -= cos(2. * x) * 1. / 20.;
        f -= sin(2. * x) * 3. / (8 * pow(x, 5.)) +
             sin(2. * x) * 1. / (10. * x * x * x);
        f += sin(2. * x) * 9. / (20. * x) - sin(2. * x) * x / 10.;
        //f -= Ci.GetCosInt(2. * x) * (1. + pow(x, 5.) / 5.); // integral-cosine
      }
      break;
    default:
      THROW(fatal_error, "Unknown ion form factor chosen");
    }
    return f;
  }
  return 0.;
}


void EPA::Initialise() {
  Settings &s = Settings::GetMainSettings();
  RegisterDefaults();
  m_aqed      = s["EPA"]["AlphaQED"].Get<double>();
  m_pref      = sqr(m_charge)*m_aqed/M_PI;
  m_approx    = s["EPA"]["Approximation"].Get<size_t>();  
  m_analytic  = s["EPA"]["AnalyticFF"].Get<bool>();  
  m_q2max     = ExtractParameter(s,"Q2Max");
  m_q2min     = ExtractParameter(s,"Q2Min");
  m_theta_max = ExtractParameter(s,"ThetaMax");
  m_pt2max    = sqr(m_energy*m_theta_max);
  m_pt2min    = ExtractParameter(s,"PT2Min");
  m_xmin      = ExtractParameter(s,"xMin");
  m_xmax      = ExtractParameter(s,"xMax");
  m_bmin      = ExtractParameter(s,"bMin"); 
  m_bmax      = ExtractParameter(s,"bMax"); 
  m_nxbins    = s["EPA"]["xBins"].Get<int>();
  m_nbbins    = s["EPA"]["bBins"].Get<int>();

  //InitFormFactor(s);
  //InitTables();
  //WriteDebugFiles(s);
  EPA_Spectra_Plotter plotter(this,string("Spectra"));
  plotter(99);
}

void EPA::RegisterDefaults() const {
  Settings &s = Settings::GetMainSettings();
  s["EPA"]["Q2Max"].SetDefault(3.0);
  s["EPA"]["Q2Min"].SetDefault(-1.);
  s["EPA"]["xMax"].SetDefault(1.);
  s["EPA"]["xMin"].SetDefault(0.);
  s["EPA"]["xBins"].SetDefault(200);
  // impact parameters in fm.  need to make this more elegant 
  s["EPA"]["bMin"].SetDefault(1.e-3);
  s["EPA"]["bMax"].SetDefault(1.e4);
  s["EPA"]["bBins"].SetDefault(100);
  s["EPA"]["PT2Min"].SetDefault(0.0);
  s["EPA"]["Form_Factor"].SetDefault(m_beam.FormFactor());
  s["EPA"]["AlphaQED"].SetDefault(0.0072992701);
  s["EPA"]["ThetaMax"].SetDefault(0.3);
  s["EPA"]["Approximation"].SetDefault(1);
  s["EPA"]["AnalyticFF"].SetDefault(true);
  s["EPA"]["Debug"].SetDefault(false);
  s["EPA"]["Debug_Files"].SetDefault("EPA_debugOutput");
}

double EPA::ExtractParameter(Settings &s,const std::string & tag) {
  std::vector<double> parms = s["EPA"][tag].GetVector<double>();
  if (parms.size()!=1 && parms.size()!=2)
    THROW(fatal_error, "Specify either one or two values for 'EPA:"+tag+"'.  Will exit.");
  double parm = (m_dir > 0) ? parms.front() : parms.back(); 
  if (tag=="PTMin" && parm>1.0) {
    /* pt2min > 1 - according to approximation of
       'qmi' calculation in CalculateWeight */
    THROW(critical_error, "Too big p_T cut 'EPA:"+tag+"'.  Will exit.");
  }
  return parm;
}

void EPA::InitFormFactor(Settings &s) {
  std::vector<int> formfactors = s["EPA"]["Form_Factor"].GetVector<int>();
  if (formfactors.size()!=1 && formfactors.size()!=2)
    THROW(fatal_error,
          "Specify either one or two values for `EPA:Form_Factor'.");
  int formfactor = (m_dir > 0) ? formfactors.front() : formfactors.back();
  switch (formfactor) {
  case  0: 
    m_type     = EPA_ff_type::point;
    m_analytic = true;
    p_ff       = new EPA_Point(m_beam,m_energy);
    break;
  case  1:
    m_type     = EPA_ff_type::Gauss;
    p_ff       = new EPA_Gauss(m_beam,m_energy);
    break;
  case  2:
    m_type     = EPA_ff_type::dipole;
    p_ff       = new EPA_Dipole(m_beam,m_energy);
    break;
  case 11:
    m_type     = EPA_ff_type::hcs;
    THROW(fatal_error,
          "Form factor not yet implemented "+ToString(int(m_type)));
    break;
  case 12:
    m_type     = EPA_ff_type::smooth_hcs;
    THROW(fatal_error,
          "Form factor not yet implemented "+ToString(int(m_type)));
    break;
  case 13:
    m_type     = EPA_ff_type::WoodSaxon;
    THROW(fatal_error,
          "Form factor not yet implemented "+ToString(int(m_type)));
    break;
  default:
    THROW(fatal_error,
          "unspecified EPA form factor: "+ToString(formfactor));
  }
}

void EPA::WriteDebugFiles(Settings &s) {
  if (s["EPA"]["Debug"].Get<bool>()) {
    std::vector<std::string> files{
      s["EPA"]["Debug_Files"].GetVector<std::string>()};
    if (files.size() != 1 && files.size() != 2)
      THROW(fatal_error,
            "Specify either one or two values for `EPA:Debug_Files'.");
    std::string filename{(m_dir > 0) ? files.front() : files.back()};
    std::string num(m_dir > 0 ? "1" : "2");
    filename += num + ".log";
    this->selfTest(filename);
  }
}

void EPA::selfTest(const std::string& filename) {
  std::ofstream debugOutput;
  debugOutput.open(filename.c_str());

  debugOutput << "# EPA::selfTest() starting ..." << std::endl;

  // select output format
  debugOutput.setf(std::ios::scientific, std::ios::floatfield);
  debugOutput.precision(10);

  double x_omega = .1e-2;
  const int atomicNumber = m_beam.GetAtomicNumber();
  const double radius = 1.2 / .197 * pow(atomicNumber, 1. / 3.);
  double omega0, gamma;
  gamma = m_energy / m_beam.Mass();
  // gamma = m_energy * atomicNumber / m_beam.Mass();
  //  energy is defined as sqrt[s_NN], N=nucleon
  //  but recalculated already in the Beam_Spectra_Handler
  omega0 = gamma / radius;

  // write parameters
  debugOutput << "# Form Factor: " << int(m_type) << std::endl;
  debugOutput << "# A= " << atomicNumber << std::endl;
  debugOutput << "# R= " << radius << std::endl;
  debugOutput << "# E= " << m_energy << std::endl;
  debugOutput << "# Z= " << m_charge << std::endl;
  debugOutput << "# M_Ion=" << m_beam.Mass() << std::endl;
  debugOutput << "# gamma= " << gamma << std::endl;
  debugOutput << "# omega0= " << omega0 << std::endl;

  // write spectrum
  while (x_omega < 5) {
    x_omega *= 1.005;
    CalculateWeight(x_omega * omega0 / m_energy, 0); // m_weight = n(x)
    debugOutput << x_omega << "\t" << x_omega * m_weight / m_energy
                << std::endl;
  }

  debugOutput << "# EPA::selfTest() finished" << std::endl << std::endl;
  debugOutput.close();
  return;
}

void EPA::InitdN_by_dx() {
  double xmin = Max(m_xmin, 1.e-6), xmax = Min(m_xmax,1.-1.e-12);
  p_N_x       = new OneDim_Table(axis(m_nxbins,  xmin,  xmax,axis_mode::log));
  KperpIntegrand   ktint(p_ff);
  Gauss_Integrator gauss(&ktint);
  msg_Out()<<METHOD<<"(x in ["<<xmin<<", "<<xmax<<"], "
	   <<"type = "<<int(m_type)<<" for "<<m_beam<<").\n";
  axis xaxis  = p_N_x->GetAxis();
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    double x = xaxis.x(i);
    ktint.SetX(x);
    double q2min    = sqr(m_mass*x)/(1.-x), Q2min = q2min+sqr(m_mass*x);
    double q2max    = m_q2max,              Q2max = q2max+sqr(m_mass*x);
    // Integrate with Gauss-Chebyshev due to log structure of integrand
    double integral = gauss.Integrate(q2min,q2max,1.e-3,2);
    double n_x      = m_pref/x * integral;
    double check    = ( M_PI * m_pref/x *
			( log(Q2max/Q2min) -
			  sqr(m_mass*x)*(1./Q2min-1./Q2max) ) );
    double check2   = (*this)(x);
    //( M_PI * pref/x *
    //			((1.+sqr(1.-x))/2. * log(q2max/q2min)-
    //			 sqr(x*m_mass)*(1./q2min-1./q2max)) );
    msg_Out()<<"   n(x = "<<x<<", q^2 in ["<<q2min<<", "<<q2max<<"]) = "
	     <<n_x<<" (reduced = "<<integral<<") vs. "<<check<<", "
	     <<"ratio = "<<(n_x/check)<<", ratio2 = "<<(check2/check)<<"\n";
  }
  exit(1);
}

void EPA::InitTables() {
  InitdN_by_dx();
  double xmin = Max(m_xmin, 1.e-6), xmax = Min(m_xmax,1.);
  msg_Out()<<METHOD<<"(x in ["<<xmin<<", "<<xmax<<"], "
	   <<"b in ["<<m_bmin<<", "<<m_bmax<<"], "
	   <<"type = "<<ToString(int(m_type))<<" for "<<m_beam<<").\n";
  p_N_xb   = new TwoDim_Table(axis(m_nxbins,  xmin,  xmax,axis_mode::log),
			      axis(m_nbbins,m_bmin,m_bmax,axis_mode::log));
  //EPA_FF_Base * ff = new EPA_Gauss(m_beam,m_energy);
  EPA_FF_Base * ff = new EPA_Point(m_beam,m_energy);
  KperpIntegrand ktint(ff);
  ktint.SetMass(m_mass=Flavour(kf_p_plus).Mass());
  ktint.SetMode(1);
  Bessel_Integrator bessel(&ktint,1);
  axis xaxis = p_N_xb->GetAxis(0), yaxis = p_N_xb->GetAxis(1);
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    double x = xaxis.x(i);
    ktint.SetX(x);
    for (size_t j=0;j<yaxis.m_nbins;j++) {
      double bT       = yaxis.x(j); 
      ktint.SetBT(bT);
      double ktmin    = 0.;
      double integral = bessel(ktmin/bT);
      if (integral<1.e-8) integral = 0.;
      double n_xb     = m_pref/x * sqr(integral);
      p_N_xb->Fill(i,j,n_xb);
      msg_Out()<<"   - 1/(xm) K("<<std::setw(12)<<std::setprecision(6)<<bT<<"*"
	       <<std::setw(12)<<std::setprecision(6)<<(x*m_mass)<<") = 1/(xm) * "
	       <<"K("<<std::setw(12)<<std::setprecision(6)<<(x*m_mass*bT)<<") = "
	       <<std::setw(12)<<std::setprecision(6)<<(integral/(x*m_mass))<<" vs. "
	       <<std::setw(12)<<std::setprecision(6)<<SF.Kn(1,x*m_mass*bT)<<", "
	       <<"ratio = "<<std::setw(12)<<std::setprecision(6)
	       <<(integral/(x*m_mass*SF.Kn(1,x*m_mass*bT)))<<".\n";
    }
  }
  BIntegrand bint(p_N_xb);
  Gauss_Integrator gaussB(&bint);
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    double x     = xaxis.x(i), x2 = x*x; 
    double vol   = M_PI*(sqr(m_bmax)-sqr(m_bmin));
    double q2min = sqr(m_mass * x) / (1. - x);
    double q2max = m_q2max;
    double check = (2.*M_PI*m_pref/x *
		    ((1.+sqr(1.-x)*log(q2max/q2min))-
		     2.*sqr(x*m_mass)*(1./q2min-1./q2max))/4.);
    /*
    switch (m_type) {
    case EPA_ff_type::Gauss:
      break;
    case EPA_ff_type::dipole:
      check = (*this)(x);
      break;
    case EPA_ff_type::point:
    default:
      break;
    }
    */
    bint.SetX(x);
    double GeV_fm = rpa->hBar()*rpa->c()*1.e12;
    double bmin   = sqrt(1./q2max);  //(x*m_mass*m_beam.Radius()/GeV_fm); 
    double bmax   = sqrt(1./q2min); //sqrt(sqr(x*m_mass*m_beam.Radius()/GeV_fm)+1.); 
    double n_x    = gaussB.Integrate(bmin,bmax,1.e-3);
    msg_Out()<<METHOD<<"(x = "<<xaxis.x(i)<<", R = "<<m_beam.Radius()<<", "
	     <<i<<" of "<<xaxis.m_nbins<<"): "
	     <<"n(x) = "<<n_x<<" vs. form("<<int(m_type)<<") = "<<check<<", "
	     <<"ratio = "<<(n_x/check)<<".\n";    
  }
  exit(1);
}

double KperpIntegrand::operator()(double xk) {
  switch (m_mode) {
  case 1: {
    // This is the mode used for the integration
    // dk_T^2 k_T^2 F[k_T^2+(x m)^2]/[k_T^2+(x m)^2] J_1(b_T k_T)
    // with xk = b_T * k_T  dk_T = dxk/b_T
    double kT = xk/m_bT, kT2 = sqr(kT), kT2tilde = kT2+sqr(m_x*m_mass);
    return 1./m_bT * kT2/kT2tilde * (*p_ff)(kT2tilde);
  }
  case 0:
  default:
    break;
  }
  // This is the mode used for the integration
  // d^2k_T k_T^2 { F[k_T^2+(x m)^2]/[k_T^2+(x m)^2] }^2 =
  // pi dk_T^2 k_T^2 { F[k_T^2+(x m)^2]/[k_T^2+(x m)^2] }^2 =
  // with xk = k_T^2
  double kT2 = xk, kT2tilde = kT2+sqr(m_x*m_mass);
  return M_PI * kT2 * sqr((*p_ff)(kT2tilde)/kT2tilde);
}

void EPA::TestIntegration() {
  msg_Out()<<METHOD<<" checks for Bessel functions:\n";
  for (size_t i=0;i<21;i++) {
    double x = i<10 ? double(i)/10. : pow(10.,(i-10)/2.);
    msg_Out()<<"   BesselK_1("<<x<<") = "
	     <<std::cyl_bessel_k(1.,x)<<" & "
	     <<"BesselJ_1("<<x<<") = "<<SF.Jn(1.,x)<<" & "
	     <<"BesselK_1("<<x<<") = "<<SF.Kn(1.,x)<<"\n";
  }
  msg_Out()<<METHOD<<" checks for point form factor:\n";
  EPA_FF_Base * ff = new EPA_Point(m_beam,m_energy);
  for (size_t i=0;i<20;i++) {
    double Q2 = i<10 ? double(i)/10 : pow(double(i-9),4);
    msg_Out()<<"   FF("<<Q2<<") = "<<(*ff)(Q2)<<"\n";
  }
  KperpIntegrand ktint(ff);
  ktint.SetMass(m_beam.IsLepton() ? m_mass : Flavour(kf_p_plus).Mass(true));
  double ktmin = 0.;    //sqrt(q2min);
  double ktmax = 1.e6; //sqrt(q2max); 
  Gauss_Integrator gauss(&ktint);
  axis xaxis      = p_N_xb->GetAxis(0), yaxis = p_N_xb->GetAxis(1);
  for (size_t i=0;i<xaxis.m_nbins;i++) {
    for (size_t j=0;j<yaxis.m_nbins;j++) {
      double x = xaxis.x(i), bT = yaxis.x(j);
      ktint.SetX(x);
      ktint.SetBT(bT);
      double integral = 0.; //gauss.Integrate(ktmin,ktmax,1.e-3);
      for (size_t N=0;N<10000000;N++) {
	double kT = ran->Get()*(ktmax-ktmin);
	integral += ktint(kT)*(ktmax-ktmin);
      }
      integral /= 10000000;
      double n_xb     = m_pref/x * sqr(integral);
      double check    = m_pref/x * sqr(x*m_mass*SF.Kn(1,x*m_mass*bT));
      double naive    = m_pref/x * 1./sqr(bT);
      msg_Out()<<METHOD<<"(x = "<<xaxis.x(i)<<", b = "<<yaxis.x(j)<<"): "
	       <<"n(x,b) = "<<n_xb<<" vs. form("<<int(m_type)<<") = "<<check<<", "
	       <<"naive = "<<naive<<", ratio = "<<(naive/check)<<".\n";    
    }
  }
  exit(1);
}



/*
double EPA::CosInt::GetCosInt(double X) {
  if (X < 0.) THROW(fatal_error,"method called with negative X");
  ATOOLS::Gauss_Integrator integrator(this);
  return integrator.Integrate(X, 100000., 1.e-4, 1);
}

*/
