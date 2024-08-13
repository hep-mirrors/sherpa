#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/EPA_FF.H"
#include "BEAM/Spectra/EPA_Spectra_Plotter.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;
using namespace std;

EPA::EPA(const Flavour beam, const double energy,
	 const double pol, const int dir) :
  Beam_Base(beamspectrum::EPA, beam, energy, pol, dir),
  m_type(EPA_ff_type::point), p_ff(nullptr), 
  m_mass(m_beam.Mass(true)), m_charge(m_beam.Charge()),
  m_gamma(m_energy/m_mass), m_plotting(0)
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
}

bool EPA::CalculateWeight(double x, double q2) {
  m_x = x; m_q2 = q2;
  m_weight = (*this)(m_x);
  if (IsNan(m_weight)) msg_Out()<<"Boink! "<<METHOD<<"(x = "<<x<<") yields NaN.\n";
  return true;
}

void EPA::FixPosition() {
  double R = p_ff->SelectB(m_x), phi = 2.*M_PI*ran->Get();
  m_position = R * Vec4D(0., cos(phi), sin(phi), 0.);
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i) {
  if (i==0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab-out;
    m_x          = out[0]/m_lab[0];
    m_q2         = dabs(out.Abs2());
  }
}

void EPA::Initialise() {
  Settings &s = Settings::GetMainSettings();
  RegisterDefaults();
  m_aqed      = s["EPA"]["AlphaQED"].Get<double>();
  m_pref      = sqr(m_charge)*m_aqed/M_PI;
  m_approx    = s["EPA"]["Approximation"].Get<size_t>();  
  m_analytic  = s["EPA"]["AnalyticFF"].Get<bool>();  
  m_plotting  = s["EPA"]["PlotSpectra"].Get<bool>();  
  m_q2max     = ExtractParameter(s,"Q2Max");
  m_q2min     = ExtractParameter(s,"Q2Min");
  m_theta_max = ExtractParameter(s,"ThetaMax");
  m_pt2max    = sqr(m_energy*m_theta_max);
  m_pt2min    = ExtractParameter(s,"PT2Min");
  m_xmin      = ExtractParameter(s,"xMin");
  m_xmax      = ExtractParameter(s,"xMax");
  m_bmin      = ExtractParameter(s,"bMin"); 
  m_bmax      = ExtractParameter(s,"bMax");
  m_WSd       = ExtractParameter(s,"WoodSaxon_d");
  m_mu        = ExtractParameter(s,"MagneticMu");
  m_Lambda2   = ExtractParameter(s,"Lambda2");
  m_Q02       = ExtractParameter(s,"Q02");
  m_nxbins    = s["EPA"]["xBins"].Get<int>();
  m_nbbins    = s["EPA"]["bBins"].Get<int>();
  
  InitFormFactor(s);
  //WriteDebugFiles(s);
  if (m_plotting>0) {
    EPA_Spectra_Plotter plotter(this,string("Spectra"));
    plotter(m_plotting);
    Tests();
  }
}

void EPA::RegisterDefaults() const {
  Settings &s = Settings::GetMainSettings();
  s["EPA"]["Q2Max"].SetDefault(3.0);
  s["EPA"]["Q2Min"].SetDefault(-1.);
  s["EPA"]["xMax"].SetDefault(1.);
  s["EPA"]["xMin"].SetDefault(0.);
  s["EPA"]["xBins"].SetDefault(0);
  s["EPA"]["bMin"].SetDefault(0.);
  s["EPA"]["bMax"].SetDefault(1.e12);
  s["EPA"]["bBins"].SetDefault(0);
  s["EPA"]["PT2Min"].SetDefault(0.0);
  s["EPA"]["Form_Factor"].SetDefault(size_t(m_beam.IsIon() ?     EPA_ff_type::WoodSaxon :
					    m_beam.IsNucleon() ? EPA_ff_type::dipole :
					    m_beam.IsMeson() ?   EPA_ff_type::dipole :
					    EPA_ff_type::point) );
  s["EPA"]["MagneticMu"].SetDefault(m_beam.IsNucleon() ? 2.79 : 0.);
  s["EPA"]["Lambda2"].SetDefault(0.71);
  s["EPA"]["Q02"].SetDefault(2.*rpa->hBar_c()/m_beam.Radius());
  s["EPA"]["WoodSaxon_d"].SetDefault(0.5);
  s["EPA"]["AlphaQED"].SetDefault(0.0072992701);
  s["EPA"]["ThetaMax"].SetDefault(0.3);
  s["EPA"]["Approximation"].SetDefault(1);
  s["EPA"]["AnalyticFF"].SetDefault(true);
  s["EPA"]["PlotSpectra"].SetDefault(0);
  s["EPA"]["Debug"].SetDefault(false);
  s["EPA"]["Debug_Files"].SetDefault("EPA_debugOutput");
}

double EPA::ExtractParameter(Settings &s,const std::string & tag) {
  std::vector<double> parms = s["EPA"][tag].GetVector<double>();
  if (parms.size()!=1 && parms.size()!=2)
    THROW(fatal_error,
	  "Specify either one or two values for 'EPA:"+tag+"'.  Will exit.");
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
    m_type = EPA_ff_type::point;
    p_ff   = new EPA_Point(m_beam);
    break;
  case  1:
    m_type = EPA_ff_type::Gauss;
    p_ff   = new EPA_Gauss(m_beam);
    p_ff->SetParam("Q02",m_Q02);
    p_ff->SetParam("Mu2",sqr(m_mu));
    break;
  case  2:
    m_type = EPA_ff_type::dipole;
    p_ff   = new EPA_Dipole(m_beam);
    p_ff->SetParam("Lambda2",m_Lambda2);
    p_ff->SetParam("Mu2",sqr(m_mu));
    break;
  case 13:
    m_type = EPA_ff_type::WoodSaxon;
    p_ff   = new EPA_WoodSaxon(m_beam);
    p_ff->SetParam("WSd",m_WSd);
    break;
  default:
    THROW(fatal_error,
          "unspecified EPA form factor: "+ToString(formfactor));
  }
  p_ff->SetSwitch("approximation",m_approx);
  p_ff->SetSwitch("analytic",m_analytic);
  p_ff->SetQ2Range(m_q2min,m_q2max);
  p_ff->SetPT2Range(m_pt2min,m_pt2max);
  p_ff->FillTables(m_nxbins,m_nbbins);
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

void EPA::Tests() {
  EPA_FF_Base * ff;
  m_theta_max = M_PI/180.;
  m_q2max     = 1.e99;
  // Testing the electron spectrum
  m_beam = Flavour(kf_e);
  m_mass = m_beam.Mass(true);
  ff     = new EPA_Point(m_beam);
  ff->SetPT2Max(sqr(m_energy*m_theta_max));
  ff->SetSwitch("Approximation",2);
  ff->SetSwitch("Analytic",1);
  for (size_t i=0;i<100;i++) {
    double x = double(i)/1000;
    CalculateWeight(x,0.);
    msg_Out()<<"N_point("<<x<<") = "<<OldWeight(x,0.)<<" vs. "<<m_weight<<"\n";
  }
  delete ff;
  // Testing the proton spectrum
  m_q2max = 16.;
  m_beam  = Flavour(kf_p_plus);
  m_mass  = m_beam.Mass(true);
  m_type  = EPA_ff_type::dipole;
  ff      = new EPA_Dipole(m_beam);
  ff->SetQ2Max(m_q2max);
  ff->SetSwitch("Approximation",1);
  ff->SetSwitch("Analytic",1);
  for (size_t i=0;i<100;i++) {
    double x = double(i)/1000;
    CalculateWeight(x,0.);
    msg_Out()<<"N_dipole("<<x<<") = "<<OldWeight(x,0.)<<" vs. "<<m_weight<<"\n";
  }
  delete ff;
  m_beam  = Flavour(kf_lead208);
  m_mass  = m_beam.Mass(true);
  m_type  = EPA_ff_type::Gauss;
  double xmin = 1.e-5/m_beam.GetAtomicNumber(), xmax = 1./m_beam.GetAtomicNumber();
  double bmin = 1.e-4*m_beam.Radius(), bmax = 1.e6*m_beam.Radius(); 
  ff      = new EPA_Gauss(m_beam);
  ff->SetQ2Max(m_q2max);
  ff->FillTables();
  ff->SetSwitch("Approximation",1);
  for (size_t i=0;i<100;i++) {
    double x = double(i)/1000/m_beam.GetAtomicNumber();
    msg_Out()<<"-----------------------------------------------------\n";
    CalculateWeight(x,0.);
    msg_Out()<<"N_gauss(x/A = "<<(x*m_beam.GetAtomicNumber())<<") = "
	     <<OldWeight(x,0.)<<" vs. "<<m_weight<<"\n";
  }
  delete ff;
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//
//  The stuff below can go when we are happy and merge into master.
//
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void EPA::selfTest(std::string filename) {
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

double EPA::OldWeight(double x, double q2) {
  // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
  //                           omega = E-E' - energy of emitted photon
  const double alpha = m_aqed;
  double f = 0.;
  m_x = x;
  m_Q2 = q2;
  if (x > 1. - m_mass / 2. / m_energy) return 0.;
  if (m_beam.Kfcode() == kf_e && true) {
    // Maximal angle for the scattered electron
    // compare hep-ph/9610406 and hep-ph/9310350
    double q2min = sqr(m_mass * m_x) / (1 - m_x);
    double q2max = Min(q2min +
		       sqr(m_energy) * (1 - m_x) * sqr(m_theta_max),m_q2max);
    f = alpha / M_PI / 2 / m_x *
      ((1 + sqr(1 - m_x)) * log(q2max / q2min) -
       2 * sqr(m_mass * m_x) * (1 / q2min - 1 / q2max));
  } else if (m_beam.Kfcode() == kf_e && false) {
    // V.M. Budnev et al., Phys. Rep. C15(1974)181, first term in eq. 6.17b
    double q2min = sqr(m_mass * m_x) / (1 - m_x);
    f = alpha / M_PI / 2 * (1 + sqr(1 - m_x)) / m_x * log(m_q2max / q2min);
    msg_Debugging() << METHOD << "(x = " << m_x << ", q^2 = " << q2
                    << ") = " << f << ", "
                    << "energy = " << m_energy << ", "
                    << "mass = " << m_mass << ".\n";
    return 1;
  } else if (m_beam.Kfcode() == kf_p_plus) {
    const double qz = 0.71;
    double qmi, qma;
    qma = m_q2max / qz;
    // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
    //                           omega = E-E' - energy of emitted photon
    qmi = m_mass * m_mass * x * x / (1 - x) / qz;
    qmi += m_pt2min / (1 - x) / qz;
    
    f = alpha / M_PI * (phi(x, qma) - phi(x, qmi)) * (1 - x) / x;
    
    f *= m_charge * m_charge;
  } else if (m_beam.IsIon()) { // n(x)
    const int atomicNumber = m_beam.GetAtomicNumber();
    const double radius = 1.2 / .197 * pow(atomicNumber, 1. / 3.);
    double omega0, gamma;
    gamma = m_energy / m_beam.Mass();
    // gamma = m_energy * atomicNumber / m_beam.Mass();
    //  energy is defined as sqrt[s_NN], N=nucleon
    //  but recalculated already in the Beam_Spectra_Handler
    omega0 = gamma / radius;
    /*
      std::cout << "radius=" << radius << std::endl;
      std::cout << "omega0=" << omega0 << std::endl;
      std::cout << "gamma=" << gamma << std::endl;
    */
    /*
      f = 2 * alpha * m_charge * m_charge / M_PI / (m_x * omega0);
      f *= phi(m_x, m_Q2);
    */
    f = 2 * alpha * m_charge * m_charge / M_PI / m_x;
    // since CalculateWeight() is dn=N(x)*dx/x and not dn=N(omega)*domega/omega
    // f = 2 * alpha * m_charge * m_charge / M_PI / (m_x * m_energy);

    msg_Out()<<" x = "<<(m_x*m_energy/omega0)<<"\n";
    f *= phi(m_x * m_energy / omega0, m_Q2); // phi(x_omega, m_Q2)
  }
  if (f < 0) f = 0.;
  return f;
}



double EPA::phi(double x, double qq) {
  if (m_beam.Kfcode() == kf_p_plus) {
    const double a = 7.16;
    const double b = -3.96;
    const double c = .028;
    double y, qq1, f;
    double term1, term2, term3;
    qq1 = 1 + qq;
    y = x * x / (1 - x);
    f = term1 = (1 + a * y) * (-log(qq1 / qq) + 1 / qq1 + 1 / (2 * qq1 * qq1) +
			       1 / (3 * qq1 * qq1 * qq1));
    f -= term2 = (1 - b) * y / (4 * qq * qq1 * qq1 * qq1);
    f += term3 = c * (1 + y / 4) *
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
    switch (m_type) {
    // switch (2) {
    case EPA_ff_type::point: // point-like form factor
      f = log(1. + (1. / (x * x))) / 2. + 1. / (1. + (1. / (x * x))) / 2. -
          1. / 2.;
      break;
    case EPA_ff_type::Gauss: // gaussian shaped nucleus
      f = (1. + x * x / (q0 * q0 * radius * radius));
      f *= ExpIntegral(1, x * x / (q0 * q0 * radius * radius));
      f -= exp(-x * x / (q0 * q0 * radius * radius));
      f /= 2.;
      msg_Out()<<METHOD<<"("<<x<<", q0 = "<<q0<<", R = "<<radius<<") = "<<f<<"\n";
      break;
      /*
	case 1: // homogeneously charged sphere
	f += 3. / (16. * pow(x, 6.));
	f += 3. / (8. * pow(x, 4.));
	f -= cos(2. * x) * 3. / (16 * pow(x, 6.)) +
	cos(2. * x) * 7. / (40. * x * x);
	f -= cos(2. * x) * 1. / 20.;
	f -= sin(2. * x) * 3. / (8 * pow(x, 5.)) +
	sin(2. * x) * 1. / (10. * x * x * x);
	f += sin(2. * x) * 9. / (20. * x) - sin(2. * x) * x / 10.;
	f -= Ci.GetCosInt(2. * x) * (1. + pow(x, 5.) / 5.); // integral-cosine
	break;
	case 3: // homogeneously charged sphere
	// (smooth function at low and high x)
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
        f -= Ci.GetCosInt(2. * x) * (1. + pow(x, 5.) / 5.); // integral-cosine
	}
	break;
      */
    default:
      THROW(fatal_error, "Unknown ion form factor chosen");
    }
    return (double)f;
  }
  return 0.;
}

