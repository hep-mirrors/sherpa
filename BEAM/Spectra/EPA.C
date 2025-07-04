#include "BEAM/Spectra/EPA.H"

#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;

EPA::EPA(const Flavour _beam, const double _energy, const double _pol,
         const int _dir) :
  Beam_Base(beamspectrum::EPA, _beam, _energy, _pol, _dir),
  m_mass(_beam.Mass(true)), m_charge(_beam.Charge()), m_gamma(_energy / m_mass),
  m_minR(Max(1.e-6,_beam.Radius())), m_maxR(1.e6)
{
  Settings &s = Settings::GetMainSettings();
  RegisterDefaults();
  m_Nbunches   = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_photon);
  m_bunches[1] = m_beam;
  // for point-like particles (i.e. leptons) we use the "classical"
  // lepton radius given by 1/alpha lambda_l/(2 pi)
  // with the Compton wavelength lambda_l
  if (m_bunches[1].IsLepton())
    m_minR = ( rpa->hBar_c()/m_bunches[1].Mass(true) / (2.*M_PI/137.) );
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0.,0.,0.,0.);
  m_on         = true;

  std::vector<double> q2Max{s["EPA"]["Q2Max"].GetVector<double>()};
  if (q2Max.size() != 1 && q2Max.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `EPA:Q2Max'.");
  m_q2Max = (_dir > 0) ? q2Max.front() : q2Max.back();

  std::vector<double> q2Min{s["EPA"]["Q2Min"].GetVector<double>()};
  if (q2Min.size() != 1 && q2Min.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `EPA:Q2Min'.");
  m_q2Min = (_dir > 0) ? q2Min.front() : q2Min.back();

  std::vector<double> xmin{s["EPA"]["xMin"].GetTwoVector<double>()};
  m_xmin = (_dir > 0) ? xmin.front() : xmin.back();

  std::vector<double> xmax{s["EPA"]["xMax"].GetTwoVector<double>()};
  m_xmax = (_dir > 0) ? xmax.front() : xmax.back();

  std::vector<double> pt_min{s["EPA"]["PTMin"].GetVector<double>()};
  if (pt_min.size() != 1 && pt_min.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `EPA:PTMin'.");
  m_pt_min = (_dir > 0) ? pt_min.front() : pt_min.back();

  m_aqed = s["EPA"]["AlphaQED"].Get<double>();

  m_theta_max = s["EPA"]["ThetaMax"].Get<double>();

  m_lo_epa = s["EPA"]["Use_old_WW"].Get<bool>();

  std::vector<int> formfactors{s["EPA"]["Form_Factor"].GetVector<int>()};
  if (formfactors.size() != 1 && formfactors.size() != 2)
    THROW(fatal_error,
          "Specify either one or two values for `EPA:Form_Factor'.");
  m_formfactor = (_dir > 0) ? formfactors.front() : formfactors.back();

  if (m_pt_min > 1.0) {
    /* pt_min > 1 - according to approximation of
       'qmi' calculation in CalculateWeight */
    THROW(critical_error, "Too big p_T cut ( " + ToString(m_pt_min) + ")");
  }
  if (s["EPA"]["Debug"].Get<bool>()) {
    std::vector<std::string> files{
        s["EPA"]["Debug_Files"].GetVector<std::string>()};
    if (files.size() != 1 && files.size() != 2)
      THROW(fatal_error,
            "Specify either one or two values for `EPA:Debug_Files'.");
    std::string filename{(_dir > 0) ? files.front() : files.back()};
    std::string num(_dir > 0 ? "1" : "2");
    filename += num + ".log";
    this->selfTest(filename);
  }
}

void EPA::RegisterDefaults() const {
  Settings &s = Settings::GetMainSettings();
  s["EPA"]["Q2Max"].SetDefault(3.0);
  s["EPA"]["Q2Min"].SetDefault(-1.);
  s["EPA"]["xMax"].SetDefault(1.);
  s["EPA"]["xMin"].SetDefault(0.);
  s["EPA"]["PTMin"].SetDefault(0.0);
  s["EPA"]["Form_Factor"].SetDefault(m_beam.FormFactor());
  s["EPA"]["AlphaQED"].SetDefault(0.0072992701);
  s["EPA"]["ThetaMax"].SetDefault(0.3);
  s["EPA"]["Use_old_WW"].SetDefault(false);
  s["EPA"]["Debug"].SetDefault(false);
  s["EPA"]["Debug_Files"].SetDefault("EPA_debugOutput");
}

void EPA::FixPosition() {
  // This is a bit of a poor-man's choice for a point-like source,
  // with a minmimal distance m_minR, either given by the actual physical
  // radius or the classical charge radius of a point-like particle
  // n(omega,b) = alpha/pi^2 omega/gamma [K1(omega b/gamma)]^2,
  // where omega = x E with E the energy of the source and gamma its
  // Lorentz-boost E/m.
  //
  // All distances here given in fm and masses in GeV
  /*
  double ratio = m_maxR/m_minR, logratio = log(ratio), R, weight, k1, k10;
  do {
    if (ran->Get()< logratio/(0.5+logratio)) {
      R      = m_minR * pow(ratio,ran->Get());
      weight = (m_vecouts[0][0]*R/rpa->hBar_c());
    }
    else {
      R      = m_minR * sqrt(ran->Get());
      weight = m_vecouts[0][0] * sqr(m_minR)/R / rpa->hBar_c();
    }
    weight *= ( sqr(bessel_k1(R*m_vecouts[0][0]/m_gamma/rpa->hBar_c()) /
		    bessel_k1(Min(m_minR,R)*Min(1.,m_vecouts[0][0])/
			      m_gamma/rpa->hBar_c())) );
    if (false && weight>1.) {
      msg_Out()<<METHOD<<"(omega = "<<m_vecouts[0][0]<<", "
	       <<"R = "<<R<<" vs Rmin = "<<m_minR<<"):\n"
	       <<"*   K1("<<(R*m_vecouts[0][0]/m_gamma/rpa->hBar_c())<<") = "
	       <<(k1=bessel_k1(R*m_vecouts[0][0]/m_gamma/rpa->hBar_c()))<<" / "
	       <<(k10=bessel_k1(Min(m_minR,R)*Min(1.,m_vecouts[0][0])/
				m_gamma/rpa->hBar_c()))
	       <<" = "<<(k1/k10)<<" --> weight = "<<weight<<"\n\n";
    }
  } while (weight<ran->Get());
  */
  double R   = ran->Get()*2.;
  double phi = 2.*M_PI*ran->Get();
  m_position = R * Vec4D(0., cos(phi), sin(phi), 0.);
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D &out, const size_t & i) {
  if (i==0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab-out;
  }
}

double EPA::CosInt::GetCosInt(double X) {
  if (X < 0.) THROW(fatal_error,"method called with negative X");
  ATOOLS::Gauss_Integrator integrator(this);
  return integrator.Integrate(X, 100000., 1.e-4, 1);
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
    CosInt Ci;
    // do form factor dependent calculation
    switch (m_formfactor) {
    case 0: // point-like form factor
      f = log(1. + (1. / (x * x))) / 2. + 1. / (1. + (1. / (x * x))) / 2. -
          1. / 2.;
      break;
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
    case 2: // gaussian shaped nucleus
      f = (1. + x * x / (q0 * q0 * radius * radius));
      f *= ExpIntegral(1, x * x / (q0 * q0 * radius * radius));
      f -= exp(-x * x / (q0 * q0 * radius * radius));
      f /= 2.;
      break;
    case 3: // homogeneously charged sphere (smooth function at low and high x)
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
    default:
      THROW(fatal_error, "Unknown ion form factor chosen");
    }
    return f;
  }
  return 0.;
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
  debugOutput << "# Form Factor: " << m_formfactor << std::endl;
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

bool EPA::CalculateWeight(double x, double q2) {
  // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
  //                           omega = E-E' - energy of emitted photon
  const double alpha = m_aqed;
  m_x = x;
  m_Q2 = q2;
  if (x > 1. - m_mass / 2. / m_energy || x > m_xmax || x < m_xmin) {
    m_weight = 0.0;
    return true;
  }
  if (m_beam.IsLepton() && !m_lo_epa) {
    // Maximal angle for the scattered electron
    // compare hep-ph/9610406 and hep-ph/9310350
    double q2min = sqr(m_mass * m_x) / (1 - m_x);
    double q2max = Min(q2min + sqr(m_energy) * (1 - m_x) * sqr(m_theta_max),m_q2Max);
    q2min = Max(sqr(m_mass * m_x) / (1 - m_x), m_q2Min);
    double f = alpha / M_PI / 2 / m_x *
               ((1 + sqr(1 - m_x)) * log(q2max / q2min) -
                2 * sqr(m_mass * m_x) * (1 / q2min - 1 / q2max));
    if (f < 0)
      f = 0.;
    m_weight = f;
    msg_Debugging() << METHOD << "(x = " << m_x << ", q^2 = " << q2
                    << ") = " << f << ", "
                    << "energy = " << m_energy << ", "
                    << "mass = " << m_mass << ".\n";
    return true;
  } else if (m_beam.IsLepton() && m_lo_epa) {
    // V.M. Budnev et al., Phys. Rep. C15(1974)181, first term in eq. 6.17b
    double q2min = sqr(m_mass * m_x) / (1 - m_x);
    double f =
        alpha / M_PI / 2 * (1 + sqr(1 - m_x)) / m_x * log(m_q2Max / q2min);
    if (f < 0)
      f = 0.;
    m_weight = f;
    msg_Debugging() << METHOD << "(x = " << m_x << ", q^2 = " << q2
                    << ") = " << f << ", "
                    << "energy = " << m_energy << ", "
                    << "mass = " << m_mass << ".\n";
    return true;
  } else if (m_beam.Kfcode() == kf_p_plus) {
    const double qz = 0.71;
    double f, qmi, qma;
    qma = m_q2Max / qz;
    // x = omega/E = (E-E')/E  ; E,E' - incoming and outgoing protons energy
    //                           omega = E-E' - energy of emitted photon
    qmi = m_mass * m_mass * x * x / (1 - x) / qz;
    qmi += m_pt_min * m_pt_min / (1 - x) / qz;

    f = alpha / M_PI * (phi(x, qma) - phi(x, qmi)) * (1 - x) / x;
    f *= m_charge * m_charge;
    if (f < 0)
      f = 0.;
    m_weight = f;
    return true;
  } else if (m_beam.IsIon()) { // n(x)
    const int atomicNumber = m_beam.GetAtomicNumber();
    const double radius = 1.2 / .197 * pow(atomicNumber, 1. / 3.);
    double f, omega0, gamma;
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
    f *= phi(m_x * m_energy / omega0, m_Q2); // phi(x_omega, m_Q2)
    // x_omega=m_x*m_energy/omega0

    m_weight = f;
    return true;
  }
  return false;
}

Beam_Base *EPA::Copy() { return new EPA(*this); }
