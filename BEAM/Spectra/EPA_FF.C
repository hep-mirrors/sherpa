#include "BEAM/Spectra/EPA_FF.H"

#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"

#include <cmath>
#include <iostream>
#include <numbers>
#include <string>
using namespace BEAM;
using namespace ATOOLS;
using string = std::string;

////////////////////////////////////////////////////////////////////////////////
//
// Base class for different "form factors" that will enter the EPA spectra.
// By and large they incorporate any information about internal structure.
//
// When dealing with the impact parameters b note
// - that we obtain the particle radii in fm and that we also use fm for any
//   potential external setting of impact parameters;
// - that we will give back impact parameters in fm as well; but
// - that the impact parameters b and radii R usually show up in conjunction
//   with photon energies, transverse momenta or particle masses, all in GeV.
//   The products enter functions such as exponentials or Bessel functions,
//   hence we have to divide by hbar*c = 0.197 GeV fm to convert distances into
//   units of 1/GeV.
////////////////////////////////////////////////////////////////////////////////

EPA_FF_Base::EPA_FF_Base(const ATOOLS::Flavour& beam, const int dir)
    ://////////////////////////////////////////////////////////////////////////////
     //
     // Initialisation of relevant beam parameters:
     // note that the particle radius is in fm and transformed into 1/GeV
     //
     //////////////////////////////////////////////////////////////////////////////
      m_beam(beam), m_mass(beam.Mass(true)), m_mass2(ATOOLS::sqr(m_mass)),
      m_R(beam.Radius() / rpa->hBar_c()), m_q2min(-1.),
      m_q2max(1.), m_pt2max(-1.),
      m_Zsquared(beam.IsIon() ? sqr(m_beam.GetAtomicNumber()) : 1.), m_b(0.),
      p_N_xb(nullptr)
{
  const auto& s    = Settings::GetMainSettings()["EPA"];
  size_t      b    = dir > 0 ? 0 : 1;
  m_q2max          = s["Q2Max"].GetTwoVector<double>()[b];
  m_q2min          = s["Q2Min"].GetTwoVector<double>()[b];
  m_nxbins         = s["xBins"].GetTwoVector<int>()[b];
  m_nbbins         = s["bBins"].GetTwoVector<int>()[b];
  m_xmin           = s["xMin"].GetTwoVector<double>()[b];
  m_xmax           = s["xMax"].GetTwoVector<double>()[b];
  m_bmin           = s["bMin"].GetTwoVector<double>()[b];
  m_b_pl_threshold = s["bThreshold"].GetTwoVector<double>()[b] * m_R;
  m_bmax           = s["bMax"].GetTwoVector<double>()[b];

  // Pre-calculate the normalisation of the  distribution \in [bmin, bmax]
  m_norm_distribution = 0.5 * std::log((ATOOLS::sqr(m_bmax) + 1.) / (ATOOLS::sqr(m_bmin) + 1.));

  if (m_bmin <= 0. || m_bmin > m_bmax)
    THROW(invalid_input, "Unphysical input for EPA impact parameter. ");
  if (m_xmin < 0. || m_xmin > m_xmax)
    THROW(invalid_input, "Unphysical input for EPA x-limits. ");
  if (m_q2min > m_q2max)
    THROW(invalid_input, "Unphysical input for EPA Q2-limits. ");
}

void EPA_FF_Base::FillTables()
{
  axis xaxis(m_nxbins, m_xmin, m_xmax, axis_mode::log);
  axis baxis(m_nbbins, m_bmin * m_R, std::min(m_b_pl_threshold, m_bmax) * m_R,
             axis_mode::linear);

  //////////////////////////////////////////////////////////////////////////////
  //
  // N(x,b) is given by the square of the Fourier transform of
  // kt^2/(kt^2+m^2x^2) F(kt^2+m^2x^2), which leads to a Bessel function.
  // We assume that the units in the b axis are in 1/GeV
  //
  //////////////////////////////////////////////////////////////////////////////
  msg_Out() << METHOD << " in " << xaxis.m_nbins << " * " << baxis.m_nbins
            << " bins:\n"
            << "   x in [" << xaxis.m_xmin << ", " << xaxis.m_xmax << "], "
            << "b in [" << baxis.m_xmin << ", " << baxis.m_xmax << "], "
            << "from R = " << m_R << " 1/GeV = " << (m_R * rpa->hBar_c())
            << " fm.\n";
  p_N_xb                   = std::make_unique<TwoDim_Table>(xaxis, baxis);
  N_xb_int*         kernel = new N_xb_int(this);
  Bessel_Integrator bessel(kernel, 1);
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    for (size_t j = 0; j < baxis.m_nbins; j++) {
      kernel->SetXB(xaxis.x(i), baxis.x(j));
      // Jacobian is d^2b = b db dphi and phi can be integrated out immediately
      double value = 2 * m_Zsquared * sqr(bessel()) / xaxis.x(i) * baxis.x(j);
      p_N_xb->Fill(i, j, value);
    }
  }
  delete kernel;
}

void EPA_FF_Base::OutputToCSV(const std::string& type)
{
  double              step_x(std::log(m_xmax / m_xmin) / double(m_nxbins));
  std::vector<double> xs(m_nxbins);
  xs[0] = m_xmin;
  for (int i = 1; i < m_nxbins; ++i) { xs[i] = xs[i - 1] * std::exp(step_x); }

  double              step_b(std::log(m_bmax / m_bmin) / double(m_nbbins));
  std::vector<double> bs(m_nbbins);
  bs[0] = m_bmin * m_R;
  for (int i = 1; i < m_nbbins; ++i) { bs[i] = bs[i - 1] * std::exp(step_b); }

  double              q2max(100.), q2min(1.e-5);
  int                 nq2steps(1000);
  double              step_q2(std::log(q2max / q2min) / double(nq2steps));
  std::vector<double> q2s(nq2steps);
  q2s[0] = q2min;
  for (int i = 1; i < nq2steps; ++i) {
    q2s[i] = q2s[i - 1] * std::exp(step_q2);
  }
  std::string prefix = m_beam.IDName() + "_" + type + "_";

  if (p_N_xb) {
    std::ofstream outfile_Nxb(prefix + "N_x_b.csv");
    outfile_Nxb << "x,b,N" << std::endl;
    for (auto& x : xs) {
      for (auto& b : bs) {
        outfile_Nxb << x << "," << b << "," << (*p_N_xb)(x, b) << std::endl;
      }
    }
    outfile_Nxb.close();
  } else {
    std::ofstream outfile_Nxb(prefix + "N_x_b.csv");
    outfile_Nxb << "x,b,N" << std::endl;
    for (auto& x : xs) { outfile_Nxb << x << ",0," << N(x) << std::endl; }
    outfile_Nxb.close();
  }

  std::ofstream outfile_FFq2(prefix + "FF_q2.csv");
  outfile_FFq2 << "q2,FF" << std::endl;
  for (auto& q2 : q2s) outfile_FFq2 << q2 << "," << FF(q2) << std::endl;
  outfile_FFq2.close();
}

////////////////////////////////////////////////////////////////////////////////
// Point-like form factors and related functions in different approximations,
// all based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_Point::EPA_Point(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  // for point-like particles (i.e. leptons) we use the "classical"
  // lepton radius given by 1/alpha lambda_l/(2 pi)
  // with the Compton wavelength lambda_l
  m_b = rpa->hBar_c() / beam.Mass(true) / (2. * M_PI / 137.);
}

double EPA_Point::N(const double& x)
{
  // Budnev et al., Phys. Rep. C15 (1974) 181, Eq. (6.17b)
  // This is in units of [1]
  // Maximal angle for the scattered electron
  // compare hep-ph/9610406 and hep-ph/9310350
  double q2min = Q2min(x);
  double q2max = Q2max(x);
  return 1. / 2 / x *
         ((1 + sqr(1 - x)) * log(q2max / q2min) -
          2 * sqr(m_mass * x) * (1 / q2min - 1 / q2max));
}

////////////////////////////////////////////////////////////////////////////////
// Point-like form factors and related functions in different approximations,
// all based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_PointApprox::EPA_PointApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  // for point-like particles (i.e. leptons) we use the "classical"
  // lepton radius given by 1/alpha lambda_l/(2 pi)
  // with the Compton wavelength lambda_l
  m_b = rpa->hBar_c() / beam.Mass(true) / (2. * M_PI / 137.);
}

double EPA_PointApprox::N(const double& x)
{
  // Budnev et al., Phys. Rep. C15 (1974) 181, Eq. (6.17b)
  // This is in units of [1]
  // Maximal angle for the scattered electron
  // compare hep-ph/9610406 and hep-ph/9310350
  double q2min = Q2min(x);
  double q2max = Q2max(x);
  // V.M. Budnev et al., Phys. Rep. C15(1974)181, first term in eq. 6.17b
  return 1. / 2 * (1 + sqr(1 - x)) / x * log(q2max / q2min);
}

////////////////////////////////////////////////////////////////////////////////
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181
// c.f. eq. D.7 and Table 8
////////////////////////////////////////////////////////////////////////////////

EPA_Proton::EPA_Proton(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_mu2         = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];
  if (m_beam.Kfcode() != kf_p_plus)
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
}

double EPA_Proton::N(const double& x)
{
  auto phi = [this](double x, double z) {
    double y   = x * x / (1. - x);
    double a   = (1. + m_mu2) / 4. + 4. * m_mass2 / m_Q02;
    double b   = 1. - 4. * m_mass2 / m_Q02;
    double c   = (m_mu2 - 1.) * std::pow(b, -4);
    double zp1 = 1. + z;

    double term1 = (1. + a * y) *
                   (-std::log(1. + 1. / z) + 1. / zp1 +
                    1. / 2. * std::pow(zp1, -2) + 1. / 3. * std::pow(zp1, -3));
    double term2 = (1. - b) * y / 4. / z * std::pow(zp1, -3);
    double term3 = c * (1. + y / 4.) *
                   (std::log((zp1 - b) / zp1) + b / zp1 +
                    std::pow(b, 2) / 2. * std::pow(zp1, -2) +
                    std::pow(b, 3) / 3. * std::pow(zp1, -3));

    return term1 + term2 + term3;
  };

  double q2min = Q2min(x);
  return (1. - x) / x * (phi(x, m_q2max / m_Q02) - phi(x, q2min / m_Q02));
}

////////////////////////////////////////////////////////////////////////////////
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181,
// c.f. eq. D.6 and Table 8, but here with the
// approximation Q2 -> 0 in the form factor
////////////////////////////////////////////////////////////////////////////////

EPA_ProtonApprox::EPA_ProtonApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_mu2         = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  if (m_beam.Kfcode() != kf_p_plus)
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
}

double EPA_ProtonApprox::N(const double& x)
{
  double q2min = Q2min(x);
  return (1. - x + m_mu2 * x * x / 2.) / x * log(m_q2max / q2min) -
         (1. - x) / x * (1. - q2min / m_q2max);
}

////////////////////////////////////////////////////////////////////////////////
// Gaussian form factors replacing the dipole form factors of Budnev et al.,
// Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_Gauss::EPA_Gauss(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_Q02(1.)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];
  FillTables();
}

////////////////////////////////////////////////////////////////////////////////
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181,
// c.f. eq. D.4 and Table 8, but here with the
// approximation Q2 -> 0 in the form factor
////////////////////////////////////////////////////////////////////////////////

EPA_HCS::EPA_HCS(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  FillTables();
}

////////////////////////////////////////////////////////////////////////////////
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_Dipole::EPA_Dipole(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];

  FillTables();
}

////////////////////////////////////////////////////////////////////////////////
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_DipoleApprox::EPA_DipoleApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];

  FillTables();
}

////////////////////////////////////////////////////////////////////////////////
// Form factor for ions in a Wood-Saxon potential.  A few comments:
//
// - Output for the density modifed to recover the atomic number A on
//   integration, and in units of fm^-3 or GeV^3. Internally we normalise it to
//   unity, because we will later multiply the square of the form factor with
//   the square of the charge Z.
// - Radius in 1/GeV, therefore "nucelar skin" m_d also in 1/GeV
////////////////////////////////////////////////////////////////////////////////

double EPA_WoodSaxon::IntegrateWithAdaptiveRange(
        std::function<double(double)> integrand, double initial_rmax, double tolerance)
{
  auto functor = new Lambda_Functor(&integrand);
  Gauss_Integrator gauss(functor);
  double rmin = 0., rmax = initial_rmax;
  double total_result = gauss.Integrate(rmin, rmax, tolerance);
  double segment_increment = 0.;

  // Adaptively extend the integration range until the tail contributes negligibly.
  do {
    rmin = rmax;
    rmax *= 2.;
    segment_increment = gauss.Integrate(rmin, rmax, tolerance);
    total_result += segment_increment;
  } while (std::abs(segment_increment) > tolerance * std::abs(total_result) + 1.e-12);

  delete functor;

  return total_result;
}

EPA_WoodSaxon::EPA_WoodSaxon(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_d(0.5 / rpa->hBar_c()), m_rho0(0.)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_d           = s["WoodSaxon_d"].GetTwoVector<double>()[b] / rpa->hBar_c();
  if (!m_beam.IsIon())
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());

  m_rho0 = CalculateDensity();
  InitFFTable(1.e-12, 1.e4);
  FillTables();
}

double EPA_WoodSaxon::CalculateDensity()
{
  auto rho_integrand_lambda = [this](double r) {
    return r * r / (1. + std::exp((r - m_R) / m_d));
  };

  double res = IntegrateWithAdaptiveRange(rho_integrand_lambda, m_R, 1.e-6);

  return 1.0 / res / 4. / M_PI;
}

void EPA_WoodSaxon::InitFFTable(const double& q2min, const double& q2max)
{
  p_FF_Q2 = std::make_unique<OneDim_Table>(axis(100000, q2min, q2max, axis_mode::log));

  for (size_t i = 0; i < p_FF_Q2->GetAxis().m_nbins; i++) {
    const double q = std::sqrt(p_FF_Q2->GetAxis().x(i));

    // Define the form factor integrand as a lambda inside the loop.
    auto ff_integrand = [this, q](double r) {
      if (q * r < 1.e-4) { // Numerically stable limit for q*r -> 0
        return r * r / (1. + std::exp((r - m_R) / m_d));
      }
      return std::sin(q * r) / q * r / (1. + std::exp((r - m_R) / m_d));
    };

    double ff_integral = IntegrateWithAdaptiveRange(ff_integrand, m_R, 1.e-6);

    double form_factor = 4. * M_PI * m_rho0 * ff_integral;
    p_FF_Q2->Fill(i, form_factor);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Class for the Ion FF in the Electric Dipole approximation using
// the analytical expression from 2207.03012, eq. 3.3
////////////////////////////////////////////////////////////////////////////////

EPA_IonApprox::EPA_IonApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  m_b_pl_threshold = m_bmin * m_R;
  p_N_xb.reset();
}

////////////////////////////////////////////////////////////////////////////////
// Dummy test class for checking the Bessel_Integrator
////////////////////////////////////////////////////////////////////////////////

EPA_Test::EPA_Test(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  FillTables();
}

void EPA_Test::FillTables()
{
  axis xaxis(100, 1.e-4, 1., axis_mode::log);
  axis baxis(100, 1.e-6, 1.e4, axis_mode::log);
  // Fill_Nxb_Table(xaxis, baxis);
}
