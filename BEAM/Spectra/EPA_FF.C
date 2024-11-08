#include "BEAM/Spectra/EPA_FF.H"

#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Special_Functions.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"

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
      m_R(beam.Radius() / rpa->hBar_c()), m_q2min(0.), m_q2max(1.),
      m_pt2max(-1.), p_N_xb(nullptr), p_Inv_xb(nullptr),
      m_approx(false)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_approx      = s["Approximation"].Get<bool>();
  m_q2max       = s["Q2Max"].GetTwoVector<double>()[b];
  m_q2min       = s["Q2Min"].GetTwoVector<double>()[b];
  m_nxbins      = s["xBins"].GetTwoVector<int>()[b];
  m_nbbins      = s["bBins"].GetTwoVector<int>()[b];
  m_xmin        = s["xMin"].GetTwoVector<double>()[b];
  m_xmax        = s["xMax"].GetTwoVector<double>()[b];
  m_bmin        = s["bMin"].GetTwoVector<double>()[b];
  m_bmax        = s["bMax"].GetTwoVector<double>()[b];
}

double EPA_FF_Base::SelectB(double x)
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // Internally, we assume b is in units of 1/GeV
  //
  //////////////////////////////////////////////////////////////////////////////
  return (*p_Inv_xb)(x, ran->Get());
}

void EPA_FF_Base::FillTables(const size_t& nx, const size_t& nb)
{
  axis xaxis(nx, m_xmin, m_xmax, axis_mode::log);
  axis baxis(nb, m_bmin * m_R, m_bmax * m_R, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
}

void EPA_FF_Base::Fill_Nxb_Table(axis& xaxis, axis& baxis)
{
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
  m_xmin                   = xaxis.m_xmin;
  m_xmax                   = xaxis.m_xmax;
  m_bmin                   = baxis.m_xmin;
  m_bmax                   = baxis.m_xmax;
  p_N_xb                   = new TwoDim_Table(xaxis, baxis);
  N_xb_int*         kernel = new N_xb_int(this);
  Bessel_Integrator bessel(kernel, 1);
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    for (size_t j = 0; j < baxis.m_nbins; j++) {
      kernel->SetXB(xaxis.x(i), baxis.x(j));
      double value = sqr(bessel());
      p_N_xb->Fill(i, j, value);
    }
  }
  delete kernel;
}

void EPA_FF_Base::Fill_Invxb_Table()
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // The distribution in impact parameter is given by
  //   d^2b N(x,b) = 2 Pi db b N(x,b)
  // so we have to solve for B (with ran being a random number):
  //   int_0^B db b N(x,b) = ran * int_0^infinity db b N(x,b) = ran * integral
  // We therefore have to produce a cumulative
  // histogram of 1/integral db b N(x,b)
  // and invert it to solve the equation above.
  //
  //////////////////////////////////////////////////////////////////////////////
  axis   xaxis = p_N_xb->GetAxis(0);
  size_t rbins = 5. * p_N_xb->GetAxis(1).m_nbins;
  axis   raxis = axis(rbins, 0., 1., axis_mode::linear);
  msg_Out() << METHOD << " in " << xaxis.m_nbins << " bins:\n"
            << "   x in [" << m_xmin << ", " << m_xmax << "], "
            << "add " << rbins << " bins for random numbers\n";
  p_Inv_xb = new TwoDim_Table(xaxis, raxis);
  OneDim_Table* cum_b;
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    double total = 0., b;
    cum_b        = p_N_xb->Cumulative(1., 1, i, total);
    cum_b->Rescale(1. / total);
    for (size_t j = 0; j < rbins; j++) {
      b = cum_b->Inverse(raxis.x(j));
      p_Inv_xb->Fill(i, j, b);
    }
    delete cum_b;
  }
}

void EPA_FF_Base::TestInvxb()
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // Testing the distribution in b N(x,b), comparing 1 M MC trials and
  // the correct numerical solution.
  //
  //////////////////////////////////////////////////////////////////////////////
  for (size_t i = 1; i < 4; i++) {
    Histogram*  histo = new Histogram(0, 0., 50., 50000);
    std::string name = string("Spectra/Bdist_x") + ToString(i) + string(".dat");
    double      x    = pow(0.1, i);
    msg_Out() << METHOD << "(x = " << x << ")\n";
    for (long int j = 0; j < 10000000; j++) {
      double b = (*p_Inv_xb)(x, ran->Get());
      histo->Insert(b, 1);
    }
    histo->Output(name);
    delete histo;
    Histogram*  histo1 = new Histogram(0, 0., 50., 50000);
    std::string name1 =
            string("Spectra/Bcalc_x") + ToString(i) + string(".dat");
    for (size_t j = 0; j < 50000; j++) {
      double b = (double(j) + 0.5) / 1000.;
      histo1->Insert(b, b * N(x, b));
    }
    histo1->Output(name1);
    delete histo1;
  }
}

void EPA_FF_Base::OutputToCSV(const std::string& type)
{
  double              step_x(std::log(m_xmax / m_xmin) / double(m_nxbins));
  std::vector<double> xs(m_nxbins);
  xs[0] = m_xmin;
  for (int i = 1; i < m_nxbins; ++i) { xs[i] = xs[i - 1] * std::exp(step_x); }

  double              step_b(std::log(m_bmax / m_bmin) / double(m_nbbins));
  std::vector<double> bs(m_nbbins);
  bs[0] = m_bmin;
  for (int i = 1; i < m_nbbins; ++i) { bs[i] = bs[i - 1] * std::exp(step_b); }

  double              q2max(1000.), q2min(1.e-3);
  int                 nq2steps(1000);
  double              step_q2(std::log(q2max / q2min) / double(nq2steps));
  std::vector<double> q2s(nq2steps);
  q2s[0] = q2min;
  for (int i = 1; i < nq2steps; ++i) {
    q2s[i] = q2s[i - 1] * std::exp(step_q2);
  }
  std::string prefix = m_beam.IDName() + "_" + type + "_";
  if (m_approx) prefix += "approx_";

  std::ofstream outfile_Nxq2(prefix + "N_x_q2.csv");
  outfile_Nxq2 << "x,q2,N" << std::endl;
  for (auto& x : xs) {
    for (auto& q2 : q2s) {
      outfile_Nxq2 << x << "," << q2 << "," << this->operator()(x, q2)
                   << std::endl;
    }
  }
  outfile_Nxq2.close();

  std::ofstream outfile_Nxb(prefix + "N_x_b.csv");
  outfile_Nxb << "x,b,N" << std::endl;
  for (auto& x : xs) {
    for (auto& b : bs) {
      outfile_Nxb << x << "," << b << "," << N(x, b) << std::endl;
    }
  }
  outfile_Nxb.close();

  std::ofstream outfile_FFq2(prefix + "FF_q2.csv");
  outfile_FFq2 << "q2,FF" << std::endl;
  for (auto& q2 : q2s) { outfile_FFq2 << q2 << "," << FF(q2) << std::endl; }
  outfile_FFq2.close();
}

double N_xb_int::operator()(double y)
{
  // TODO change this
  //////////////////////////////////////////////////////////////////////////////
  //
  // Integration argument y here is bT*qT as mandated by the Bessel function:
  // - argument of form factor Q^2 = qT^2+x^2m^2 with qT^2 = y^2/bT^2
  // - we assume that m_b is in units of 1/GeV, qT is in GeV, overall results
  //   are in GeV.
  //
  //////////////////////////////////////////////////////////////////////////////
  double qT = y / m_b , qT2 = sqr(qT), xm2 = sqr(m_x * p_ff->Mass());
  double res = (*p_ff)(m_x, (qT2 + xm2) / (1. - m_x)) / m_b;
  return res;
}

////////////////////////////////////////////////////////////////////////////////
//
// Point-like form factors and related functions in different approximations,
// all based on Budnev et al., Phys. Rep. C15 (1974) 181.
//
////////////////////////////////////////////////////////////////////////////////

EPA_Point::EPA_Point(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  // FillTables(m_nxbins, m_nbbins);
}

double EPA_Point::operator()(const double& x, const double& Q2)
{
  const double q2min = ATOOLS::sqr(m_mass * x) / (1 - x);
  if (Q2 < q2min) return 0.;
  double wt = (1. + sqr(1. - x)) / 2.;
  if (!m_approx) wt -= (1. - x) * q2min / Q2;
  return wt / x / Q2;
}

void EPA_Point::FillTables(const size_t& nx, const size_t& nb)
{
  /*axis xaxis(nx, 1.e-5, 1., axis_mode::log);
  axis baxis(nb, 1.e-3, 1.e6, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();*/
}

double EPA_Point::N(const double& x, const double& b)
{
  // Budnev et al., Phys. Rep. C15 (1974) 181, Eq. (6.17b)
  // This is in units of [1]
  // Maximal angle for the scattered electron
  // compare hep-ph/9610406 and hep-ph/9310350
  double q2min = Q2min(x);
  double q2max = Q2max(x);
  double f;
  if (!m_approx)
    f = 1. / 2 / x *
        ((1 + sqr(1 - x)) * log(q2max / q2min) -
         2 * sqr(m_mass * x) * (1 / q2min - 1 / q2max));
  else// V.M. Budnev et al., Phys. Rep. C15(1974)181, first term in eq. 6.17b
    f = 1. / 2 * (1 + sqr(1 - x)) / x * log(q2max / q2min);
  if (f < 0) f = 0.;
  return f;
}

double EPA_Point::SelectB(double x) { return 0.; }

////////////////////////////////////////////////////////////////////////////////
//
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181
// c.f. eq. D.4 and Table 8
//
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

  FillTables(m_nxbins, m_nbbins);
}

double EPA_Proton::operator()(const double& x, const double& Q2)
{
  const double q2min(sqr(m_mass * x) / (1. - x));
  const double C(m_mu2);
  const double D(4. * m_mass2 + Q2 * m_mu2 / (4. * m_mass2 + Q2));
  return (sqr(x) / 2. * C + (1. - x) * (1. - q2min / Q2) * D) * sqr(FF(Q2)) /
         x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181,
// c.f. eq. D.4 and Table 8, but here with the
// approximation Q2 -> 0 in the form factor
//
////////////////////////////////////////////////////////////////////////////////

EPA_ProtonApprox::EPA_ProtonApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_mu2         = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  if (m_beam.Kfcode() != kf_p_plus)
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());

  FillTables(m_nxbins, m_nbbins);
}

double EPA_ProtonApprox::operator()(const double& x, const double& Q2)
{
  const double q2min(sqr(m_mass * x) / (1. - x));
  const double C(m_mu2);
  const double D(4. * m_mass2 + Q2 * m_mu2 / (4. * m_mass2 + Q2));
  // FF = 1, hence omitted
  return (sqr(x) / 2. * C + (1. - x) * (1. - q2min / Q2) * D) / x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Gaussian form factors replacing the dipole form factors of Budnev et al.,
// Phys. Rep. C15 (1974) 181.
//
////////////////////////////////////////////////////////////////////////////////

EPA_Gauss::EPA_Gauss(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_Q02(1.),
      m_Zsquared(beam.Kfcode() == kf_p_plus ? 1 : sqr(m_beam.GetAtomicNumber()))
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];
  FillTables(m_nxbins, m_nbbins);
}

double EPA_Gauss::operator()(const double& x, const double& Q2)
{
  const double q2min(sqr(m_mass * x) / (1. - x));
  // approximation with C = 0
  const double D(m_Zsquared * sqr(FF(Q2)));
  return (1. - x) * (1. - q2min / Q2) * D / x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Proton form factor as defined in Budnev et al., Phys. Rep. C15 (1974) 181,
// c.f. eq. D.4 and Table 8, but here with the
// approximation Q2 -> 0 in the form factor
//
////////////////////////////////////////////////////////////////////////////////

EPA_HCS::EPA_HCS(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir),
      m_Zsquared(beam.Kfcode() == kf_p_plus ? 1 : sqr(m_beam.GetAtomicNumber()))
{
  FillTables(m_nxbins, m_nbbins);
}

double EPA_HCS::operator()(const double& x, const double& Q2)
{
  // Taken from Vidovic et al., PhysRevC.47.2308
  // form factor is 3*j_1(Q*R)/(Q*R) with j_1 the spherical bessel function
  // there the approximation with C = 0 is implicit
  const double q2min(sqr(m_mass * x) / (1. - x));
  const double QR(sqrt(Q2) * m_R);
  const double D(m_Zsquared * sqr(FF(QR)));
  return ((1. - x) * (1. - q2min / Q2) * D) / x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
//
////////////////////////////////////////////////////////////////////////////////

EPA_Dipole::EPA_Dipole(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir),
      m_Zsquared(beam.Kfcode() == kf_p_plus ? 1 : sqr(m_beam.GetAtomicNumber()))
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];

  FillTables(m_nxbins, m_nbbins);
}

double EPA_Dipole::operator()(const double& x, const double& Q2)
{
  const double q2min(sqr(m_mass * x) / (1. - x));
  const double D(m_Zsquared * sqr(FF(Q2)));
  return ((1. - x) * (1. - q2min / Q2) * D) / x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
//
////////////////////////////////////////////////////////////////////////////////

EPA_DipoleApprox::EPA_DipoleApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir),
      m_Zsquared(beam.Kfcode() == kf_p_plus ? 1 : sqr(m_beam.GetAtomicNumber()))
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];

  FillTables(m_nxbins, m_nbbins);
}

double EPA_DipoleApprox::operator()(const double& x, const double& Q2)
{
  // neglecting the Q2 dependence of the form factor and approximating
  // it with Q2 = 0
  const double q2min(sqr(m_mass * x) / (1. - x));
  const double D(m_Zsquared);
  return ((1. - x) * (1. - q2min / Q2) * D) / x / Q2;
}

////////////////////////////////////////////////////////////////////////////////
//
// Form factor for ions in a Wood-Saxon potential.  A few comments:
//
// - Output for the density modifed to recover the atomic number A on
//   integration, and in units of fm^-3 or GeV^3. Internally we normalise it to
//   unity, because we will later multiply the square of the form factor with
//   the square of the charge Z.
// - Radius in 1/GeV, therefore "nucelar skin" m_d also in 1/GeV
//
////////////////////////////////////////////////////////////////////////////////

EPA_WoodSaxon::EPA_WoodSaxon(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_d(0.5 / rpa->hBar_c()), p_FF_Q2(nullptr),
      p_N(nullptr), m_rho0(0.)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_d           = s["WoodSaxon_d"].GetTwoVector<double>()[b];
  if (!m_beam.IsIon())
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());

  FillTables(m_nxbins, m_nbbins);
}

void EPA_WoodSaxon::FillTables(const size_t& nx, const size_t& nb)
{
  m_rho0 = CalculateDensity();
  InitFFTable(1.e-12, 1.e4);
  InitNTable(1.e-10, 1.);
  double xmin = 1.e-9 / m_beam.GetMassNumber();
  double xmax = 1. / m_beam.GetMassNumber();
  double bmin = 1.e-3 * m_R, bmax = 1.e6 * m_R;
  axis   xaxis(nx, xmin, xmax, axis_mode::log);
  axis   baxis(nb, bmin, bmax, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
}

void EPA_WoodSaxon::InitFFTable(const double& q2min, const double& q2max)
{
  p_FF_Q2 = new OneDim_Table(axis(100000, q2min, q2max, axis_mode::log));
  WS_potential*    ws = new WS_potential(m_R, m_d);
  Gauss_Integrator gauss(ws);
  for (size_t i = 0; i < p_FF_Q2->GetAxis().m_nbins; i++) {
    ws->SetQ(sqrt(p_FF_Q2->GetAxis().x(i)));
    double rmin = 0., rmax = m_R;
    double res = m_rho0 * gauss.Integrate(rmin, rmax, 1.e-6, 0);
    double inc(0.);
    do {// TODO what's the logic here?
      rmin = rmax;
      rmax *= 2.;
      res += inc = m_rho0 * gauss.Integrate(rmin, rmax, 1.e-6, 0);
    } while (rmax < 4. * m_R || dabs(inc / res) > 1.e-6);
    p_FF_Q2->Fill(i, res);
  }
  delete ws;
}

void EPA_WoodSaxon::InitNTable(const double& xmin, const double& xmax)
{
  p_N = new OneDim_Table(axis(10000, xmin, xmax, axis_mode::log));
  N_argument*      nx = new N_argument(this);
  Gauss_Integrator gauss(nx);
  for (size_t i = 0; i < p_N->GetAxis().m_nbins; i++) {
    double x     = p_N->GetAxis().x(i);
    double q2min = Q2min(x);
    double q2max = Q2min(x) + m_pt2max;
    nx->SetX(x);
    double res = gauss.Integrate(q2min, q2max, 1.e-6, 0);
    if (!(i % 1000))
      msg_Out() << METHOD << "(" << q2min << ", " << q2max << "): "
                << "N(" << x << ") = " << res << "\n";
    p_N->Fill(i, res);
  }
}

double EPA_WoodSaxon::CalculateDensity()
{
  Rho_argument*    rho = new Rho_argument(m_R, m_d);
  Gauss_Integrator gauss(rho);
  double           rmin = 0., rmax = m_R;
  double           res = gauss.Integrate(rmin, rmax, 1.e-3), inc = 0.;
  do {// TODO what's the logic here?
    rmin = rmax;
    rmax *= 2.;
    res += inc = gauss.Integrate(rmin, rmax, 1.e-6);
  } while (inc > 1.e-6 * res);
  return 1. / res;
}

double EPA_WoodSaxon::operator()(const double& x, const double& Q2)
{
  return (1. - x) * (1. - Q2min(x) / Q2) * (*p_FF_Q2)(Q2);
}

////////////////////////////////////////////////////////////////////////////////
//
// Dummy test class for checking the Bessel_Integrator
//
////////////////////////////////////////////////////////////////////////////////

EPA_Test::EPA_Test(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  FillTables(m_nxbins, m_nbbins);
}

void EPA_Test::FillTables(const size_t& nx, const size_t& nb)
{
  axis xaxis(100, 1.e-4, 1., axis_mode::log);
  axis baxis(100, 1.e-6, 1.e4, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
  TestInvxb();
}

double EPA_Test::operator()(const double& x, const double& Q2)
{
  return 1. / x / sqrt((1. - x) * Q2);
}
