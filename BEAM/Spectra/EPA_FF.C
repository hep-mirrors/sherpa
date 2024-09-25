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

ATOOLS::Special_Functions ATOOLS::SF;

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
      m_beam(beam), m_mass(m_beam.Mass(true)), m_mass2(ATOOLS::sqr(m_mass)),
      m_R(m_beam.Radius() / rpa->hBar_c()), m_q2min(0.), m_q2max(1.e99),
      m_pt2max(-1.), p_Nred_x(nullptr), p_N_xb(nullptr), p_Inv_xb(nullptr),
      m_approx(false)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_approx      = s["Approximation"].Get<bool>();
  m_q2max       = s["Q2Max"].GetTwoVector<double>()[b];
  m_q2min       = s["Q2Min"].GetTwoVector<double>()[b];
  m_nxbins      = s["xBins"].GetTwoVector<int>()[b];
  m_nbbins      = s["bBins"].GetTwoVector<int>()[b];
}

double EPA_FF_Base::SelectB(const double& x)
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // Internally, we assume b is in units of 1/GeV, we return it in fm here.
  //
  //////////////////////////////////////////////////////////////////////////////
  double xeff = m_xmin;
  if (xeff < m_xmin) xeff = m_xmin * (1. + 1.e-6);
  if (xeff > m_xmax) xeff = m_xmax * (1. - 1.e-6);
  double b = (*p_Inv_xb)(xeff, ran->Get());
  return b * rpa->hBar_c();
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
      msg_Out() << METHOD << ": filling bin (" << i << ", " << j << ")\n";
      kernel->SetXB(xaxis.x(i), baxis.x(j));
      double value = sqr(bessel()) / M_PI;
      p_N_xb->Fill(i, j, value);
    }
  }
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

double N_xb_int::operator()(double y)
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // Integration argument y here is bT*qT as mandated by the Bessel function:
  // - argument of form factor Q^2 = qT^2+x^2m^2 with qT^2 = y^2/bT^2
  // - we assume that m_b is in units of 1/GeV, qT is in GeV, overall results
  //   are in GeV.
  //
  //////////////////////////////////////////////////////////////////////////////
  double qT = y / m_b, qT2 = sqr(qT), xm2 = sqr(m_x * p_ff->Mass());
  double res = (1. / m_b * qT2 / (qT2 + xm2) * (*p_ff)(m_x, qT2 + xm2));
  return res;
}

void EPA_FF_Base::Fill_Nredx_Table()
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // The b axis is defined by the particle radius taken in 1/GeV
  //
  //////////////////////////////////////////////////////////////////////////////
  axis xaxis              = p_N_xb->GetAxis(0);
  p_Nred_x                = new OneDim_Table(xaxis);
  Nred_x_int*      kernel = new Nred_x_int(this);
  Gauss_Integrator gauss(kernel);
  double           bmin = m_R, bmax = p_N_xb->GetAxis(1).m_xmax;
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    kernel->SetX(xaxis.x(i));
    double value = gauss.Integrate(bmin, bmax, 1.e-5);
    p_Nred_x->Fill(i, value);
  }
}

double Nred_x_int::operator()(double b)
{
  //////////////////////////////////////////////////////////////////////////////
  //
  // Result is in 1/GeV, will yield 1/GeV^2 upon integration
  //
  //////////////////////////////////////////////////////////////////////////////
  return 2. * M_PI * b * (*p_ff->GetN_xb())(m_x, b);
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
  FillTables(m_nxbins, m_nbbins);
}

double EPA_Point::operator()(const double& x, const double& Q2)
{
  double wt = (1. + sqr(1. - x)) / 2.;
  if (!m_approx) wt -= (1. - x) * Q2min(x) / Q2;
  return wt;
}

void EPA_Point::FillTables(const size_t& nx, const size_t& nb)
{
  axis xaxis(nx, 1.e-9, 1., axis_mode::log);
  axis baxis(nb, 1.e-3, 1.e6, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
}

void EPA_Point::Fill_Nxb_Table(axis& xaxis, axis& baxis)
{
  m_xmin = xaxis.m_xmin;
  m_xmax = xaxis.m_xmax;
  m_bmin = baxis.m_xmin;
  m_bmax = baxis.m_xmax;
  p_N_xb = new TwoDim_Table(xaxis, baxis);
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    for (size_t j = 0; j < baxis.m_nbins; j++) {
      p_N_xb->Fill(i, j, N(xaxis.x(i), baxis.x(j)));
    }
  }
}

double EPA_Point::N(const double& x)
{
  // Budnev et al., Phys. Rep. C15 (1974) 181, Eq. (6.17b)
  // This is in units of [1]
  double q2min = Q2min(x), q2max = Q2max(x);
  if (x <= 0. || x >= 1. || q2max <= q2min) return 0.;
  double wt = log(q2max / q2min) * (1. + sqr(1. - x)) / 2.;
  if (!m_approx) wt -= sqr(m_mass * x) * (1. / q2min - 1. / q2max);
  return wt;
}

double EPA_Point::SelectB(const double& x)
{
  double bmin = 1.e-12, bmax = 1.e12;
  if (m_R >= 1.e-12) bmin = 1.e-6 * m_R, bmax = 1.e12 * m_R;
  double b = bmin, maxval = b * N(x, b);
  size_t trials = 10000;
  while (trials > 0) {
    b = bmin * pow(bmax / bmin, ran->Get());
    if (b * N(x, b) / maxval > ran->Get()) break;
    --trials;
  }
  return b * rpa->hBar_c();
}

double EPA_Point::N(const double& x, const double& b)
{
  // Result is in GeV^2, inherited from the mass in the argument.
  // Integration over the (2D) impact parameter plane (in units of 1/GeV^2)
  // will yield a result in units of [1].
  double arg = m_mass * x * SF.Kn(1, m_mass * x * b);
  return sqr(arg) / M_PI;
}

double EPA_Point::ReducedN(const double& x)
{
  // Result is in units of [1].
  double mxR = m_mass * x * m_R, mxR2 = sqr(mxR);
  double K0 = SF.Kn(0, mxR), K1 = SF.Kn(1, mxR);
  return mxR2 * (sqr(K0) - sqr(K1)) + 2. * mxR * K0 * K1;
}

////////////////////////////////////////////////////////////////////////////////
//
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
// In its present form, this form factor makes sense for nucleons only.
//
////////////////////////////////////////////////////////////////////////////////

EPA_Dipole::EPA_Dipole(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_Lambda2(0.71)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_mu2         = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  m_Lambda2     = s["Lambda2"].GetTwoVector<double>()[b];
  if (!m_beam.IsNucleon())
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
  m_mu2 = m_beam.Charge() != 0. ? 2.79 * 2.79 : 0.;

  // a, b, c coeffients from Budnev et al., Eq. (D.7)
  m_aDip = (1. + m_mu2) / 4. + 4. * m_mass2 / m_Lambda2;// should be  7.16
  m_bDip = 1. - 4. * m_mass2 / m_Lambda2;               // should be -3.96
  m_cDip = (m_mu2 - 1.) / pow(m_bDip, 4.);              // should be  0.028

  FillTables(m_nxbins, m_nbbins);
}

void EPA_Dipole::FillTables(const size_t& nx, const size_t& nb)
{
  axis xaxis(nx, 1.e-9, 1., axis_mode::log);
  axis baxis(nb, 1.e-3 * m_R, 1.e6 * m_R, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
  Fill_Nredx_Table();
  TestInvxb();
}

double EPA_Dipole::operator()(const double& x, const double& Q2)
{
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), prefC = m_mu2, prefD = 1.;
  // taking into account Q^2-dependence of form factors by over-riding their
  // value at Q^2 = 0.
  if (!m_approx) {
    prefC = m_mu2 / sqr(1. + Q2 / m_Lambda2);
    prefD = (4. * m_mass2 + Q2 * m_mu2) / (4. * m_mass2 + Q2) /
            sqr(1. + Q2 / m_Lambda2);
  }
  return (sqr(x) / 2. * prefC + (1. - x) * (1. + q2min / Q2) * prefD);
}

const double EPA_Dipole::phi(const double& y, const double& arg) const
{
  // phi_i(x) from Budnev et al., Eq. (D.7)
  double one_arg = 1. + arg;
  double sum1    = 1. / one_arg + 1. / (2. * one_arg * one_arg) +
                1. / (3. * one_arg * one_arg * one_arg);
  double sum2 = (m_bDip / one_arg + m_bDip * m_bDip / (2. * one_arg * one_arg) +
                 m_bDip * m_bDip * m_bDip / (3. * one_arg * one_arg * one_arg));
  double term1 = (1. + m_aDip * y) * (-log(one_arg / arg) + sum1);
  double term2 = (1. - m_bDip) * y / (4. * arg * pow(one_arg, 3));
  double term3 =
          m_cDip * (1. + y / 4.) * (log((one_arg - m_bDip) / one_arg) + sum2);
  return (term1 - term2 + term3);
}

double EPA_Dipole::N(const double& x)
{
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), q2max = Q2max(x);
  if (q2max <= q2min) return 0.;
  // taking into account the Q^2-dependence of form factors ...
  if (!m_approx) {
    double y = sqr(x) / (1. - x);
    return (1. - x) * (phi(y, q2max / m_Lambda2) - phi(y, q2min / m_Lambda2));
  }
  // ... or ignoring it.
  return ((1. - x + m_mu2 * sqr(x) / 2.) * log(q2max / q2min) -
          (1. - x) * (1. - q2min / q2max));
}

double EPA_Dipole::N(const double& x, const double& b)
{
  return (*p_N_xb)(x, b);
}

double EPA_Dipole::ReducedN(const double& x)
{
  return (*p_Nred_x)(x);
}

////////////////////////////////////////////////////////////////////////////////
//
// Gaussian form factors replacing the dipole form factors of Budnev et al.,
// Phys. Rep. C15 (1974) 181.
//
////////////////////////////////////////////////////////////////////////////////

EPA_Gauss::EPA_Gauss(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_Q02(1.)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_mu2         = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  m_Q02         = s["Q02"].GetTwoVector<double>()[b];
  // TODO ist m_Q02 in 1/GeV oder fm? Set defaults for mu2 and Q02 in EPA.C!
  if (m_beam == ATOOLS::Flavour(kf_p_plus)) {
    // m_mu2 = 2.79 * 2.79;
    // m_Q02 = 0.71;
  } else if (m_beam.IsIon()) {
    // m_mu2    = 0.;
    // m_Q02    = sqr(2. / m_R);
    m_pt2max = sqr(1. / m_R);
  }

  FillTables(m_nxbins, m_nbbins);
}

void EPA_Gauss::FillTables(const size_t& nx, const size_t& nb)
{
  double xmin = 1.e-9 / m_beam.GetAtomicNumber();
  double xmax = 1. / m_beam.GetAtomicNumber();
  double bmin = 1.e-3 * m_R, bmax = 1.e6 * m_R;
  axis   xaxis(nx, xmin, xmax, axis_mode::log);
  axis   baxis(nb, bmin, bmax, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
  Fill_Nredx_Table();
}

double EPA_Gauss::operator()(const double& x, const double& Q2)
{
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  // but modifying the dipole form to a Gaussian form
  double q2min = Q2min(x), prefC = m_mu2, prefD = 1.;
  // taking into account Q^2-dependence of form factors by over-riding their
  // value at Q^2 = 0 ...
  if (!m_approx) {
    prefC = m_mu2 * exp(-Q2 / m_Q02);
    prefD = exp(-Q2 / m_Q02);
  }
  // ... or ignoring it
  return (sqr(x) / 2. * prefC + (1. - x) * (1. + q2min / Q2) * prefD);
}

double EPA_Gauss::N(const double& x)
{
  // c.f. V.M. Budnev et al., Phys. Rep. C15(1974)181, Eq. (D.7)
  double q2min = Q2min(x), q2max = Q2max(x), term1, term2;
  // taking into account the Q^2-dependence of form factors and using that
  // Ei(-x) = -IncompleteGamma(0,x)
  if (!m_approx) {
    term1 = SF.IncompleteGamma(0, q2min / m_Q02) -
            SF.IncompleteGamma(0, q2max / m_Q02);
    term2 = term1 - q2min * (exp(-q2max / m_Q02) / q2max -
                             exp(-q2min / m_Q02) / q2min + 1. / m_Q02 * term1);
  }
  // ... or ignoring it.
  else {
    term1 = log(q2max / q2min);
    term2 = term1 + q2min * (1. / q2max - 1. / q2min);
  }
  double res = m_mu2 * sqr(x) * term1 + (1. - x) * term2;
  return res;
}

double EPA_Gauss::N(const double& x, const double& b)
{
  return (*p_N_xb)(x, b);
}

double EPA_Gauss::ReducedN(const double& x)
{
  return (*p_Nred_x)(x);
}

////////////////////////////////////////////////////////////////////////////////
//
// Form factor for ions in a Wood-Saxon potential.  A few comments:
//
// - Output for the density modifed to recover the atomic number A on
// integration,
//   and in units of fm^-3 or GeV^3. Internally we normalise it to unity,
//   because we will later multiply the square of the form factor with the
//   square of the charge Z.
// - Radius in 1/GeV, therefore "nucelar skin" m_d also in 1/GeV
//
////////////////////////////////////////////////////////////////////////////////

EPA_WoodSaxon::EPA_WoodSaxon(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_d(0.5 / rpa->hBar_c()), p_FF_Q2(nullptr),
      p_N(nullptr)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = dir > 0 ? 0 : 1;
  m_d           = s["WoodSaxon_d"].GetTwoVector<double>()[b];
  if (!m_beam.IsIon())
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
  m_pt2max = sqr(1. / m_R);

  FillTables(m_nxbins, m_nbbins);
}

void EPA_WoodSaxon::FillTables(const size_t& nx, const size_t& nb)
{
  m_rho0 = CalculateDensity();
  InitFFTable(1.e-12, 1.e4);
  InitNTable(1.e-10, 1.);
  double xmin = 1.e-9 / m_beam.GetAtomicNumber();
  double xmax = 1. / m_beam.GetAtomicNumber();
  double bmin = 1.e-3 * m_R, bmax = 1.e6 * m_R;
  axis   xaxis(nx, xmin, xmax, axis_mode::log);
  axis   baxis(nb, bmin, bmax, axis_mode::log);
  Fill_Nxb_Table(xaxis, baxis);
  Fill_Invxb_Table();
  Fill_Nredx_Table();
}

void EPA_WoodSaxon::InitFFTable(const double& q2min, const double& q2max)
{
  p_FF_Q2 = new OneDim_Table(axis(100000, q2min, q2max, axis_mode::log));
  WS_potential*    ws = new WS_potential(m_R, m_d);
  Gauss_Integrator gauss(ws);
  for (size_t i = 0; i < p_FF_Q2->GetAxis().m_nbins; i++) {
    ws->SetQ(sqrt(p_FF_Q2->GetAxis().x(i)));
    double rmin = 0., rmax = m_R;
    double res = m_rho0 * gauss.Integrate(rmin, rmax, 1.e-6, 0), inc = 0.;
    do {
      rmin = rmax;
      rmax *= 2.;
      res += inc = m_rho0 * gauss.Integrate(rmin, rmax, 1.e-6, 0);
    } while (rmax < 4. * m_R || dabs(inc / res) > 1.e-6);
    p_FF_Q2->Fill(i, res);
  }
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
  do {
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

double EPA_WoodSaxon::N(const double& x)
{
  return (*p_N)(x);
}

double EPA_WoodSaxon::N(const double& x, const double& b)
{
  return (*p_N_xb)(x, b);
}

double EPA_WoodSaxon::ReducedN(const double& x)
{
  return (*p_Nred_x)(x);
}
