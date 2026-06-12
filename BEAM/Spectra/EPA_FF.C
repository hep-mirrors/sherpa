#include "BEAM/Spectra/EPA_FF.H"

#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Filon_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Table_Cache.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_File.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Shell_Tools.H"

#include <cmath>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <vector>
using namespace BEAM;
using namespace ATOOLS;

namespace {
  // RAII wrapper that opens the Sherpa result directory as a My_File zip-DB for
  // the duration of the EPA cache I/O. While the DB is open, the cache tables
  // are stored inside <root>.zip -- My_Out_File::Close prefix-matches the write
  // path (root + "/EPA/<file>.tab") against the registered DB and routes it into
  // the archive -- alongside the matrix-element results, instead of as loose
  // files. OpenDB/CloseDB are MPI-collective (they broadcast), so this guard
  // must be constructed and destroyed on every rank in lock-step; every rank
  // builds the same spectra with the same settings, so that holds. CloseDB (in
  // the destructor) flushes the archive back to disk.
  struct ResultDB_Guard {
    std::string m_db;
    bool m_active;
    ResultDB_Guard(std::string db, bool active)
        : m_db(std::move(db)), m_active(active)
    {
      if (m_active) My_In_File::OpenDB(m_db);
    }
    ~ResultDB_Guard()
    {
      if (m_active) My_In_File::CloseDB(m_db);
    }
    ResultDB_Guard(const ResultDB_Guard&) = delete;
    ResultDB_Guard& operator=(const ResultDB_Guard&) = delete;
  };
} // namespace

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
    : //////////////////////////////////////////////////////////////////////////////
      //
      // Initialisation of relevant beam parameters:
      // note that the particle radius is in fm and transformed into 1/GeV
      //
      //////////////////////////////////////////////////////////////////////////////
      m_beam(beam), m_A(beam.IsIon() ? beam.GetMassNumber() : 1),
      m_mass(beam.Mass(true) / m_A), m_mass2(m_mass * m_mass),
      m_R(beam.Radius() / rpa->hBar_c()), m_q2min(-1.), m_q2max(1.),
      m_pt2max(-1.),
      m_Zsquared(beam.IsIon() ? sqr(m_beam.GetAtomicNumber()) : 1.), m_b(0.),
      p_N_xb(nullptr)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t b = dir > 0 ? 0 : 1;
  m_q2max = s["Q2Max"].GetTwoVector<double>()[b];
  m_q2min = s["Q2Min"].GetTwoVector<double>()[b];
  m_nxbins = s["xBins"].GetTwoVector<int>()[b];
  m_nbbins = s["bBins"].GetTwoVector<int>()[b];
  m_xmin = s["xMin"].GetTwoVector<double>()[b];
  m_xmax = s["xMax"].GetTwoVector<double>()[b];
  m_bmin = s["bMin"].GetTwoVector<double>()[b];
  m_b_pl_threshold = s["bThreshold"].GetTwoVector<double>()[b] * m_R;
  m_bmax = s["bMax"].GetTwoVector<double>()[b];
  m_chi_cut = s["chiMax"].GetTwoVector<double>()[b];

  if (m_bmin <= 0. || m_bmin > m_bmax)
    THROW(invalid_input, "Unphysical input for EPA impact parameter. ");
  if (m_chi_cut <= 0.)
    THROW(invalid_input, "Unphysical input for EPA chiMax (must be > 0). ");
  if (m_xmin <= 0. || m_xmin > m_xmax || m_xmax > 1.)
    THROW(invalid_input, "Unphysical input for EPA x-limits. ");
  if (m_q2min > m_q2max)
    THROW(invalid_input, "Unphysical input for EPA Q2-limits. ");

  // Resolve the RESULT_DIRECTORY root (e.g. <config>/Results). A relative
  // RESULT_DIRECTORY is taken w.r.t. the config path, exactly as
  // Matrix_Element_Handler does, so m_respath matches the path the ME handler
  // opens as a zip-DB and the EPA cache lands inside the very same
  // <RESULT_DIRECTORY>.zip (BuildNxbTable opens that DB around the cache I/O).
  //
  // Caching is only possible when result files are written at all: if the user
  // disabled GENERATE_RESULT_DIRECTORY, or RESULT_DIRECTORY resolves to empty,
  // leave m_respath empty, which LoadOrBuild/BuildNxbTable treat as "caching
  // disabled" (no DB opened, no read, no write). GENERATE_RESULT_DIRECTORY is
  // read with its own default so this does not depend on Matrix_Element_Handler
  // having registered it yet.
  m_cache = s["CacheTables"].SetDefault(true).Get<bool>();
  Settings& mainsettings = Settings::GetMainSettings();
  const bool storeresults =
      mainsettings["GENERATE_RESULT_DIRECTORY"].Get<bool>();
  std::string res =
      ShortenPathName(mainsettings["RESULT_DIRECTORY"].Get<std::string>());
  if (storeresults && !res.empty()) {
    if (res[0] != '/' && mainsettings.GetPath() != "")
      res = mainsettings.GetPath() + "/" + res;
    m_respath = res;
  }
  // else: m_respath stays empty -> caching disabled for this beam
}

std::string EPA_FF_Base::BaseCacheId(const std::string& tag) const
{
  std::ostringstream s;
  s.precision(17);
  s << "v" << s_cache_version << "|" << tag << "|A" << m_A << "|Z2"
    << m_Zsquared << "|m" << m_mass << "|R" << m_R << "|x[" << m_nxbins << ","
    << m_xmin << "," << m_xmax << "]|b[" << m_nbbins << "," << m_bmin << ","
    << m_b_pl_threshold << "," << m_bmax << "]";
  return s.str();
}

std::string EPA_FF_Base::CacheFileName(const std::string& tag,
                                       const std::string& key) const
{
  return m_beam.IDName() + "_" + tag + "_" +
         std::to_string(std::hash<std::string>{}(key)) + ".tab";
}

template <class Table>
void EPA_FF_Base::LoadOrBuild(const std::string& tag, const std::string& key,
                             std::unique_ptr<Table>& target,
                             const std::function<void()>& fill)
{
  const bool docache = CachingActive();
  // The tables live as entries EPA/<file>.tab inside <m_respath>.zip; the caller
  // (BuildNxbTable) keeps that zip-DB open for the duration, so Table_Cache's
  // My_In_File/My_Out_File resolve against it (see ResultDB_Guard) rather than
  // touching loose files.
  const std::string path = m_respath + "/EPA/" + CacheFileName(tag, key);
  // Table_Cache::Load is collective (rank 0 reads and broadcasts the content),
  // so every rank sees the same content and reaches the same hit/miss decision
  // -- no shared filesystem and no manual rank coordination needed.
  if (docache) {
    auto t = Table_Cache::Load<Table>(path, key);
    if (t) {
      target = std::move(t);
      return;
    }
  }
  // Miss: every rank runs fill() (which populates target) exactly as the
  // pre-cache code did; Table_Cache::Save writes the cache on rank 0 only, so
  // the next run is a load.
  fill();
  if (docache) Table_Cache::Save(path, key, *target);
}

void EPA_FF_Base::BuildNxbTable()
{
  // Open the Sherpa result directory as a zip-DB for the whole build, so both
  // this N(x,b) table and any prerequisite (the Woods-Saxon form-factor table,
  // built by PrepareFill) are stored inside <RESULT_DIRECTORY>.zip rather than
  // as loose files. This is the single place the DB is opened: BuildFFQTable
  // runs inside PrepareFill below, so one scope covers both tables with no
  // nested OpenDB. On a cache hit the open DB also makes the entries readable.
  const bool docache = CachingActive();
  ResultDB_Guard db(m_respath + "/", docache);
  // PrepareFill() builds any prerequisite (e.g. the Woods-Saxon form-factor
  // table) before FillTables() fills the N(x,b) grid.
  LoadOrBuild<TwoDim_Table>("nxb", CacheId(), p_N_xb, [this] {
    PrepareFill();
    FillTables();
  });
}

void EPA_FF_Base::FillTables()
{
  axis xaxis(m_nxbins, m_xmin, m_xmax, axis_mode::log);
  axis baxis(m_nbbins, m_bmin * m_R, std::min(m_b_pl_threshold, m_bmax * m_R),
             axis_mode::log);

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
  p_N_xb = std::make_unique<TwoDim_Table>(xaxis, baxis);
  N_xb_int* kernel = new N_xb_int(this);
  Bessel_Integrator bessel(kernel, 1);
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    for (size_t j = 0; j < baxis.m_nbins; j++) {
      msg_Debugging() << METHOD << ": Filling table for x = " << xaxis.x(i)
                      << ", and b = " << baxis.x(j) << "\n";
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
  msg_Out() << METHOD << ": Writing output for " << m_beam.IDName() << " and "
            << type << "\n";

  double step_x(std::log(m_xmax / m_xmin) / double(m_nxbins));
  std::vector<double> xs(m_nxbins);
  xs[0] = m_xmin;
  for (int i = 1; i < m_nxbins; ++i) {
    xs[i] = xs[i - 1] * std::exp(step_x);
  }

  double step_b(std::log(m_bmax / m_bmin) / double(m_nbbins));
  std::vector<double> bs(m_nbbins);
  bs[0] = m_bmin * m_R;
  for (int i = 1; i < m_nbbins; ++i) {
    bs[i] = bs[i - 1] * std::exp(step_b);
  }

  double q2max(1e4), q2min(1.e-12);
  int nq2steps(10000);
  double step_q2(std::log(q2max / q2min) / double(nq2steps));
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
    for (auto& x : xs) {
      outfile_Nxb << x << ",0," << N(x, 0.) << std::endl;
    }
    outfile_Nxb.close();
  }

  std::ofstream outfile_FFq2(prefix + "FF_q2.csv");
  outfile_FFq2 << "q2,FF" << std::endl;
  for (auto& q2 : q2s)
    outfile_FFq2 << q2 << "," << FF(q2) << std::endl;
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
  m_b = rpa->hBar_c() / m_mass / (2. * M_PI / 137.);
}

double EPA_Point::N(const double& x, const double& ran)
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
  m_b = rpa->hBar_c() / m_mass / (2. * M_PI / 137.);
}

double EPA_PointApprox::N(const double& x, const double& ran)
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
  size_t b = dir > 0 ? 0 : 1;
  m_mu2 = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  m_Q02 = s["Q02"].GetTwoVector<double>()[b];
  if (m_beam.Kfcode() != kf_p_plus)
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
  m_b = 1.001 * m_R;
}

double EPA_Proton::N(const double& x, const double& ran)
{
  auto phi = [this](double x, double z) {
    double y = x * x / (1. - x);
    double a = (1. + m_mu2) / 4. + 4. * m_mass2 / m_Q02;
    double b = 1. - 4. * m_mass2 / m_Q02;
    double c = (m_mu2 - 1.) * std::pow(b, -4);
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
  size_t b = dir > 0 ? 0 : 1;
  m_mu2 = sqr(s["MagneticMu"].GetTwoVector<double>()[b]);
  if (m_beam.Kfcode() != kf_p_plus)
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());
  m_b = 1.001 * m_R;
}

double EPA_ProtonApprox::N(const double& x, const double& ran)
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
  size_t b = dir > 0 ? 0 : 1;
  m_Q02 = s["Q02"].GetTwoVector<double>()[b];
  BuildNxbTable();
}

std::string EPA_Gauss::CacheId() const
{
  std::ostringstream s;
  s.precision(17);
  s << "|Q02" << m_Q02;
  return BaseCacheId("Gauss") + s.str();
}

////////////////////////////////////////////////////////////////////////////////
// Form factor of a homogeneously charged sphere (HCS) of radius R, i.e. the
// Fourier transform of a uniform charge density:
//   F(qR) = 3 (sin(qR) - qR cos(qR)) / (qR)^3,
// with q = sqrt(Q^2) and R = m_R (the beam radius in 1/GeV), so the argument
// qR is dimensionless and F(0) = 1.
////////////////////////////////////////////////////////////////////////////////

EPA_HCS::EPA_HCS(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  BuildNxbTable();
}

std::string EPA_HCS::CacheId() const { return BaseCacheId("HCS"); }

////////////////////////////////////////////////////////////////////////////////
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_Dipole::EPA_Dipole(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t b = dir > 0 ? 0 : 1;
  m_Q02 = s["Q02"].GetTwoVector<double>()[b];

  BuildNxbTable();
}

std::string EPA_Dipole::CacheId() const
{
  std::ostringstream s;
  s.precision(17);
  s << "|Q02" << m_Q02;
  return BaseCacheId("Dipole") + s.str();
}

////////////////////////////////////////////////////////////////////////////////
// Dipole form factors and related functions in different approximations, all
// based on Budnev et al., Phys. Rep. C15 (1974) 181.
////////////////////////////////////////////////////////////////////////////////

EPA_DipoleApprox::EPA_DipoleApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  BuildNxbTable();
}

std::string EPA_DipoleApprox::CacheId() const
{
  return BaseCacheId("DipoleApprox");
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

EPA_WoodSaxon::EPA_WoodSaxon(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir), m_d(0.55 / rpa->hBar_c()),
      m_R_WS(6.49 / rpa->hBar_c()), m_rnodes(1024), m_rmaxfactor(16.)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t b = dir > 0 ? 0 : 1;
  m_R_WS = s["WoodsSaxon_R"].GetTwoVector<double>()[b] / rpa->hBar_c();
  m_d = s["WoodsSaxon_d"].GetTwoVector<double>()[b] / rpa->hBar_c();
  // Validate before the int->size_t assignment: a negative rNodes would wrap to
  // a huge size_t and make InitFFTable attempt a runaway allocation.
  const int rnodes = s["WoodsSaxon_rNodes"].GetTwoVector<int>()[b];
  if (rnodes < 2)
    THROW(invalid_input, "EPA WoodsSaxon_rNodes must be >= 2.");
  m_rnodes = static_cast<size_t>(rnodes);
  m_rmaxfactor = s["WoodsSaxon_rMaxFactor"].GetTwoVector<double>()[b];
  if (m_rmaxfactor <= 0.)
    THROW(invalid_input, "EPA WoodsSaxon_rMaxFactor must be > 0.");
  if (!m_beam.IsIon())
    THROW(fatal_error, "Wrong form factor for " + m_beam.IDName());

  // Try the N(x,b) cache first: on a hit the (expensive) form-factor table and
  // the density are never needed at run time and are skipped entirely. On a
  // miss PrepareFill() builds the form-factor table. The debug CSV output
  // (OutputSpectra/All) reaches FF(), which lazily builds p_FF_Q if an N(x,b)
  // hit skipped it, so no extra build is needed here.
  BuildNxbTable();
}

std::string EPA_WoodSaxon::WSParamString() const
{
  std::ostringstream s;
  s.precision(17);
  s << "|RWS" << m_R_WS << "|d" << m_d << "|q[" << m_q_n << "," << m_q_min << ","
    << m_q_max << "]|rN" << m_rnodes << "|rMax" << m_rmaxfactor;
  return s.str();
}

std::string EPA_WoodSaxon::CacheId() const
{
  return BaseCacheId("WoodSaxon") + WSParamString();
}

void EPA_WoodSaxon::BuildFFQTable()
{
  if (p_FF_Q) return; // already built
  const std::string key =
      "v" + std::to_string(s_cache_version) + "|WSFF" + WSParamString();
  LoadOrBuild<OneDim_Table>("wsff", key, p_FF_Q, [this] { InitFFTable(); });
}

void EPA_WoodSaxon::InitFFTable()
{
  // Filon sine transform of the Woods-Saxon density:
  //   F(q)  = rho0/q * \int_0^rmax sin(q r) r rho(r) dr,
  //   rho(r) = 1/(1+exp((r-R_WS)/d)),  rho0 = 1 / \int_0^rmax r^2 rho dr,
  // so that F(0) = 1 by construction. The shared r-grid only has to resolve
  // rho(r); the oscillation is integrated analytically, so F(q) is accurate
  // over the whole q-range (unlike a plain quadrature, which would need the
  // grid to resolve sin(q r) itself).
  const double rmax = m_R_WS + m_rmaxfactor * m_d;
  size_t intervals = m_rnodes < 2 ? 2 : m_rnodes;
  if (intervals % 2) ++intervals; // Filon-Simpson needs an even interval count
  const size_t nnodes = intervals + 1;
  const double h = rmax / double(intervals);

  msg_Out() << METHOD << ": building Woods-Saxon form factor (" << m_q_n
            << " q-bins, " << intervals << " r-intervals to r = " << rmax
            << " 1/GeV).\n";

  std::vector<double> g(nnodes); // g(r) = r * rho(r)
  double norm = 0.;              // \int r^2 rho dr via composite Simpson
  for (size_t j = 0; j < nnodes; ++j) {
    const double r = double(j) * h;
    const double rho = 1. / (1. + std::exp((r - m_R_WS) / m_d));
    g[j] = r * rho;
    const double w = (j == 0 || j == nnodes - 1) ? 1. : (j % 2 ? 4. : 2.);
    norm += w * (r * r * rho);
  }
  norm *= h / 3.;
  const double rho0 = 1. / norm;

  Filon_Integrator filon(std::move(g), 0., rmax);
  p_FF_Q = std::make_unique<OneDim_Table>(
      axis(m_q_n, m_q_min, m_q_max, axis_mode::linear));
  for (size_t i = 0; i <= static_cast<size_t>(m_q_n); ++i) {
    const double q = p_FF_Q->GetAxis().x(i);
    const double ff = (q == 0.) ? 1. : rho0 * filon.SineTransform(q) / q;
    p_FF_Q->Fill(i, ff);
  }
}

////////////////////////////////////////////////////////////////////////////////
// Approximation of the Woods-Saxon the Woods-Saxon distribution as a hard
// sphere, with radius RA, convoluted with a Yukawa potential with
// range a = 0.7 fm
////////////////////////////////////////////////////////////////////////////////

EPA_WoodSaxonApprox::EPA_WoodSaxonApprox(const ATOOLS::Flavour& beam,
                                         const int dir)
    : EPA_FF_Base(beam, dir)
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t b = dir > 0 ? 0 : 1;
  m_R_WS = s["WoodsSaxon_R"].GetTwoVector<double>()[b] / rpa->hBar_c();
  m_a = s["WoodsSaxonApprox_a"].GetTwoVector<double>()[b] / rpa->hBar_c();

  BuildNxbTable();
}

std::string EPA_WoodSaxonApprox::CacheId() const
{
  std::ostringstream s;
  s.precision(17);
  s << "|RWS" << m_R_WS << "|a" << m_a;
  return BaseCacheId("WoodSaxonApprox") + s.str();
}

////////////////////////////////////////////////////////////////////////////////
// Class for the Ion FF in the Electric Dipole approximation using
// the analytical expression from 2207.03012, eq. 3.3
////////////////////////////////////////////////////////////////////////////////

EPA_IonApprox::EPA_IonApprox(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  m_b_pl_threshold = m_R;
  EPA_IonApprox::FillTables();
}

void EPA_IonApprox::FillTables()
{
  // This is NOT used for the calculation of the flux in the integration,
  // instead only filled for debugging purposes using the "OutputAllSpectra"
  // setting
  axis xaxis(m_nxbins, m_xmin, m_xmax, axis_mode::log);
  axis baxis(m_nbbins, m_bmin * m_R, m_bmax * m_R, axis_mode::log);

  msg_Out() << METHOD << " in " << xaxis.m_nbins << " * " << baxis.m_nbins
            << " bins:\n"
            << "   x in [" << xaxis.m_xmin << ", " << xaxis.m_xmax << "], "
            << "b in [" << baxis.m_xmin << ", " << baxis.m_xmax << "], "
            << "from R = " << m_R << " 1/GeV = " << (m_R * rpa->hBar_c())
            << " fm.\n";
  p_N_xb = std::make_unique<TwoDim_Table>(xaxis, baxis);
  for (size_t i = 0; i < xaxis.m_nbins; i++) {
    for (size_t j = 0; j < baxis.m_nbins; j++) {
      double chi = xaxis.x(i) * m_mass * baxis.x(j);
      double val = 2 * m_Zsquared * baxis.x(j) * xaxis.x(i) * m_mass2 *
                   ATOOLS::sqr(ATOOLS::SF.Kn(1, chi));
      p_N_xb->Fill(i, j, val);
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
// Class for the Ion FF in the Electric Dipole approximation integrated in b
////////////////////////////////////////////////////////////////////////////////

EPA_IonApproxIntegrated::EPA_IonApproxIntegrated(const ATOOLS::Flavour& beam,
                                                 const int dir)
    : EPA_FF_Base(beam, dir)
{
  m_b_pl_threshold = m_R;
  p_N_xb.reset();
  m_b = 1.001 * m_R;
}

double EPA_IonApproxIntegrated::N(const double& x, const double& ran)
{
  double chi = x * m_mass * m_R * Max(1., m_bmin);
  return 2 * m_Zsquared / x *
         (chi * SF.Kn(0, chi) * SF.Kn(1, chi) -
          chi * chi / 2. * (sqr(SF.Kn(1, chi)) - sqr(SF.Kn(0, chi))));
}

////////////////////////////////////////////////////////////////////////////////
// Dummy test class for checking the Bessel_Integrator
////////////////////////////////////////////////////////////////////////////////

EPA_Test::EPA_Test(const ATOOLS::Flavour& beam, const int dir)
    : EPA_FF_Base(beam, dir)
{
  FillTables();
}
