#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/EPA_FF.H"
#include "BEAM/Spectra/EPA_Spectra_Plotter.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"

#include <fstream>
#include <string>

using namespace BEAM;
using namespace ATOOLS;
using string = std::string;

EPA::EPA(const Flavour& beam, const double energy, const double pol,
         const int dir)
    : Beam_Base(beamspectrum::EPA, beam, energy, pol, dir),
      m_type(EPA_ff_type::point), p_ff(nullptr), m_mass(m_beam.Mass(true)),
      m_charge(m_beam.Charge()), m_plotting(0)
{
  if (m_beam.Charge() == 0.)
    THROW(fatal_error,
          "No photon flux for uncharged particles. Can not enable EPA. ")
  Initialise();
  m_Nbunches = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_photon);
  m_bunches[1] = m_beam;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0., 0., 0., 0.);
  m_on         = true;
}

bool EPA::CalculateWeight(double x, double q2)
{
  m_x      = x;
  m_q2     = q2;
  m_weight = (x >= m_xmin && x <= m_xmax)
                     ? ATOOLS::Max(0., m_pref / x * p_ff->N(x))
                     : 0.;
  if (IsNan(m_weight))
    msg_Out() << "Boink! " << METHOD << "(x = " << x << ") yields NaN.\n";
  return true;
}

void EPA::FixPosition()
{
  double R = p_ff->SelectB(m_x), phi = 2. * M_PI * ran->Get();
  m_position = R * Vec4D(0., cos(phi), sin(phi), 0.);
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D& out, const size_t& i)
{
  if (i == 0) {
    m_vecouts[0] = out;
    m_vecouts[1] = m_lab - out;
    m_x          = out[0] / m_lab[0];
    m_q2         = dabs(out.Abs2());
  }
}

void EPA::Initialise()
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  RegisterDefaults();
  size_t b   = m_dir > 0 ? 0 : 1;
  m_aqed     = s["AlphaQED"].Get<double>();
  m_pref     = sqr(m_charge) * m_aqed / M_PI;
  m_plotting = s["PlotSpectra"].Get<bool>();
  m_pt2max   = !m_beam.IsIon()
                       ? sqr(m_energy * s["ThetaMax"].GetTwoVector<double>()[b])
                       : sqr(1. / m_beam.Radius());
  m_xmin     = s["xMin"].GetTwoVector<double>()[b];
  m_xmax     = s["xMax"].GetTwoVector<double>()[b];
  m_bmin     = s["bMin"].GetTwoVector<double>()[b];
  m_bmax     = s["bMax"].GetTwoVector<double>()[b];

  if (m_plotting) {
    Tests();
    /*EPA_Spectra_Plotter* plotter =
            new EPA_Spectra_Plotter(this, string("Spectra"));
    (*plotter)(99);
    delete plotter;*/
    THROW(normal_exit, "Tests done.");
  }

  m_type = static_cast<EPA_ff_type>(s["Form_Factor"].GetTwoVector<int>()[b]);
  switch (m_type) {
    case EPA_ff_type::point: p_ff = new EPA_Point(m_beam, m_dir); break;
    case EPA_ff_type::Gauss: p_ff = new EPA_Gauss(m_beam, m_dir); break;
    case EPA_ff_type::dipole: p_ff = new EPA_Dipole(m_beam, m_dir); break;
    case EPA_ff_type::WoodSaxon: p_ff = new EPA_WoodSaxon(m_beam, m_dir); break;
    default: THROW(not_implemented, "unknown EPA form factor. ");
  }
  p_ff->SetPT2Max(m_pt2max);
}

void EPA::RegisterDefaults() const
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  s["Q2Max"].SetDefault(3.0);
  s["Q2Min"].SetDefault(-1.);
  s["xMax"].SetDefault(1.);
  s["xMin"].SetDefault(1.e-6);
  s["xBins"].SetDefault(12);
  s["bMin"].SetDefault(0.);
  s["bMax"].SetDefault(1.e12);
  s["bBins"].SetDefault(10);
  s["Form_Factor"].SetDefault(size_t(m_beam.IsIon() ? EPA_ff_type::WoodSaxon
                                     : m_beam.IsNucleon() ? EPA_ff_type::dipole
                                     : m_beam.IsMeson()   ? EPA_ff_type::dipole
                                                        : EPA_ff_type::point));
  s["MagneticMu"].SetDefault(m_beam.IsNucleon() ? 2.79 : 0.);
  // TODO check the default for ions below
  s["Q02"].SetDefault(m_beam.IsNucleon() ? 0.71 : sqr(2. / m_beam.Radius()));
  s["WoodSaxon_d"].SetDefault(0.5);
  s["AlphaQED"].SetDefault(0.0072992701);
  s["ThetaMax"].SetDefault(0.3);
  s["Approximation"].SetDefault(false);
  s["PlotSpectra"].SetDefault(false);
}

void EPA::Tests()
{
  // Point
  auto* ff_e = new EPA_Point(Flavour(kf_e), 0);
  ff_e->OutputToCSV("point");
  auto* ff_p_point = new EPA_Point(Flavour(kf_p_plus), 0);
  ff_p_point->OutputToCSV("point");

  // Dipole
  auto* ff_p_dip = new EPA_Dipole(Flavour(kf_p_plus), 0);
  ff_p_dip->OutputToCSV("dipole");

  // Gauss
  auto* ff_p_gauss = new EPA_Gauss(Flavour(kf_p_plus), 0);
  ff_p_gauss->OutputToCSV("gauss");
  auto* ff_lead_gauss = new EPA_Gauss(Flavour(kf_lead208), 0);
  ff_lead_gauss->OutputToCSV("gauss");
  auto* ff_calc_gauss = new EPA_Gauss(Flavour(kf_calcium40), 0);
  ff_calc_gauss->OutputToCSV("gauss");

  // Woods-Saxon
  auto* ff_lead_ws = new EPA_WoodSaxon(Flavour(kf_lead208), 0);
  ff_lead_ws->OutputToCSV("ws");
  auto* ff_calc_ws = new EPA_WoodSaxon(Flavour(kf_calcium40), 0);
  ff_calc_ws->OutputToCSV("ws");

  return;
}
