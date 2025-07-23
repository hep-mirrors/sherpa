#include "BEAM/Spectra/EPA.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Settings.H"

#include <string>

using namespace BEAM;
using namespace ATOOLS;
using string = std::string;

EPA::EPA(const Flavour& beam, const double energy, const double pol,
         const int dir)
    : Beam_Base(beamspectrum::EPA, beam, energy, pol, dir),
      m_fftype(EPA_ff_type::point), p_ff(nullptr), m_mass(beam.Mass(true)),
      m_aqed(1. / 127), m_pref(0.), m_q2(0.), m_pt2max(-1.), m_xmin(0.),
      m_xmax(1.), m_plotting(0)
{
  if (m_beam.Charge() == 0.)
    THROW(fatal_error,
          "No photon flux for uncharged particles. Can not enable EPA. ")
  m_Nbunches = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = Flavour(kf_photon);
  m_bunches[1] = m_beam;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0., 0., 0., 0.);
  m_on         = true;
  m_pref       = m_aqed / M_PI;

  RegisterDefaults();
  Initialise();
  if (m_plotting) {
    Tests();
    THROW(normal_exit, "Tests done.");
  }
}

bool EPA::CalculateWeight(double x, double q2)
{
  m_x      = x;
  m_q2     = q2;
  m_weight = (x > m_xmin && x < m_xmax) ? ATOOLS::Max(0., m_pref * p_ff->N(x))
                                        : 0.;
  if (IsNan(m_weight))
    msg_Out() << "Boink! " << METHOD << "(x = " << x << ") yields NaN.\n";
  return true;
}

void EPA::FixPosition()
{
  double phi = 2. * M_PI * ran->Get();
  m_position = p_ff->ImpactParameter() * Vec4D(0., cos(phi), sin(phi), 0.);
}

void EPA::SetOutMomentum(const ATOOLS::Vec4D& out)
{
  m_vecouts[0] = out;
  m_vecouts[1] = m_lab - out;
  m_q2         = dabs(out.Abs2());
}

void EPA::Initialise()
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  size_t      b = m_dir > 0 ? 0 : 1;
  m_aqed        = s["AlphaQED"].Get<double>();
  m_plotting    = s["PlotSpectra"].Get<bool>();
  m_pt2max      = !m_beam.IsIon()
                          ? sqr(m_energy * s["ThetaMax"].GetTwoVector<double>()[b])
                          : sqr(rpa->hBar_c() / m_beam.Radius());
  m_xmin        = s["xMin"].GetTwoVector<double>()[b];
  m_xmax        = s["xMax"].GetTwoVector<double>()[b];

  m_fftype = static_cast<EPA_ff_type>(s["Form_Factor"].GetTwoVector<int>()[b]);
  switch (m_fftype) {
    case EPA_ff_type::point:
      p_ff = new EPA_Point(m_beam, m_dir); break;
    case EPA_ff_type::pointApprox:
      p_ff = new EPA_PointApprox(m_beam, m_dir); break;
    case EPA_ff_type::proton:
      p_ff = new EPA_Proton(m_beam, m_dir); break;
    case EPA_ff_type::protonApprox:
      p_ff = new EPA_ProtonApprox(m_beam, m_dir); break;
    case EPA_ff_type::Gauss:
      p_ff = new EPA_Gauss(m_beam, m_dir); break;
    case EPA_ff_type::hcs:
      p_ff = new EPA_HCS(m_beam, m_dir); break;
    case EPA_ff_type::dipole:
      p_ff = new EPA_Dipole(m_beam, m_dir); break;
    case EPA_ff_type::dipoleApprox:
      p_ff = new EPA_DipoleApprox(m_beam, m_dir); break;
    case EPA_ff_type::ionApprox:
      p_ff = new EPA_IonApprox(m_beam, m_dir); break;
    case EPA_ff_type::WoodSaxon:
      p_ff = new EPA_WoodSaxon(m_beam, m_dir); break;
    default: THROW(not_implemented, "unknown EPA form factor. ");
  }
  p_ff->SetPT2Max(m_pt2max);
}

void EPA::RegisterDefaults() const
{
  const auto& s = Settings::GetMainSettings()["EPA"];
  s["Q2Max"].SetDefault(1.);
  s["Q2Min"].SetDefault(-1.);
  s["xMax"].SetDefault(1.);
  s["xMin"].SetDefault(m_beam.IsIon() ? 1.e-10 : 1.e-6);
  s["xBins"].SetDefault(50);
  s["bMin"].SetDefault(0.1);
  s["bThreshold"].SetDefault(4.);
  s["bMax"].SetDefault(1.e2);
  s["bBins"].SetDefault(50);
  s["Form_Factor"].SetDefault(size_t(m_beam.IsIon() ? EPA_ff_type::WoodSaxon
                                     : m_beam.IsNucleon() ? EPA_ff_type::dipole
                                     : m_beam.IsMeson()   ? EPA_ff_type::dipole
                                                        : EPA_ff_type::point));
  s["MagneticMu"].SetDefault(m_beam.IsNucleon() ? 2.79 : 0.);
  s["Q02"].SetDefault(m_beam.IsNucleon()
                              ? 0.71
                              : sqr(2. / m_beam.Radius() * rpa->hBar_c()));
  s["WoodSaxon_d"].SetDefault(0.5);
  s["AlphaQED"].SetDefault(0.0072992701);
  s["ThetaMax"].SetDefault(0.3);
  s["PlotSpectra"].SetDefault(false);
}

void EPA::Tests()
{
  // Test
  //auto* ff_test = new EPA_Test(Flavour(kf_photon), 0);
  //ff_test->OutputToCSV("test");
  //delete ff_test;

  // Lepton
  auto* ff_e = new EPA_Point(Flavour(kf_e), 0);
  ff_e->OutputToCSV("point");
  delete ff_e;
  auto* ff_eApprox = new EPA_PointApprox(Flavour(kf_e), 0);
  ff_eApprox->OutputToCSV("pointApprox");
  delete ff_eApprox;

  // Proton
  // ======
  // Proton Sachs
  auto* ff_p_proton = new EPA_Proton(Flavour(kf_p_plus), 0);
  ff_p_proton->OutputToCSV("proton");
  delete ff_p_proton;
  // Proton Sachs Approx
  auto* ff_p_protonApprox = new EPA_ProtonApprox(Flavour(kf_p_plus), 0);
  ff_p_protonApprox->OutputToCSV("protonApprox");
  delete ff_p_protonApprox;
  // Gauss
  auto* ff_p_gauss = new EPA_Gauss(Flavour(kf_p_plus), 0);
  ff_p_gauss->OutputToCSV("gauss");
  delete ff_p_gauss;
  // HCS
  auto* ff_p_hcs = new EPA_HCS(Flavour(kf_p_plus), 0);
  ff_p_hcs->OutputToCSV("hcs");
  delete ff_p_hcs;
  // Dipole
  auto* ff_p_dip = new EPA_Dipole(Flavour(kf_p_plus), 0);
  ff_p_dip->OutputToCSV("dipole");
  delete ff_p_dip;
  // DipoleApprox
  auto* ff_p_dipApprox = new EPA_DipoleApprox(Flavour(kf_p_plus), 0);
  ff_p_dipApprox->OutputToCSV("dipoleApprox");
  delete ff_p_dipApprox;
  // Ion Approx
  auto* ff_p_ionApprox = new EPA_IonApprox(Flavour(kf_p_plus), 0);
  ff_p_ionApprox->OutputToCSV("ionApprox");
  delete ff_p_ionApprox;

  std::vector<kf_code> ions({kf_lead208, kf_calcium40});
  for (kf_code ion : ions) {
    // Gauss
    auto* ff_ion_gauss = new EPA_Gauss(Flavour(ion), 0);
    ff_ion_gauss->OutputToCSV("gauss");
    delete ff_ion_gauss;
    // HCS
    auto* ff_ion_hcs = new EPA_HCS(Flavour(ion), 0);
    ff_ion_hcs->OutputToCSV("hcs");
    delete ff_ion_hcs;
    // Dipole
    auto* ff_ion_dip = new EPA_Dipole(Flavour(ion), 0);
    ff_ion_dip->OutputToCSV("dipole");
    delete ff_ion_dip;
    // DipoleApprox
    auto* ff_ion_dipApprox = new EPA_DipoleApprox(Flavour(ion), 0);
    ff_ion_dipApprox->OutputToCSV("dipoleApprox");
    delete ff_ion_dipApprox;
    // Woods-Saxon
    auto* ff_ion_ws = new EPA_WoodSaxon(Flavour(ion), 0);
    ff_ion_ws->OutputToCSV("ws");
    delete ff_ion_ws;
    // Ion Approx
    auto* ff_ionApprox = new EPA_IonApprox(Flavour(ion), 0);
    ff_ionApprox->OutputToCSV("ionApprox");
    delete ff_ionApprox;
  }
}
