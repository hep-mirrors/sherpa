#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Parameters.H"
#include "BEAM/Spectra/DM_beam.H"
#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/Laser_Backscattering.H"
#include "BEAM/Spectra/Monochromatic.H"
#include "BEAM/Spectra/Pomeron.H"
#include "BEAM/Spectra/Reggeon.H"

using namespace ATOOLS;
using namespace BEAM;

using string = std::string;

std::ostream& BEAM::operator<<(std::ostream& ostr, const beammode bmode)
{
  switch (bmode) {
    case beammode::relic_density: return ostr << "Relic Density calculation";
    case beammode::collider: return ostr << "Collider";
    case beammode::DM_annihilation: return ostr << "Dark Matter annihilation";

    default: break;
  }
  return ostr << "Undefined";
}

std::ostream& BEAM::operator<<(std::ostream& ostr, const collidermode cmode)
{
  switch (cmode) {
    case collidermode::monochromatic: return ostr << "no spectra";
    case collidermode::spectral_1: return ostr << "spectrum for 1";
    case collidermode::spectral_2: return ostr << "spectrum for 2";
    case collidermode::both_spectral: return ostr << "spectra for both";
    default: break;
  }
  return ostr << "Undefined";
}

std::ostream& BEAM::operator<<(std::ostream& ostr, const beamspectrum spect)
{
  switch (spect) {
    case beamspectrum::monochromatic: return ostr << "Monochromatic";
    case beamspectrum::EPA: return ostr << "Equivalent Photons";
    case beamspectrum::Pomeron: return ostr << "Pomeron";
    case beamspectrum::Reggeon: return ostr << "Reggeon";
    case beamspectrum::laser_backscattering:
      return ostr << "Laser Backscattering";
    case beamspectrum::DM: return ostr << "Dark Matter";

    default: break;
  }
  return ostr << "Undefined";
}

Beam_Parameters::Beam_Parameters() : m_settings(Settings::GetMainSettings())
{
  RegisterDefaults();
}

Beam_Parameters::~Beam_Parameters() = default;

void Beam_Parameters::RegisterDefaults()
{
  RegisterDefaultBeams();
  RegisterDarkMatterDefaults();
  RegisterLaserDefaults();
  RegisterEPADefaults();
  RegisterPomeronDefaults();
  RegisterReggeonDefaults();
}

Beam_Base* Beam_Parameters::InitSpectrum(const size_t& num)
{
  switch (GetSpectrum(num)) {
    case beamspectrum::monochromatic: return InitializeMonochromatic(num);
    case beamspectrum::Gaussian:
      THROW(fatal_error, "Gaussian beam spectrum not yet implemented");
    case beamspectrum::laser_backscattering:
      return InitializeLaserBackscattering(num);
    case beamspectrum::simple_Compton: return InitializeSimpleCompton(num);
    case beamspectrum::EPA: return InitializeEPA(num);
    case beamspectrum::Pomeron: return InitializePomeron(num);
    case beamspectrum::Reggeon: return InitializeReggeon(num);
    case beamspectrum::DM: return InitializeDM_beam(num);

    default: break;
  }
  msg_Error() << "Warning in Beam_Initialization::SpecifySpectra :\n"
              << "   No beam spectrum specified for beam " << num + 1
              << "\n   Will initialize monochromatic beam.\n";
  return InitializeMonochromatic(num);
}

Beam_Base* Beam_Parameters::InitializeMonochromatic(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  double  beam_energy =
          Max(m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num],
              beam_particle.Mass());
  double beam_polarization =
          m_settings["BEAM_POLARIZATIONS"].GetTwoVector<double>()[num];
  return new Monochromatic(beam_particle, beam_energy, beam_polarization,
                           1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializeLaserBackscattering(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  if ((beam_particle != Flavour(kf_e)) &&
      (beam_particle != Flavour(kf_e).Bar())) {
    msg_Error() << "Error in Beam_Initialization::SpecifySpectra :\n"
                << "   Tried to initialize Laser_Backscattering for "
                << beam_particle << ".\n";
    return nullptr;
  }
  double beam_energy = m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num];
  double beam_polarization =
          m_settings["BEAM_POLARIZATIONS"].GetTwoVector<double>()[num];
  double laser_energy       = m_settings["E_LASER"].GetTwoVector<double>()[num];
  double laser_polarization = m_settings["P_LASER"].GetTwoVector<double>()[num];
  int    mode               = m_settings["LASER_MODE"].Get<int>();
  bool   angles             = m_settings["LASER_ANGLES"].Get<bool>();
  bool   nonlin             = m_settings["LASER_NONLINEARITY"].Get<bool>();
  return new Laser_Backscattering(beam_particle, beam_energy, beam_polarization,
                                  laser_energy, laser_polarization, mode,
                                  (int) angles, (int) nonlin, 1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializeSimpleCompton(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  if ((beam_particle != Flavour(kf_e)) &&
      (beam_particle != Flavour(kf_e).Bar())) {
    msg_Error() << "Error in Beam_Initialization::SpecifySpectra :\n"
                << "   Tried to initialize Simple_Compton for " << beam_particle
                << ".\n";
    return nullptr;
  }
  double beam_energy = m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num];
  double beam_polarization =
          m_settings["BEAM_POLARIZATIONS"].GetTwoVector<double>()[num];
  double laser_energy       = m_settings["E_LASER"].GetTwoVector<double>()[num];
  double laser_polarization = m_settings["P_LASER"].GetTwoVector<double>()[num];
  return new Laser_Backscattering(beam_particle, beam_energy, beam_polarization,
                                  laser_energy, laser_polarization, -1, 0, 0,
                                  1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializeEPA(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  double beam_energy = m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num];
  if (beam_particle.IsIon()) beam_energy *= beam_particle.GetMassNumber();
  double beam_polarization =
          m_settings["BEAM_POLARIZATIONS"].GetTwoVector<double>()[num];
  return new EPA(beam_particle, beam_energy, beam_polarization, 1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializePomeron(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  if (beam_particle.Kfcode() != kf_p_plus) {
    msg_Error() << "Error in Beam_Initialization::SpecifySpectra:\n"
                << "   Tried to initialize Pomeron for " << beam_particle
                << ".\n"
                << "   This option is not available.\n";
    return nullptr;
  }
  double beam_energy = m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num];
  return new Pomeron(beam_particle, beam_energy, 0., 1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializeReggeon(int num)
{
  Flavour beam_particle = GetFlavour("BEAMS", num);
  if (beam_particle.Kfcode() != kf_p_plus) {
    msg_Error() << "Error in Beam_Initialization::SpecifySpectra:\n"
                << "   Tried to initialize Reggeon for " << beam_particle
                << ".\n"
                << "   This option is not available.\n";
    return nullptr;
  }
  double beam_energy = m_settings["BEAM_ENERGIES"].GetTwoVector<double>()[num];
  return new Reggeon(beam_particle, beam_energy, 0., 1 - 2 * num);
}

Beam_Base* Beam_Parameters::InitializeDM_beam(int num)
{
  Flavour       beam_particle = GetFlavour("BEAMS", num);
  double        temperature   = m_settings["DM_TEMPERATURE"].Get<double>();
  DM_type::code formfactor    = static_cast<DM_type::code>(
          m_settings["DM_ENERGY_DISTRIBUTION"].Get<int>());
  bool relativistic = m_settings["DM_RELATIVISTIC"].Get<bool>();
  return new DM_beam(beam_particle, temperature, formfactor, relativistic,
                     1 - 2 * num);
}

const Flavour Beam_Parameters::GetFlavour(const std::string& tag,
                                          const size_t&      pos)
{
  std::vector<int> beam{m_settings[tag].GetVector<int>()};
  if (beam.size() != 1 && beam.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `BEAMS'.");
  int flav{(pos == 0) ? beam.front() : beam.back()};
  InitializeFlav((kf_code) abs(flav));
  Flavour beam_particle = Flavour((kf_code) abs(flav));
  if (flav < 0) beam_particle = beam_particle.Bar();
  return beam_particle;
}

void Beam_Parameters::RegisterDefaultBeams()
{
  // NOTE: for backwards-compatibility we allow the use of <setting_name>_i
  // with i=1,2 as alternatives for BEAMS, BEAM_SPECTRA, and BEAM_ENERGIES. We
  // do not advertise this in the manual, it's only to make conversion of run
  // cards less error-prone.
  const auto defmode = string("Collider");
  const auto beammode =
          m_settings["BEAM_MODE"].SetDefault(defmode).Get<string>();
  const auto defbeam = 0;
  const auto beam1   = m_settings["BEAM_1"].SetDefault(defbeam).Get<int>();
  const auto beam2   = m_settings["BEAM_2"].SetDefault(defbeam).Get<int>();
  m_settings["BEAMS"].SetDefault({beam1, beam2});

  string     defspectrum{"Monochromatic"};
  const auto spectrum1 =
          m_settings["BEAM_SPECTRUM_1"].SetDefault(defspectrum).Get<string>();
  const auto spectrum2 =
          m_settings["BEAM_SPECTRUM_2"].SetDefault(defspectrum).Get<string>();
  m_settings["BEAM_SPECTRA"].SetDefault({spectrum1, spectrum2});

  const auto defenergy = 0.0;
  const auto energy1 =
          m_settings["BEAM_ENERGY_1"].SetDefault(defenergy).Get<double>();
  const auto energy2 =
          m_settings["BEAM_ENERGY_2"].SetDefault(defenergy).Get<double>();
  m_settings["BEAM_ENERGIES"].SetDefault({energy1, energy2});

  m_settings["BEAM_POLARIZATIONS"].SetDefault({0.0, 0.0});
}

void Beam_Parameters::RegisterDarkMatterDefaults()
{
  m_settings["DM_TEMPERATURE"].SetDefault(1.).Get<double>();
  m_settings["DM_ENERGY_DISTRIBUTION"]
          .SetDefault(int(DM_type::Boltzmann))
          .Get<int>();
  m_settings["DM_beam_weighted"].SetDefault(true).Get<bool>();
  m_settings["DM_RELATIVISTIC"].SetDefault(true).Get<bool>();
  m_settings["RELIC_DENSITY_EMAX"].SetDefault(1.e6).Get<double>();
}

void Beam_Parameters::RegisterEPADefaults()
{
  // EPA defaults are declared in the class itself as the form factor is
  // beam-flavour dependent.
}

void Beam_Parameters::RegisterPomeronDefaults()
{
  m_settings["Pomeron"]["tMax"].SetDefault(1.);
  m_settings["Pomeron"]["xMax"].SetDefault(1.);
  m_settings["Pomeron"]["xMin"].SetDefault(0.);
  // taken from H1 2006 Fit B, hep-ex/0606004
  m_settings["Pomeron"]["B"].SetDefault(5.5);
  m_settings["Pomeron"]["Alpha_intercept"].SetDefault(1.111);
  m_settings["Pomeron"]["Alpha_slope"].SetDefault(0.06);
}

void Beam_Parameters::RegisterReggeonDefaults()
{
  m_settings["Reggeon"]["tMax"].SetDefault(1.);
  m_settings["Reggeon"]["xMax"].SetDefault(1.);
  m_settings["Reggeon"]["xMin"].SetDefault(0.);
  // taken from H1 2006 Fit B, hep-ex/0606004
  m_settings["Reggeon"]["B"].SetDefault(1.6);
  m_settings["Reggeon"]["Alpha_intercept"].SetDefault(0.5);
  m_settings["Reggeon"]["Alpha_slope"].SetDefault(0.3);
  m_settings["Reggeon"]["n"].SetDefault(1.4e-3);
}

void Beam_Parameters::RegisterLaserDefaults()
{
  m_settings["E_LASER"].SetDefault(0.0);
  m_settings["P_LASER"].SetDefault(0.0);
  m_settings["LASER_MODE"].SetDefault(true);
  m_settings["LASER_ANGLES"].SetDefault(false);
  m_settings["LASER_NONLINEARITY"].SetDefault(false);
}

bool Beam_Parameters::SpecifyMode()
{
  string mode = m_settings["BEAM_MODE"].Get<string>();
  if (mode == string("Relic_Density")) m_beammode = beammode::relic_density;
  else if (mode == string("Collider"))
    m_beammode = beammode::collider;
  else if (mode == string("DM_Annihilation"))
    m_beammode = beammode::DM_annihilation;
  else
    m_beammode = beammode::unknown;
  return (m_beammode != beammode::unknown);
}

bool Beam_Parameters::SpecifySpectra()
{
  std::vector<string> beam_spectra{
          m_settings["BEAM_SPECTRA"].GetVector<string>()};
  if (beam_spectra.empty() || beam_spectra.size() > 2)
    THROW(fatal_error, "Specify either one or two values for `BEAM_SPECTRA'.");
  for (short int num = 0; num < 2; num++) {
    string bs{(num == 0) ? beam_spectra.front() : beam_spectra.back()};
    if (bs == "Monochromatic" || bs == "None")
      m_beamspec[num] = beamspectrum::monochromatic;
    else if (bs == "Gaussian")
      m_beamspec[num] = beamspectrum::Gaussian;
    else if (bs == "Laser_Backscattering")
      m_beamspec[num] = beamspectrum::laser_backscattering;
    else if (bs == "Simple_Compton")
      m_beamspec[num] = beamspectrum::simple_Compton;
    else if (bs == "EPA")
      m_beamspec[num] = beamspectrum::EPA;
    else if (bs == "Pomeron")
      m_beamspec[num] = beamspectrum::Pomeron;
    else if (bs == "Reggeon")
      m_beamspec[num] = beamspectrum::Reggeon;
    else if (bs == "DM_beam")
      m_beamspec[num] = beamspectrum::DM;
    else
      m_beamspec[num] = beamspectrum::unknown;
  }
  return (m_beamspec[0] != beamspectrum::unknown &&
          m_beamspec[1] != beamspectrum::unknown);
}

void Beam_Parameters::InitializeFlav(kf_code flav)
{
  constexpr double A = 0.93149410372;
  if (s_kftable.find(flav) == s_kftable.end()) {
    bool initialize_diquarks(false);
    if (flav == kf_p_plus) {
      AddParticle(kf_p_plus, 0.938272, 0.8783, .0, 3, 1, 1, 1, "P+", "P^{+}");
      initialize_diquarks = true;
    } else if (flav == kf_n) {
      AddParticle(kf_n, 0.939566, 0.8783, 7.424e-28, 0, 1, 1, 1, "n", "n");
      initialize_diquarks = true;
    } else if (flav == kf_e) {
      AddParticle(kf_e, 0.000511, .0, .0, -3, 0, 1, 0, 1, 1, 0, "e-", "e+",
                  "e^{-}", "e^{+}");
    } else if (flav == kf_photon) {
      AddParticle(kf_photon, .0, .0, .0, 0, 0, 2, -1, 1, 1, 0, "P", "P", "P",
                  "P");
    } else if (flav == kf_deuterium) {
      AddParticle(kf_deuterium, 2.014 * A, 2.1, 0., 3, 0, true, 1, "deuterium",
                  "H$_2$");
    } else if (flav == kf_helium4) {
      AddParticle(kf_helium4, 4.0026033 * A, 1.9049, 0., 6, 0, true, 1,
                  "helium4", "He$_4$");
    } else if (flav == kf_carbon12) {
      AddParticle(kf_carbon12, 12.0000000 * A, 2.7473, 0., 18, 0, true, 1,
                  "carbon12", "C$_{12}$");
    } else if (flav == kf_calcium40) {
      AddParticle(kf_calcium40, 39.9625912 * A, 4.1039, 0., 60, 0, true, 1,
                  "calcium40", "Ca$_{40}$");
    } else if (flav == kf_silver107) {
      AddParticle(kf_silver107, 106.905097 * A, 5.6970, 0., 141, 0, true, 1,
                  "silver107", "Ag$_{107}$");
    } else if (flav == kf_gold197) {
      AddParticle(kf_gold197, 196.966552 * A, 6.9823, 0., 237, 0, true, 1,
                  "gold197", "Au$_{197}$");
    } else if (flav == kf_lead206) {
      AddParticle(kf_lead206, 205.9744653 * A, 7.0871, 0., 246, 0, true, 1,
                  "lead206", "Pb$_{206}$");
    } else if (flav == kf_lead207) {
      AddParticle(kf_lead207, 206.9758969 * A, 7.0986, 0., 246, 0, true, 1,
                  "lead207", "Pb$_{207}$");
    } else if (flav == kf_lead208) {
      AddParticle(kf_lead208, 207.9766521 * A, 7.1100, 0., 246, 0, true, 1,
                  "lead208", "Pb$_{208}$");
    } else if (flav == kf_uranium238) {
      AddParticle(kf_uranium238, 238.0507900 * A, 7.4366, 0., 276, 0, true, 1,
                  "uranium238", "U$_{238}$");
    } else {
      THROW(fatal_error,
            "You specified a beam particle " + ToString(flav) +
                    "which is not contained in your chosen model. Will abort.");
    }
    if (initialize_diquarks) {
      // Particle_Info(kfc,mass,radius,width,
      // icharge,strong,spin,majorana,on,stable,massive,
      //               idname,antiname,texname,antitexname);
      AddParticle(1103, 0.77133, 0., 0., -2, -3, 2, 0, 1, 1, 1, "dd_1", "dd_1b",
                  "dd_1b", "dd_1b");
      AddParticle(2101, 0.57933, 0., 0., 1, -3, 0, 0, 1, 1, 1, "ud_0", "ud_0b",
                  "ud_0b", "ud_0b");
      AddParticle(2103, 0.77133, 0., 0., 1, -3, 2, 0, 1, 1, 1, "ud_1", "ud_1b",
                  "ud_1b", "ud_1b");
      AddParticle(2203, 0.77133, 0., 0., 4, -3, 2, 0, 1, 1, 1, "uu_1", "uu_1b",
                  "uu_1b", "uu_1b");
      AddParticle(3101, 0.80473, 0., 0., -2, -3, 0, 0, 1, 1, 1, "sd_0", "sd_0b",
                  "sd_0b", "sd_0b");
      AddParticle(3103, 0.92953, 0., 0., -2, -3, 2, 0, 1, 1, 1, "sd_1", "sd_1b",
                  "sd_1b", "sd_1b");
      AddParticle(3201, 0.80473, 0., 0., 1, -3, 0, 0, 1, 1, 1, "su_0", "su_0b",
                  "su_0b", "su_0b");
      AddParticle(3203, 0.92953, 0., 0., 1, -3, 2, 0, 1, 1, 1, "su_1", "su_1b",
                  "su_1b", "su_1b");
      AddParticle(3303, 1.09361, 0., 0., -2, -3, 2, 0, 1, 1, 1, "ss_1", "ss_1b",
                  "ss_1b", "ss_1b");
      AddParticle(4101, 1.96908, 0., 0., 1, -3, 0, 0, 1, 1, 1, "cd_0", "cd_0b",
                  "cd_0b", "cd_0b");
      AddParticle(4103, 2.00808, 0., 0., 1, -3, 2, 0, 1, 1, 1, "cd_1", "cd_1b",
                  "cd_1b", "cd_1b");
      AddParticle(4201, 1.96908, 0., 0., 4, -3, 0, 0, 1, 1, 1, "cu_0", "cu_0b",
                  "cu_0b", "cu_0b");
      AddParticle(4203, 2.00808, 0., 0., 4, -3, 2, 0, 1, 1, 1, "cu_1", "cu_1b",
                  "cu_1b", "cu_1b");
      AddParticle(4301, 2.15432, 0., 0., 1, -3, 0, 0, 1, 1, 1, "cs_0", "cs_0b",
                  "cs_0b", "cs_0b");
      AddParticle(4303, 2.17967, 0., 0., 1, -3, 2, 0, 1, 1, 1, "cs_1", "cs_1b",
                  "cs_1b", "cs_1b");
      AddParticle(4403, 3.27531, 0., 0., 4, -3, 2, 0, 1, 1, 1, "cc_1", "cc_1b",
                  "cc_1b", "cc_1b");
      AddParticle(5101, 5.38897, 0., 0., -2, -3, 0, 0, 1, 1, 1, "bd_0", "bd_0b",
                  "bd_0b", "bd_0b");
      AddParticle(5103, 5.40145, 0., 0., -2, -3, 2, 0, 1, 1, 1, "bd_1", "bd_1b",
                  "bd_1b", "bd_1b");
      AddParticle(5201, 5.38897, 0., 0., 1, -3, 0, 0, 1, 1, 1, "bu_0", "bu_0b",
                  "bu_0b", "bu_0b");
      AddParticle(5203, 5.40145, 0., 0., 1, -3, 2, 0, 1, 1, 1, "bu_1", "bu_1b",
                  "bu_1b", "bu_1b");
      AddParticle(5301, 5.56725, 0., 0., -2, -3, 0, 0, 1, 1, 1, "bs_0", "bs_0b",
                  "bs_0b", "bs_0b");
      AddParticle(5303, 5.57536, 0., 0., -2, -3, 2, 0, 1, 1, 1, "bs_1", "bs_1b",
                  "bs_1b", "bs_1b");
      AddParticle(5401, 6.67143, 0., 0., 1, -3, 0, 0, 1, 1, 1, "bc_0", "bc_0b",
                  "bc_0b", "bc_0b");
      AddParticle(5403, 6.67397, 0., 0., 1, -3, 2, 0, 1, 1, 1, "bc_1", "bc_1b",
                  "bc_1b", "bc_1b");
      AddParticle(5503, 10.07354, 0., 0., -2, -3, 2, 0, 1, 1, 1, "bb_1",
                  "bb_1b", "bb_1b", "bb_1b");
    }
  }
}
