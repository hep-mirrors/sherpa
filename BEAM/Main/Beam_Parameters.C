#include "BEAM/Main/Beam_Parameters.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Spectra/Monochromatic.H"
#include "BEAM/Spectra/Laser_Backscattering.H"
#include "BEAM/Spectra/Fixed_Target.H"
#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/DM_beam.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;
using namespace std;

std::ostream& BEAM::operator<<(std::ostream& ostr, const beammode::code bmode) {
  switch (bmode) {
  case beammode::relic_density:
    return ostr<<"Relic Density calculation";
  case beammode::collider:
    return ostr<<"Collider";
  case beammode::DM_annihilation:
    return ostr<<"Dark Matter annihilation";
  case beammode::Fixed_Target:
    return ostr<<"Fixed Target annihilation";
  default:
    break;
  }
  return ostr<<"Undefined";
}

std::ostream& BEAM::operator<<(std::ostream& ostr, const collidermode::code cmode) {
  switch (cmode) {
  case collidermode::monochromatic:
    return ostr<<"no spectra";
  case collidermode::spectral_1:
    return ostr<<"spectrum for 1";
  case collidermode::spectral_2:
    return ostr<<"spectrum for 2";
  case collidermode::both_spectral:
    return ostr<<"spectra for both";
  default:
    break;
  }
  return ostr<<"Undefined";
}

std::ostream& BEAM::operator<<(std::ostream& ostr, const beamspectrum::code spect) {
  switch (spect) {
  case beamspectrum::monochromatic:
    return ostr<<"Monochromatic";
  case beamspectrum::EPA:
    return ostr<<"Equivalent Photons";
  case beamspectrum::laser_backscattering:
    return ostr<<"Laser Backscattering";
  case beamspectrum::DM:
    return ostr<<"Dark Matter";
  case beamspectrum::Fixed_Target:
    return ostr<<"Fixed Target";
  default:
    break;
  }
  return ostr<<"Undefined";
}

Beam_Parameters::Beam_Parameters() :
  m_settings(Settings::GetMainSettings())
{
  RegisterDefaults();
}

Beam_Parameters::~Beam_Parameters() {}

void Beam_Parameters::RegisterDefaults()
{
  RegisterDefaultBeams();
  RegisterDarkMatterDefaults();
  RegisterLaserDefaults();
  RegisterEPADefaults();
}

Beam_Base * Beam_Parameters::InitSpectrum(const size_t & num) {
  switch (GetSpectrum(num)) {
  case beamspectrum::monochromatic :
    return InitializeMonochromatic(num);
  case beamspectrum::Gaussian :
    THROW(fatal_error, "Gaussian beam spectrum not yet implemented");
  case beamspectrum::laser_backscattering :
    return InitializeLaserBackscattering(num);
  case beamspectrum::simple_Compton :
    return InitializeSimpleCompton(num);
  case beamspectrum::EPA :
    return InitializeEPA(num);
  case beamspectrum::DM :
    return InitializeDM_beam(num);
  case beamspectrum::Fixed_Target :
    return InitializeFixed_Target(num);
  default :
    break;
  }
  msg_Error()<<"Warning in Beam_Initialization::SpecifySpectra :"<<endl
	     <<"   No beam spectrum specified for beam "<<num+1<<endl
	     <<"   Will initialize monochromatic beam."<<endl;
  return InitializeMonochromatic(num);
}

Beam_Base * Beam_Parameters::InitializeMonochromatic(int num)
{
  Flavour beam_particle     = GetFlavour("BEAMS",num);
  double beam_energy        = Max((*this)("BEAM_ENERGIES",num), beam_particle.Mass());
  double beam_polarization  = (*this)("BEAM_POLARIZATIONS",num);
  return new Monochromatic(beam_particle,beam_energy,beam_polarization,1-2*num);
}

Beam_Base * Beam_Parameters::InitializeLaserBackscattering(int num)
{
  Flavour beam_particle     = GetFlavour("BEAMS",num);
  if ((beam_particle!=Flavour(kf_e)) && (beam_particle!=Flavour(kf_e).Bar())) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra :\n"
	       <<"   Tried to initialize Laser_Backscattering for "
	       <<beam_particle<<".\n";
    return NULL;
  }
  double beam_energy        = (*this)("BEAM_ENERGIES",num); 
  double beam_polarization  = (*this)("BEAM_POLARIZATIONS",num);
  double laser_energy       = (*this)("E_LASER",num);
  double laser_polarization = (*this)("P_LASER",num);
  int  mode   = (*this)("LASER_MODE");
  bool angles = (*this)("LASER_ANGLES");
  bool nonlin = (*this)("LASER_NONLINEARITY");
  return new Laser_Backscattering(beam_particle,
				  beam_energy,beam_polarization,
				  laser_energy,laser_polarization,
				  mode,(int)angles,(int)nonlin,1-2*num);
}

Beam_Base * Beam_Parameters::InitializeSimpleCompton(int num)
{
  Flavour beam_particle     = GetFlavour("BEAMS",num);
  if ((beam_particle!=Flavour(kf_e)) && (beam_particle!=Flavour(kf_e).Bar())) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra :\n"
	       <<"   Tried to initialize Simple_Compton for "
	       <<beam_particle<<".\n";
    return NULL;
  }
  double beam_energy        = (*this)("BEAM_ENERGIES",num); 
  double beam_polarization  = (*this)("BEAM_POLARIZATIONS",num);
  double laser_energy       = (*this)("E_LASER",num);
  double laser_polarization = (*this)("P_LASER",num);
  return new Laser_Backscattering(beam_particle,
				  beam_energy,beam_polarization,
				  laser_energy,laser_polarization,
				  -1,0,0,1-2*num);
}

Beam_Base * Beam_Parameters::InitializeEPA(int num)
{
  Flavour beam_particle     = GetFlavour("BEAMS",num);
  if (beam_particle.Kfcode()!=kf_p_plus &&
      beam_particle.Kfcode()!=kf_e &&
      !beam_particle.IsIon()) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra:\n"<<endl
               <<"   Tried to initialize EPA for "<<beam_particle<<".\n"
	       <<"   This option is not available (yet).\n";
    return NULL;
  }
  double beam_energy = (*this)("BEAM_ENERGIES",num); 
  if (beam_particle.IsIon()) beam_energy *= beam_particle.GetAtomicNumber();
  double beam_polarization = (*this)("BEAM_POLARIZATIONS",num);
  return new EPA(beam_particle,beam_energy,beam_polarization,1-2*num);
}

Beam_Base * Beam_Parameters::InitializeDM_beam(int num)
{
  Flavour beam_particle    = GetFlavour("BEAMS",num);
  double  temperature      = (*this)("DM_TEMPERATURE");
  DM_type::code formfactor = DM_type::code(Switch("DM_ENERGY_DISTRIBUTION"));
  bool    relativistic     = On("DM_RELATIVISTIC");
  return new DM_beam(beam_particle,temperature,
		     formfactor,relativistic,1-2*num);
}

Beam_Base * Beam_Parameters::InitializeFixed_Target(int num)
{
  double beam_energy        = (*this)("BEAM_ENERGIES",num); 
  double beam_polarization  = (*this)("BEAM_POLARIZATIONS",num);
  Flavour beam_particle    = GetFlavour("BEAMS",num);
  return new Fixed_Target(beam_particle,beam_energy,beam_polarization,1-2*num);
}

const Flavour Beam_Parameters::GetFlavour(const std::string & tag,const size_t & pos) {
  vector<int> beam{ m_settings[tag].GetVector<int>() };
  if (beam.size() != 1 && beam.size() != 2)
    THROW(fatal_error, "Specify either one or two values for `BEAMS'.");
  int flav{ (pos == 0) ? beam.front() : beam.back() };
  InitializeFlav((kf_code)abs(flav));
  Flavour beam_particle     = Flavour((kf_code)abs(flav));
  if (flav<0) beam_particle = beam_particle.Bar();
  return beam_particle;
}

const std::string Beam_Parameters::String(const string & tag,const int & pos) const {
  if (pos<0) return m_settings[tag].Get<string>();
  vector<string> params{ m_settings[tag].GetVector<string>() };
  if (pos>1 || pos>params.size()-1) {
    string message = string("Parameter number mismatch for tag = ")+tag+string(" at pos = ")+ToString(pos);
    THROW(fatal_error, message);
  }
  return (pos==0)?params.front():params.back();      
}

const double Beam_Parameters::operator()(const string & tag,const int & pos) const {
  if (pos<0) return m_settings[tag].Get<double>();
  vector<double> params{ m_settings[tag].GetVector<double>() };
  if (!(tag=="BEAM_ENERGIES" || tag=="BEAM_POLARIZATIONS") &&
      (pos>1 || pos>params.size()-1)) {
    string message = string("Parameter number mismatch for tag = ")+tag+string(" at pos = ")+ToString(pos);
    THROW(fatal_error, message);
  }
  return (pos==0)?params.front():params.back();      
}

const int Beam_Parameters::Switch(const string & tag,const int & pos) const {
  if (pos<0) return m_settings[tag].Get<int>();
  vector<int> params{ m_settings[tag].GetVector<int>() };
  if (pos>1 || pos>params.size()-1)  {
    string message = string("Parameter number mismatch for tag = ")+tag+string(" at pos = ")+ToString(pos);
    THROW(fatal_error, message);
  }
  return (pos==0)?params.front():params.back();      
}

const bool Beam_Parameters::On(const string & tag) const {
  return m_settings[tag].Get<bool>();
}


void Beam_Parameters::RegisterDefaultBeams() {
  // NOTE: for backwards-compatibility we allow the use of <setting_name>_i
  // with i=1,2 as alternatives for BEAMS, BEAM_SPECTRA, and BEAM_ENERGIES. We
  // do not advertise this in the manual, it's only to make conversion of run
  // cards less error-prone.
  const auto defmode  = string("Collider"); 
  const auto beammode = m_settings["BEAM_MODE"].SetDefault(defmode).Get<string>();
  const auto defbeam = 0;
  const auto beam1 = m_settings["BEAM_1"].SetDefault(defbeam).Get<int>();
  const auto beam2 = m_settings["BEAM_2"].SetDefault(defbeam).Get<int>();
  m_settings["BEAMS"].SetDefault({beam1, beam2});

  string defspectrum {"Monochromatic"};
  const auto spectrum1
    = m_settings["BEAM_SPECTRUM_1"].SetDefault(defspectrum).Get<string>();
  const auto spectrum2
    = m_settings["BEAM_SPECTRUM_2"].SetDefault(defspectrum).Get<string>();
  m_settings["BEAM_SPECTRA"].SetDefault({spectrum1, spectrum2});

  const auto defenergy = 0.0;
  const auto energy1 = m_settings["BEAM_ENERGY_1"].SetDefault(defenergy).Get<double>();
  const auto energy2 = m_settings["BEAM_ENERGY_2"].SetDefault(defenergy).Get<double>();
  m_settings["BEAM_ENERGIES"].SetDefault({energy1, energy2});

  m_settings["BEAM_POLARIZATIONS"].SetDefault({0.0, 0.0});
}

void Beam_Parameters::RegisterDarkMatterDefaults() {
  m_settings["DM_TEMPERATURE"].SetDefault(1.).Get<double>();
  m_settings["DM_ENERGY_DISTRIBUTION"].SetDefault(int(DM_type::Boltzmann)).Get<int>();
  m_settings["DM_beam_weighted"].SetDefault(true).Get<bool>();
  m_settings["DM_RELATIVISTIC"].SetDefault(true).Get<bool>();
  m_settings["RELIC_DENSITY_EMAX"].SetDefault(1.e6).Get<double>();
}

void Beam_Parameters::RegisterEPADefaults() {}

void Beam_Parameters::RegisterLaserDefaults() {
  m_settings["E_LASER"].SetDefault(0.0);
  m_settings["P_LASER"].SetDefault(0.0);
  m_settings["LASER_MODE"].SetDefault(true);
  m_settings["LASER_ANGLES"].SetDefault(false);
  m_settings["LASER_NONLINEARITY"].SetDefault(false);
}

bool Beam_Parameters::SpecifyMode() {
  bool okay(true);
  string mode = m_settings["BEAM_MODE"].Get<string>();
  if      (mode==string("Relic_Density"))
    m_beammode = beammode::relic_density; 
  else if (mode==string("Collider"))
    m_beammode = beammode::collider;
  else if (mode==string("DM_Annihilation"))
    m_beammode = beammode::DM_annihilation;
  else if (mode==string("Fixed_Target"))
    m_beammode = beammode::Fixed_Target;
  else
    m_beammode = beammode::unknown;
  return (m_beammode!=beammode::unknown);
}

bool Beam_Parameters::SpecifySpectra() {
  vector<string> beam_spectra{ m_settings["BEAM_SPECTRA"].GetVector<string>() };
  if (beam_spectra.size() == 0 || beam_spectra.size() > 2)
    THROW(fatal_error, "Specify either one or two values for `BEAM_SPECTRA'.");
  for (short int num=0;num<2;num++) {
    string bs{ (num == 0) ? beam_spectra.front() : beam_spectra.back() };
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
    else if (bs == "DM_beam")
      m_beamspec[num] = beamspectrum::DM;
    else if (bs == "Fixed_Target")
      m_beamspec[num] = beamspectrum::Fixed_Target;
    else
      m_beamspec[num] = beamspectrum::unknown;
  }
  return (m_beamspec[0]!=beamspectrum::unknown && m_beamspec[1]!=beamspectrum::unknown);
}

void Beam_Parameters::InitializeFlav(kf_code flav) {
  if (s_kftable.find(flav)==s_kftable.end()) {
    if (flav==kf_p_plus) {
      AddParticle(kf_p_plus, 0.938272, 0.8783, .0, 3, 1, 1, 1, "P+", "P^{+}");
    }
    else if (flav==kf_n) {
      AddParticle(kf_n, 0.939566, 0.8783, 7.424e-28, 0, 1, 1, 1, "n", "n");
    }
    else if (flav==kf_e) {
      AddParticle(kf_e, 0.000511, .0, .0, -3, 0, 1, 0, 1, 1, 0, "e-", "e+", "e^{-}", "e^{+}");
    }
    else if (flav==kf_photon) {
      AddParticle(22, .0, .0, .0, 0, 0, 2, -1, 1, 1, 0, "P", "P", "P", "P");
    }
    else if (flav==kf_lead208) {
      AddParticle(1000822080, 193.75, 5.5012, 246, 0, 0, "Pb208", "Pb208");
    }
    else if (flav==kf_lead207) {
      AddParticle(1000822070, 192.82, 5.4943, 246, -1, 2, "Pb207", "Pb207");
    }
    else if (flav==kf_lead206) {
      AddParticle(1000822060, 192.82, 5.4902, 246, 0, 2, "Pb206", "Pb206");
    }
    else if (flav==kf_gold197) {
      AddParticle(1000791970, 183.5, 5.4371, 237, 3, 2, "Au197", "Au197");
    }
    else if (flav==kf_calcium40) {
      AddParticle(1000200400, 37.26, 3.4776, 60, 0, 2, "Ca40", "Ca40");
    }
    else {
      THROW(fatal_error,"You specified a beam particle "+ToString(flav)+
            "which is not contained in your chosen model. Will abort.");
    }
  }
}
