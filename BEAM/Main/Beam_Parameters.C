#include "BEAM/Main/Beam_Parameters.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Spectra/Monochromatic.H"
#include "BEAM/Spectra/Laser_Backscattering.H"
#include "BEAM/Spectra/Fixed_Target.H"
#include "BEAM/Spectra/EPA.H"
#include "BEAM/Spectra/Pomeron.H"
#include "BEAM/Spectra/DM_beam.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;
using namespace std;

std::ostream& BEAM::operator<<(std::ostream& ostr, const beammode bmode) {
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

std::ostream& BEAM::operator<<(std::ostream& ostr, const collidermode cmode) {
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

std::ostream& BEAM::operator<<(std::ostream& ostr, const beamspectrum spect) {
  switch (spect) {
  case beamspectrum::monochromatic:
    return ostr<<"Monochromatic";
  case beamspectrum::EPA:
    return ostr<<"Equivalent Photons";
  case beamspectrum::Pomeron:
    return ostr<<"Pomeron";
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

Beam_Parameters::~Beam_Parameters() = default;

void Beam_Parameters::RegisterDefaults()
{
  RegisterDefaultBeams();
  RegisterDarkMatterDefaults();
  RegisterLaserDefaults();
  RegisterEPADefaults();
  RegisterPomeronDefaults();
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
  case beamspectrum::Pomeron :
    return InitializePomeron(num);
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
    return nullptr;
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
    return nullptr;
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
    return nullptr;
  }
  double beam_energy = (*this)("BEAM_ENERGIES",num);
  if (beam_particle.IsIon()) beam_energy *= beam_particle.GetAtomicNumber();
  double beam_polarization = (*this)("BEAM_POLARIZATIONS",num);
  return new EPA(beam_particle,beam_energy,beam_polarization,1-2*num);
}

Beam_Base * Beam_Parameters::InitializePomeron(int num)
{
  Flavour beam_particle     = GetFlavour("BEAMS",num);
  if (beam_particle.Kfcode()!=kf_p_plus) {
    msg_Error()<<"Error in Beam_Initialization::SpecifySpectra:\n"<<endl
                <<"   Tried to initialize Pomeron for "<<beam_particle<<".\n"
                <<"   This option is not available.\n";
    return nullptr;
  }
  double beam_energy = (*this)("BEAM_ENERGIES",num);
  return new Pomeron(beam_particle,beam_energy,0.,1-2*num);
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

void Beam_Parameters::RegisterEPADefaults() {
  // EPA defaults are declared in the class itself as the form factor is
  // beam-flavour dependent.
}

void Beam_Parameters::RegisterPomeronDefaults() {
  m_settings["Pomeron"]["tMax"].SetDefault(-1.e12);
  // taken from Goharipour:2018yov
  m_settings["Pomeron"]["A"].SetDefault(1.0);
  m_settings["Pomeron"]["B"].SetDefault(7.0);
  m_settings["Pomeron"]["Alpha_intercept"].SetDefault(1.0938);
  m_settings["Pomeron"]["Alpha_slope"].SetDefault(0.);
}

void Beam_Parameters::RegisterLaserDefaults() {
  m_settings["E_LASER"].SetDefault(0.0);
  m_settings["P_LASER"].SetDefault(0.0);
  m_settings["LASER_MODE"].SetDefault(true);
  m_settings["LASER_ANGLES"].SetDefault(false);
  m_settings["LASER_NONLINEARITY"].SetDefault(false);
}

bool Beam_Parameters::SpecifyMode() {
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
  if (beam_spectra.empty() || beam_spectra.size() > 2)
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
    else if (bs == "Pomeron")
      m_beamspec[num] = beamspectrum::Pomeron;
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
    bool initialize_diquarks(false);
    if (flav==kf_p_plus) {
      s_kftable[flav] = new Particle_Info(kf_p_plus,0.938272,0.8783,.0,3,1,1,1,
					  "P+","P^{+}");
      initialize_diquarks=true;
    }
    else if (flav==kf_n) {
      s_kftable[flav] = new Particle_Info(kf_n,0.939566,0.8783,7.424e-28,0,1,1,1,
					  "n","n");
      initialize_diquarks=true;
    }
    else if (flav==kf_e) {
      s_kftable[flav] = new Particle_Info(kf_e,0.000511,.0,.0,-3,0,1,0,1,1,0,
					  "e-","e+","e^{-}","e^{+}");
    }
    else if (flav==kf_photon) {
      s_kftable[flav] = new Particle_Info(22,.0,.0,.0,0,0,2,-1,1,1,0,
					  "P","P","P","P");
    }
    else if (flav==kf_lead208) {
      s_kftable[flav] = new Particle_Info(1000822080, 193.75, 5.5012, 246, 0, 0,
					  "Pb208", "Pb208");
    }
    else if (flav==kf_lead207) {
      s_kftable[flav] = new Particle_Info(1000822070, 192.82, 5.4943, 246, -1, 2,
					  "Pb207", "Pb207");
    }
    else if (flav==kf_lead206) {
      s_kftable[flav] = new Particle_Info(1000822060, 192.82, 5.4902, 246, 0, 2,
					  "Pb206", "Pb206");
    }
    else if (flav==kf_gold197) {
      s_kftable[flav] = new Particle_Info(1000791970, 183.5, 5.4371, 237, 3, 2,
					  "Au197", "Au197");
    }
    else if (flav==kf_calcium40) {
      s_kftable[flav] = new Particle_Info(1000200400, 37.26, 3.4776, 60, 0, 2,
					  "Ca40", "Ca40");
    }
    else {
      THROW(fatal_error,"You specified a beam particle "+ToString(flav)+
            "which is not contained in your chosen model. Will abort.");
    }
    if (initialize_diquarks) {
      s_kftable[1103] = new Particle_Info(1103,0.77133,0.,0,-2,-3,2,0,0,1,1,
					  "dd_1","dd_1b","dd_1b","dd_1b");
      s_kftable[2101] = new Particle_Info(2101,0.57933,0.,0,1,-3,0,0,0,1,1,
					  "ud_0","ud_0b","ud_0b","ud_0b");
      s_kftable[2103] = new Particle_Info(2103,0.77133,0.,0,1,-3,2,0,0,1,1,
					  "ud_1","ud_1b","ud_1b","ud_1b");
      s_kftable[2203] = new Particle_Info(2203,0.77133,0.,0,4,-3,2,0,0,1,1,
					  "uu_1","uu_1b","uu_1b","uu_1b");
      s_kftable[3101] = new Particle_Info(3101,0.80473,0.,0,-2,-3,0,0,0,1,1,
					  "sd_0","sd_0b","sd_0b","sd_0b");
      s_kftable[3103] = new Particle_Info(3103,0.92953,0.,0,-2,-3,2,0,0,1,1,
					  "sd_1","sd_1b","sd_1b","sd_1b");
      s_kftable[3201] = new Particle_Info(3201,0.80473,0.,0,1,-3,0,0,0,1,1,
					  "su_0","su_0b","su_0b","su_0b");
      s_kftable[3203] = new Particle_Info(3203,0.92953,0.,0,1,-3,2,0,0,1,1,
					  "su_1","su_1b","su_1b","su_1b");
      s_kftable[3303] = new Particle_Info(3303,1.09361,0.,0,-2,-3,2,0,0,1,1,
					  "ss_1","ss_1b","ss_1b","ss_1b");
      s_kftable[4101] = new Particle_Info(4101,1.96908,0.,0,1,-3,0,0,0,1,1,
					  "cd_0","cd_0b","cd_0b","cd_0b");
      s_kftable[4103] = new Particle_Info(4103,2.00808,0.,0,1,-3,2,0,0,1,1,
					  "cd_1","cd_1b","cd_1b","cd_1b");
      s_kftable[4201] = new Particle_Info(4201,1.96908,0.,0,4,-3,0,0,0,1,1,
					  "cu_0","cu_0b","cu_0b","cu_0b");
      s_kftable[4203] = new Particle_Info(4203,2.00808,0.,0,4,-3,2,0,0,1,1,
					  "cu_1","cu_1b","cu_1b","cu_1b");
      s_kftable[4301] = new Particle_Info(4301,2.15432,0.,0,1,-3,0,0,0,1,1,
					  "cs_0","cs_0b","cs_0b","cs_0b");
      s_kftable[4303] = new Particle_Info(4303,2.17967,0.,0,1,-3,2,0,0,1,1,
					  "cs_1","cs_1b","cs_1b","cs_1b");
      s_kftable[4403] = new Particle_Info(4403,3.27531,0.,0,4,-3,2,0,0,1,1,
					  "cc_1","cc_1b","cc_1b","cc_1b");
      s_kftable[5101] = new Particle_Info(5101,5.38897,0.,0,-2,-3,0,0,0,1,1,
					  "bd_0","bd_0b","bd_0b","bd_0b");
      s_kftable[5103] = new Particle_Info(5103,5.40145,0.,0,-2,-3,2,0,0,1,1,
					  "bd_1","bd_1b","bd_1b","bd_1b");
      s_kftable[5201] = new Particle_Info(5201,5.38897,0.,0,1,-3,0,0,0,1,1,
					  "bu_0","bu_0b","bu_0b","bu_0b");
      s_kftable[5203] = new Particle_Info(5203,5.40145,0.,0,1,-3,2,0,0,1,1,
					  "bu_1","bu_1b","bu_1b","bu_1b");
      s_kftable[5301] = new Particle_Info(5301,5.56725,0.,0,-2,-3,0,0,0,1,1,
					  "bs_0","bs_0b","bs_0b","bs_0b");
      s_kftable[5303] = new Particle_Info(5303,5.57536,0.,0,-2,-3,2,0,0,1,1,
					  "bs_1","bs_1b","bs_1b","bs_1b");
      s_kftable[5401] = new Particle_Info(5401,6.67143,0.,0,1,-3,0,0,0,1,1,
					  "bc_0","bc_0b","bc_0b","bc_0b");
      s_kftable[5403] = new Particle_Info(5403,6.67397,0.,0,1,-3,2,0,0,1,1,
					  "bc_1","bc_1b","bc_1b","bc_1b");
      s_kftable[5503] = new Particle_Info(5503,10.07354,0.,0,-2,-3,2,0,0,1,1,
					  "bb_1","bb_1b","bb_1b","bb_1b");
    }
  }
}

