#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "BEAM/Main/Beam_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "BEAM/Main/Collider_Kinematics.H"
#include "BEAM/Main/Collider_Weight.H"
#include "BEAM/Main/DM_Annihilation_Kinematics.H"
#include "BEAM/Main/DM_Annihilation_Weight.H"
#include "BEAM/Main/RelicDensity_Kinematics.H"
#include "BEAM/Main/RelicDensity_Weight.H"

using namespace ATOOLS;
using namespace BEAM;
using namespace std;

Beam_Spectra_Handler::Beam_Spectra_Handler()
    : p_kinematics(nullptr), p_weight(nullptr), m_beammode(beammode::collider),
      m_collidermode(collidermode::monochromatic), m_mode(0),
      m_polarisation(0) {
  msg_Info() << "Initializing beam spectra ...\n";
  for (size_t i = 0; i < 2; i++)
    p_BeamBase[i] = nullptr;
  // simple check for remotely sensible beam parameters
  if (!m_parameters.SpecifyMode() || !m_parameters.SpecifySpectra())
    THROW(fatal_error, "Bad parameters in BEAM_SPECTRA_HANDLER.");
  m_beammode = m_parameters.GetMode();
  // initialise the beams, kinematics & weights
  if (!InitTheBeams() || !InitTheKinematics() || !InitTheWeight())
    THROW(fatal_error, "Bad spectra in BEAM_SPECTRA_HANDLER.");
  m_on = p_kinematics->On();

  msg_Info() << "  Type: " << m_type << std::endl
             << "  Beam 1: " << p_BeamBase[0]->Beam()
             << " (enabled = " << p_BeamBase[0]->On() << ", "
             << "momentum = " << p_BeamBase[0]->InMomentum() << ")" << std::endl
             << "  Beam 2: " << p_BeamBase[1]->Beam()
             << " (enabled = " << p_BeamBase[0]->On() << ", "
             << "momentum = " << p_BeamBase[1]->InMomentum() << ")" << std::endl;
}

Beam_Spectra_Handler::~Beam_Spectra_Handler() {
  for (short int i = 0; i < 2; i++) {
    if (p_BeamBase[i]) {
      delete p_BeamBase[i];
      p_BeamBase[i] = nullptr;
    }
  }
  if (p_weight!=NULL)     { delete p_weight;     p_weight     = NULL; }
  if (p_kinematics!=NULL) { delete p_kinematics; p_kinematics = NULL; }
  }

bool Beam_Spectra_Handler::InitTheBeams() {
  for (short int i = 0; i < 2; i++) {
    p_BeamBase[i] = m_parameters.InitSpectrum(i);
    if (p_BeamBase[i] == nullptr)
      return false;
    if (p_BeamBase[i]->On())
      m_mode += i + 1;
    if (p_BeamBase[i]->PolarisationOn())
      m_polarisation += i + 1;
  }
  switch (m_mode) {
  case 1:
    m_collidermode = collidermode::spectral_1;
    break;
  case 2:
    m_collidermode = collidermode::spectral_2;
    break;
  case 3:
    m_collidermode = collidermode::both_spectral;
    break;
  default:
    break;
  }
  rpa->gen.SetBeam1(p_BeamBase[0]->Beam());
  rpa->gen.SetBeam2(p_BeamBase[1]->Beam());
  rpa->gen.SetPBeam(0, p_BeamBase[0]->InMomentum());
  rpa->gen.SetPBeam(1, p_BeamBase[1]->InMomentum());
  rpa->gen.SetPBunch(0, p_BeamBase[0]->OutMomentum());
  rpa->gen.SetPBunch(1, p_BeamBase[1]->OutMomentum());
  double ecms = (p_BeamBase[0]->InMomentum()+p_BeamBase[1]->InMomentum()).Abs();
  rpa->gen.SetEcms(ecms);
  Settings::GetMainSettings().AddGlobalTag("E_CMS", ToString(ecms));
  return true;
}

bool Beam_Spectra_Handler::InitTheKinematics() {
  switch (m_beammode) {
  case beammode::relic_density:
      m_type       = std::string("Relic Density");
      p_kinematics = new RelicDensity_Kinematics(p_BeamBase);
      break;
  case beammode::collider:
    m_type       = std::string("Collider Setup");
    p_kinematics = new Collider_Kinematics(p_BeamBase);
    break;
  case beammode::DM_annihilation:
    m_type       = std::string("DM Annihilation");
    p_kinematics = new DM_Annihilation_Kinematics(p_BeamBase);
    break;
  case beammode::unknown:
    break;
  default:
    break;
  }
  return (p_kinematics != nullptr);
}

bool Beam_Spectra_Handler::InitTheWeight() {
  switch (m_beammode) {
  case beammode::relic_density:
    p_weight = new RelicDensity_Weight(p_kinematics);
    break;
  case beammode::collider:
    p_weight = new Collider_Weight(p_kinematics);
    break;
  case beammode::DM_annihilation:
    p_weight = new DM_Annihilation_Weight(p_kinematics);
    break;
  case beammode::unknown:
    break;
  default:
    break;
  }
  return (p_weight != nullptr);
}

void Beam_Spectra_Handler::FixPositions() {
  for (short int beam=0;beam<2;beam++) p_BeamBase[beam]->FixPosition();
}

// TODO: Improve this handling for rescattering etc.
bool Beam_Spectra_Handler::CheckConsistency(const ATOOLS::Flavour *_beams,
                                            const ATOOLS::Flavour *_bunches) {
  for (int i = 0; i < 2; i++) {
    if (_beams[i]   != GetBeam(i)->Beam() ||
  _bunches[i] != GetBeam(i)->Bunch() ) return false;
  }
  return true;
}

bool Beam_Spectra_Handler::CheckConsistency(const ATOOLS::Flavour *_bunches) {
  for (int i = 0; i < 2; i++) {
    if (_bunches[i] != GetBeam(i)->Bunch()) return false;
  }
  return true;
}
