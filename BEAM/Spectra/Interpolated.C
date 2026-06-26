#include "BEAM/Spectra/Interpolated.H"

using namespace BEAM;
using namespace ATOOLS;

Interpolated_Neutrinos::
Interpolated_Neutrinos(const ATOOLS::Flavour & flav,const double & energy,const int dir) :
  Beam_Base(beamspectrum::neutrinos_from_protons,flav,energy,0,dir,0),
  p_table(nullptr)
{
  m_Nbunches = 2;
  m_bunches.resize(m_Nbunches);
  m_bunches[0] = flav;
  m_bunches[1] = flav;
  m_vecouts.resize(m_Nbunches);
  m_vecouts[0] = Vec4D(m_energy, 0., 0., m_dir * m_energy);
  m_vecouts[1] = Vec4D(0., 0., 0., 0.);
  m_on         = true;
  Initialise();
}

void Interpolated_Neutrinos::Initialise() {
  const auto& s    = Settings::GetMainSettings()["NEUTRINO_BEAM"];
  int pdgid        = s["Neutrino"].SetDefault(14).Get<int>();
  Flavour neutrino = Flavour(pdgid);
  m_bunches[0]     = neutrino;
  std::string path = s["Path"].SetDefault("./Faserv2").Get<std::string>();
  std::string file = (s["File"].SetDefault("FASERv2_"+ToString(pdgid)+".txt").
		      Get<std::string>());
  std::string fullpath = path+"/"+file; 
  std::ifstream f(fullpath.c_str());
  if (!f.good()) 
    THROW(fatal_error,"File ["+fullpath+"] does not exist.");
  p_table = new OneDim_Flexible_Table(fullpath);
  double  lumi     = s["LUMINOSITY"].SetDefault(3000.).Get<double>();
}

Beam_Base * Interpolated_Neutrinos::Copy() {
  return nullptr;
}

Flavour Interpolated_Neutrinos::Remnant() {
  return Flavour(kf_none);
}

bool Interpolated_Neutrinos::
CalculateWeight(const double x,const double Q2) {
  return true;
}
