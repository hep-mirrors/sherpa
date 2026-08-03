#include "BEAM/Spectra/Interpolated.H"

using namespace BEAM;
using namespace ATOOLS;

Interpolated_Neutrinos::
Interpolated_Neutrinos(const ATOOLS::Flavour & flav,const double & energy,const int dir) :
  Beam_Base(beamspectrum::neutrinos_from_protons,flav,energy,0,dir,0),
  p_table(nullptr), m_xmin(0.), m_xmax(1.), m_lumiscale(1.)
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
  const auto& s = Settings::GetMainSettings()["NEUTRINO_BEAM"];

  std::string fullpath = s["File"].SetDefault("").Get<std::string>();
  if (fullpath.empty())
    THROW(fatal_error,"NEUTRINO_BEAM: File must be set to the flux file's path.");
  std::ifstream f(fullpath.c_str());
  if (!f.good())
    THROW(fatal_error,"File ["+fullpath+"] does not exist.");
  p_table = new OneDim_Flexible_Table(fullpath);

  // Luminosity [fb^-1] the flux file was produced for.
  double bakedinlumi = s["FileLuminosity"].SetDefault(-1.).Get<double>();
  if (bakedinlumi<=0.)
    THROW(fatal_error,"NEUTRINO_BEAM: FileLuminosity must be set to the "
	  "integrated luminosity [fb^-1] the flux file was produced for.");
  double targetlumi = s["LUMINOSITY"].SetDefault(bakedinlumi).Get<double>();
  m_lumiscale = targetlumi/bakedinlumi;

  double tablemin = p_table->GetX()->front();
  double tablemax = p_table->GetX()->back();
  m_xmin = Max(tablemin/m_energy, 0.);
  m_xmax = Min(tablemax/m_energy, 1.);
}

Beam_Base * Interpolated_Neutrinos::Copy() {
  return new Interpolated_Neutrinos(m_beam,m_energy,m_dir);
}

Flavour Interpolated_Neutrinos::Remnant() {
  return Flavour(kf_none);
}

bool Interpolated_Neutrinos::
CalculateWeight(const double x,const double Q2) {
  m_x  = x;
  m_Q2 = Q2;
  double enu = x*m_energy;
  if (enu<p_table->GetX()->front() || enu>p_table->GetX()->back()) {
    m_weight = 0.;
    return true;
  }
  // Cancels the dx=dEnu/m_energy Jacobian from PHASIC's integration over x.
  m_weight = m_lumiscale * (*p_table)(enu) * m_energy;
  if (m_weight<0. || IsNan(m_weight)) m_weight = 0.;
  return true;
}
