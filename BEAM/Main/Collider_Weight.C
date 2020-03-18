#include "BEAM/Main/Collider_Weight.H"

using namespace BEAM;
using namespace ATOOLS;

Collider_Weight::Collider_Weight(Kinematics_Base * kinematics) :
  Weight_Base(kinematics), m_mode(collidermode::unknown)
{
  if (p_beams[0]->Type()==beamspectrum::monochromatic &&
      p_beams[0]->Type()==beamspectrum::monochromatic)
    m_mode = collidermode::monochromatic;
  else if (p_beams[0]->Type()!=beamspectrum::monochromatic &&
	   p_beams[0]->Type()==beamspectrum::monochromatic)
    m_mode = collidermode::spectral_1;
  else if (p_beams[0]->Type()==beamspectrum::monochromatic &&
	   p_beams[0]->Type()!=beamspectrum::monochromatic)
    m_mode = collidermode::spectral_2;
  else if (p_beams[0]->Type()!=beamspectrum::monochromatic &&
	   p_beams[0]->Type()!=beamspectrum::monochromatic)
    m_mode = collidermode::both_spectral;
  if (m_mode==collidermode::unknown)
    THROW(fatal_error, "Bad settings for collider mode.");
}
  
Collider_Weight::~Collider_Weight() {}

void Collider_Weight::AssignKeys(Integration_Info *const info) {
  m_sprimekey.Assign(m_keyid+std::string("s'"),5,0,info);
  m_ykey.Assign(m_keyid+std::string("y"),3,0,info);
  m_xkey.Assign(m_keyid+std::string("x"),4,0,info);
}

bool Collider_Weight::Calculate(const double & scale) {
  m_weight=1.;
  switch (m_mode) {
  case collidermode::monochromatic:
    return true;
  case collidermode::spectral_1:
    return p_beams[0]->CalculateWeight(m_xkey[2],scale);
  case collidermode::spectral_2:
    return p_beams[1]->CalculateWeight(m_xkey[5],scale);
  case collidermode::both_spectral:
    return (p_beams[0]->CalculateWeight(m_xkey[2],scale) &&
	    p_beams[0]->CalculateWeight(m_xkey[5],scale));
  }
  return false;
}

const double Collider_Weight::operator()(ATOOLS::Flavour * flin) {
  switch (m_mode) {
  case collidermode::monochromatic:
    return 1.;
  case collidermode::spectral_1:
    return p_beams[0]->Weight(flin[0]);
  case collidermode::spectral_2:
    return p_beams[1]->Weight(flin[1]);
  case collidermode::both_spectral:
    return p_beams[0]->Weight(flin[0]) * p_beams[1]->Weight(flin[1]);
  }
  return 0.;
}
