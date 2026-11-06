#include "BEAM/Main/Beam_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"

#include <cmath>

using namespace BEAM;

namespace {

  Transverse_Profile MakeTransverseProfile(beamspectrum type,
                                          const ATOOLS::Flavour& beam, int dir)
  {
    // EPA supplies an x-dependent photon-emission profile of its own. Do not
    // apply the generic parton profile or its settings to that spectrum.
    if (type == beamspectrum::EPA)
      return Transverse_Profile(Transverse_Profile_Type::Gaussian, 0.);

    auto s = ATOOLS::Settings::GetMainSettings()["BEAM_TRANSVERSE_PROFILE"];
    const size_t index = dir > 0 ? 0 : 1;
    const auto models = s["TYPE"].SetDefault(Transverse_Profile_Type::Gaussian)
                            .GetTwoVector<Transverse_Profile_Type>();
    s.DeclareVectorSettingsWithEmptyDefault({"RMS_RADIUS"});
    // Allow override of radius
    const bool has_rms = s["RMS_RADIUS"].IsSetExplicitly() ||
                        !s["RMS_RADIUS"].GetVector<double>().empty();
    const double rms_fm = has_rms
        ? s["RMS_RADIUS"].GetTwoVector<double>()[index]
        : std::sqrt(2. / 3.) * beam.Radius();
    return Transverse_Profile(models[index], rms_fm / ATOOLS::rpa->hBar_c());
  }

} // namespace

Beam_Base::Beam_Base(beamspectrum _type, const ATOOLS::Flavour& _beam,
                     const double _energy, const double _polarisation,
                     const int _dir)
    : m_type(_type), m_beam(_beam),
      m_position(ATOOLS::Vec4D(0., 0., 0., 0.)), m_Nbunches(1), m_dir(_dir),
      m_energy(_energy), m_polarisation(_polarisation), m_x(1.), m_Q2(0.),
      m_weight(1.), m_on(false),
      m_transverse_profile(MakeTransverseProfile(_type, _beam, _dir)), m_eran(0.)
{
  if (_energy == 0.) { // Fixed Target mode
    m_lab = ATOOLS::Vec4D(m_beam.Mass(true), 0., 0., 0.);
  } else {
    double mass = m_beam.Mass();
    // For heavy ions we internally only ever treat the per-nucleon momentum;
    // hence we scale down the mass from the nucleus mass to the nucleon mass.
    if (m_beam.IsIon())
      mass /= m_beam.GetMassNumber();
    double disc = 1.0 - ATOOLS::sqr(mass / m_energy);
    if (disc < 0.) {
      msg_Error() << "Error in Beam_Base :" << m_type << std::endl
                  << "   Mismatch of energy and mass of beam particle : "
                  << m_beam << " / " << m_energy << std::endl
                  << "   Will lead to termination of program." << std::endl;
      ATOOLS::Abort();
    }
    double pz = m_dir * m_energy * sqrt(disc);
    m_lab     = ATOOLS::Vec4D(m_energy, 0., 0., pz);
  }
  m_bunches.resize(m_Nbunches, m_beam);
  m_vecouts.resize(m_Nbunches, m_lab);

  msg_Info() << "  Beam " << (m_dir > 0 ? 1 : 2)
             << " transverse profile: ";
  if (m_type == beamspectrum::EPA) {
    msg_Info() << "EPA photon emission, no separate profile.\n";
  } else {
    msg_Info() << m_transverse_profile.Type()
               << ", transverse RMS radius = "
               << m_transverse_profile.RMSRadius() * ATOOLS::rpa->hBar_c()
               << " fm"
               << (m_transverse_profile.RMSRadius() == 0. ? " (pointlike)" : "")
               << "; used only by overlap geometry.\n";
  }
}
