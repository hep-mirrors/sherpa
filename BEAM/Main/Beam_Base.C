#include "BEAM/Main/Beam_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Settings.H"
#include "BEAM/Main/Woods_Saxon.H"

#include <cmath>

using namespace BEAM;

namespace {

  Transverse_Profile MakeTransverseProfile(beamspectrum type,
                                          const ATOOLS::Flavour& beam, int dir)
  {
    const size_t index = dir > 0 ? 0 : 1;
    const double hbarc = ATOOLS::rpa->hBar_c();
    msg_Info() << "  Beam " << index + 1 << " transverse profile: ";
    // EPA supplies an x-dependent photon-emission profile of its own. Do not
    // apply the generic parton profile or its settings to that spectrum.
    if (type == beamspectrum::EPA) {
      msg_Info() << "EPA photon emission, no separate profile.\n";
      return Transverse_Profile(Transverse_Profile_Type::Gaussian, 0.);
    }

    auto s = ATOOLS::Settings::GetMainSettings()["BEAM_TRANSVERSE_PROFILE"];
    const auto model = s["TYPE"].SetDefault(Transverse_Profile_Type::Gaussian)
                           .GetTwoVector<Transverse_Profile_Type>()[index];
    if (model == Transverse_Profile_Type::WoodsSaxon) {
      // the density of the nucleus, as in its photon flux and Glauber overlap
      const Woods_Saxon density = Woods_Saxon::FromSettings(beam, index);
      const Transverse_Profile profile(density);
      msg_Info() << model << ", R = " << density.R() * hbarc
                 << " fm, d = " << density.d() * hbarc
                 << " fm, transverse RMS radius = " << profile.RMSRadius() * hbarc
                 << " fm; used only by overlap geometry.\n";
      return profile;
    }

    s.DeclareVectorSettingsWithEmptyDefault({"RMS_RADIUS"});
    // Allow override of radius
    const bool has_rms = s["RMS_RADIUS"].IsSetExplicitly() ||
                        !s["RMS_RADIUS"].GetVector<double>().empty();
    const double rms_fm = has_rms
        ? s["RMS_RADIUS"].GetTwoVector<double>()[index]
        // the radius of an ion is that of a uniform sphere, <s^2> = 2 R^2 / 5,
        // that of other beams an RMS radius, <s^2> = 2 R^2 / 3
        : std::sqrt(beam.IsIon() ? 2. / 5. : 2. / 3.) * beam.Radius();
    msg_Info() << model << ", transverse RMS radius = " << rms_fm << " fm"
               << (rms_fm == 0. ? " (pointlike)" : "")
               << "; used only by overlap geometry.\n";
    return Transverse_Profile(model, rms_fm / hbarc);
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
}
