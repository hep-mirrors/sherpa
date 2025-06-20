#include "BEAM/Main/Beam_Base.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"

using namespace BEAM;

Beam_Base::Beam_Base(beamspectrum _type, const ATOOLS::Flavour& _beam,
                     const double _energy, const double _polarisation,
                     const int _dir, int mode) :
  m_type(_type), m_Nbunches(1), m_beam(_beam),
  m_position(ATOOLS::Vec4D(0.,0.,0.,0.)),
  m_dir(_dir), m_energy(_energy), m_polarisation(_polarisation),
  m_x(1.), m_Q2(0.), m_weight(1.), m_on(false) {
  double disc = mode ? 1.0 : 1.0 - ATOOLS::sqr(m_beam.Mass() / m_energy);
  if (disc < 0.) {
    msg_Error() << "Error in Beam_Base :" << m_type << std::endl
                << "   Mismatch of energy and mass of beam particle : "
                << m_beam << " / " << m_energy << std::endl
                << "   Will lead to termination of program." << std::endl;
    ATOOLS::Abort();
  }
  double pz = m_dir * m_energy * sqrt(disc);
  m_lab     = ATOOLS::Vec4D(m_energy, 0., 0., pz);
  m_bunches.resize(m_Nbunches, m_beam);
  m_vecouts.resize(m_Nbunches, m_lab);
}
