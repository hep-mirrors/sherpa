#include "BEAM/Spectra/Target_Matter.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace BEAM;

Flavour Target_Matter::SampleFlavour(const double& Z, const double& A)
{
  return (ran->Get() < Z / A) ? Flavour(kf_p_plus) : Flavour(kf_n);
}

Target_Matter::Target_Matter(const double& Z, const double& A,
                             const double& energy, const double& polarisation,
                             const int& dir)
    // Mass/momentum stay fixed for the run, so this reports as monochromatic.
    : Beam_Base(beamspectrum::monochromatic, SampleFlavour(Z, A), energy,
               polarisation, dir),
      m_Z(Z), m_A(A) {}


Beam_Base * Target_Matter::Copy()
{
  return new Target_Matter(m_Z, m_A, m_energy, m_polarisation, m_dir);
}
