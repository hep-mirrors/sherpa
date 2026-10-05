#include "BEAM/Spectra/Gaussian.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

#include <cmath>

using namespace ATOOLS;
using namespace BEAM;

Gaussian::Gaussian(const Flavour _beam, const double _energy,
                   const double _spread, const double _nsigma,
                   const double _polarisation, const int _dir)
    : Beam_Base(beamspectrum::Gaussian, _beam,
                _energy * (1. + _nsigma * _spread), _polarisation, _dir),
      m_spread(_spread), m_nsigma(_nsigma)
{
  if (m_spread <= 0.)
    THROW(fatal_error, "Gaussian beam spectrum needs a positive BEAM_SPREAD.");
  if (m_nsigma <= 0.)
    THROW(fatal_error, "Gaussian beam spectrum needs a positive BEAM_SPREAD_NSIGMA.");
  if (m_nsigma * m_spread >= 1.)
    THROW(fatal_error, "Gaussian beam spectrum: BEAM_SPREAD * BEAM_SPREAD_NSIGMA "
                       "must be < 1.");

  m_x0     = 1. / (1. + m_nsigma * m_spread);
  m_sigmax = m_spread * m_x0;
  m_xmin   = m_x0 - m_nsigma * m_sigmax;
  m_norm   = 1. / (sqrt(2. * M_PI) * m_sigmax * erf(m_nsigma / sqrt(2.)));
  m_on     = true;

  msg_Info() << "Gaussian beam spectrum for " << m_beam << ":\n"
             << "   nominal energy " << NominalEnergy() << " GeV, spread "
             << m_spread * 100. << "% (" << m_spread * NominalEnergy() * 1000.
             << " MeV), truncated at " << m_nsigma << " sigma.\n"
             << "   beam momentum built at the maximum energy " << m_energy
             << " GeV, x in [" << m_xmin << ", 1].\n";
}

Beam_Base* Gaussian::Copy()
{
  return new Gaussian(m_beam, NominalEnergy(), m_spread, m_nsigma,
                      m_polarisation, m_dir);
}

bool Gaussian::CalculateWeight(const double _x, const double _scale)
{
  m_x  = _x;
  m_Q2 = _scale;
  if (_x < m_xmin || _x > 1.) {
    m_weight = 0.;
    return false;
  }
  const double d = (_x - m_x0) / m_sigmax;
  m_weight = m_norm * exp(-0.5 * d * d);
  return true;
}
