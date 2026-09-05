#ifndef BEAM_Spectra_Gaussian_H
#define BEAM_Spectra_Gaussian_H

#include "BEAM/Main/Beam_Base.H"

namespace BEAM {
  /*!
    \class Gaussian
    \brief Gaussian beam energy spread, truncated at a finite number of sigma.

  */
  class Gaussian : public Beam_Base {
  private:
    /// Relative spread sigma_E/E of the nominal beam energy.
    double m_spread;
    /// Truncation of the Gaussian, in units of sigma.
    double m_nsigma;
    /// Centre of the x distribution, = E_nominal/E_max.
    double m_x0;
    /// Width of the x distribution, = m_spread/(1+m_nsigma*m_spread).
    double m_sigmax;
    /// Lower edge of the sampled range; the upper edge is exactly 1.
    double m_xmin;
    /// 1/(sqrt(2 pi) sigma_x erf(n/sqrt(2))), the truncated normalisation.
    double m_norm;

  public:
    Gaussian(const ATOOLS::Flavour beam, const double energy,
             const double spread, const double nsigma,
             const double polarisation, const int dir);

    Beam_Base* Copy() override;
    bool CalculateWeight(const double x, const double scale) override;

    double Xmin() const override { return m_xmin; }
    double Xmax() const override { return 1.; }
    double Peak() const override { return m_x0; }
    /*!
    */
    double Exponent() const override { return 0.; }

    /// Nominal (mean) beam energy, as opposed to Energy(), which is the maximum.
    double NominalEnergy() const { return m_energy * m_x0; }
    const double& Spread() const { return m_spread; }
  };
}
#endif
