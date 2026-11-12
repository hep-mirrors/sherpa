#ifndef BEAM_Main_Woods_Saxon_H
#define BEAM_Main_Woods_Saxon_H

#include <cstddef>

namespace ATOOLS {
  class Flavour;
}

namespace BEAM {

  /**
   * Woods-Saxon density of a nucleus with half-density radius R and skin depth
   * d (both in 1/GeV),
   *   rho(r) = rho0 / (1 + exp((r - R) / d)),
   * normalised to unity, \int d^3r rho(r) = 1. Multiply by the mass number A
   * for the nucleon density.
   */
  class Woods_Saxon {
  private:

    double m_R, m_d, m_rho0;

  public:

    Woods_Saxon(double R, double d);

    /// The density of the nucleus in beam 0 or 1, from EPA:WoodsSaxon_R and
    /// EPA:WoodsSaxon_d in fm, which default to 1.118 A^(1/3) and 0.54.
    static Woods_Saxon FromSettings(const ATOOLS::Flavour& nucleus,
                                    size_t beam);

    double R() const { return m_R; }
    double d() const { return m_d; }

    double Density(double r) const;
    /// Thickness T(s) = \int dz rho(sqrt(s^2 + z^2)), with \int d^2s T(s) = 1.
    double Thickness(double s) const;
    /// sqrt(<s^2>) in the transverse plane, with <s^2> = 2/3 <r^2>.
    double TransverseRMSRadius() const;
  };

} // namespace BEAM

#endif
