#ifndef ATOOLS_Math_Filon_Integrator_H
#define ATOOLS_Math_Filon_Integrator_H

#include "ATOOLS/Math/Function_Base.H"
#include <cstddef>
#include <vector>

namespace ATOOLS {

  //////////////////////////////////////////////////////////////////////////////
  // Filon-Simpson quadrature for oscillatory integrals
  //   I_s(k) = \int_a^b f(r) sin(k r) dr,   I_c(k) = \int_a^b f(r) cos(k r) dr.
  //
  // The smooth factor f is sampled once on a uniform grid of (2*npanels+1)
  // nodes; the oscillatory factor sin/cos(k r) is integrated analytically over
  // each panel (Abramowitz & Stegun 25.4.47). The result is therefore accurate
  // for arbitrarily large k as long as the grid resolves f itself -- unlike a
  // plain Simpson/Riemann sum, which would need to resolve the oscillation.
  //
  // Re-use the same instance for many k: f is sampled only in the constructor.
  //////////////////////////////////////////////////////////////////////////////
  class Filon_Integrator {
  private:

    double m_a, m_b, m_h;
    size_t m_npanels;        // number of Simpson double-panels (2 intervals each)
    std::vector<double> m_f; // f sampled at r_j = a + j*h, j = 0..2*npanels

    // Filon coefficients alpha/beta/gamma as functions of theta = k*h, using a
    // small-theta Taylor expansion to avoid cancellation near theta = 0.
    static void Coefficients(double theta, double& alpha, double& beta,
                             double& gamma);

  public:

    // Sample f on a uniform grid of (2*npanels+1) nodes over [a, b].
    Filon_Integrator(Function_Base* f, double a, double b, size_t npanels);
    // Take pre-sampled values (odd count >= 3, uniformly spaced over [a, b]).
    Filon_Integrator(std::vector<double> fvals, double a, double b);

    double SineTransform(double k) const;   // \int_a^b f(r) sin(k r) dr
    double CosineTransform(double k) const; // \int_a^b f(r) cos(k r) dr

    double Step() const { return m_h; }
    size_t NNodes() const { return m_f.size(); }
    double Node(size_t j) const { return m_a + double(j) * m_h; }
    double Value(size_t j) const { return m_f[j]; }
  };

} // namespace ATOOLS

#endif
