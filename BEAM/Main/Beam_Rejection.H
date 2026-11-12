#ifndef BEAM_Main_Beam_Rejection_H
#define BEAM_Main_Beam_Rejection_H

#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "BEAM/Main/Woods_Saxon.H"
#include <cmath>

namespace BEAM {
  class Beam_Rejection_Base {
  protected:

    double m_radius1, m_radius2;

  public:

    Beam_Rejection_Base(const ATOOLS::Flavour& flav1,
                        const ATOOLS::Flavour& flav2)
        : m_radius1(flav1.Radius() / ATOOLS::rpa->hBar_c()),
          m_radius2(flav2.Radius() / ATOOLS::rpa->hBar_c())
    {
    }
    virtual ~Beam_Rejection_Base() = default;

    virtual double operator()(double b) = 0;
  };

  class Radius_Rejection: public Beam_Rejection_Base {
  public:

    Radius_Rejection(const ATOOLS::Flavour& flav1, const ATOOLS::Flavour& flav2)
        : Beam_Rejection_Base(flav1, flav2)
    {
    }
    double operator()(double b) final
    {
      return b < (m_radius1 + m_radius2) ? 0. : 1.;
    }
  };

  class Proton_Proton_Rejection: public Beam_Rejection_Base {
  private:

    double m_sNN;

  public:

    Proton_Proton_Rejection(const ATOOLS::Flavour& flav1,
                            const ATOOLS::Flavour& flav2, double sNN)
        : Beam_Rejection_Base(flav1, flav2), m_sNN(sNN)
    {
    }
    double operator()(double b) final
    {
      // Taken from fit in Section 3 in arxiv:2207.03012; given in GeV^-2
      constexpr double A(9.81), B(0.211), C(0.0185);
      double b0 = A + B * std::log(m_sNN) + C * ATOOLS::sqr(std::log(m_sNN));
      return ATOOLS::sqr(1. - std::exp(-1. * ATOOLS::sqr(b) / 2 / b0));
    }
  };

  class Proton_Nucleon_Rejection: public Beam_Rejection_Base {
  public:

    Proton_Nucleon_Rejection(const ATOOLS::Flavour& flav1,
                             const ATOOLS::Flavour& flav2)
        : Beam_Rejection_Base(flav1, flav2)
    {
    }
    double operator()(double b) final
    {
      THROW(not_implemented, "Proton-Nucleon rejection not yet implemented");
    }
  };

  // Probability that two identical nuclei do not interact hadronically at
  // impact parameter b, in the optical Glauber model (arxiv:1710.07098):
  //   P(b) = exp[-sigma_NN T_AA(b)],
  //   T_AA(b) = \int d^2s T_A(s) T_A(|b - s|),
  //   T_A(s) = \int dz rho(sqrt(s^2 + z^2)),
  // with rho the Woods-Saxon density normalised to A, and sigma_NN following
  // from sNN (in GeV^2).
  class Nucleus_Nucleus_Rejection: public Beam_Rejection_Base {
  private:

    double m_sigmaNN, m_bmax;
    // log(T_AA) is nearly linear in b where P(b) rises from zero to one
    ATOOLS::OneDim_Table m_lnTAA;

  public:

    Nucleus_Nucleus_Rejection(const ATOOLS::Flavour& flav1,
                              const ATOOLS::Flavour& flav2,
                              const Woods_Saxon& density, double sNN);
    double TAA(double b) const
    {
      return b < m_bmax ? std::exp(m_lnTAA(b)) : 0.;
    }
    double operator()(double b) final { return std::exp(-m_sigmaNN * TAA(b)); }
  };
} // namespace BEAM
#endif
