#ifndef BEAM_Spectra_EPA_FF_H
#define BEAM_Spectra_EPA_FF_H

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Special_Functions.H"
#include "ATOOLS/Phys/Flavour.H"
#include <ATOOLS/Org/Run_Parameter.H>
#include <cmath>

namespace ATOOLS {
  class Settings;
}

namespace BEAM {
  enum class EPA_ff_type {
    point        = 0,
    pointApprox  = 1,
    proton       = 2,
    protonApprox = 3,
    ionApproxInt = 4,
    Gauss        = 10,
    hcs          = 11,
    dipole       = 12,
    dipoleApprox = 13,
    WoodSaxon    = 14,
    WoodSaxonApprox = 15,
    ionApprox    = 16,
    Test         = 99
  };

  class EPA_FF_Base {
  private:
    double m_norm_distribution;

  protected:
    ATOOLS::Flavour m_beam;
    size_t          m_A;
    double          m_mass, m_mass2, m_R;
    double          m_q2min, m_q2max, m_pt2max;
    double          m_xmin, m_xmax;
    double          m_bmin,
            m_bmax;//< bmin, bmax only relative to radius
    double m_Zsquared, m_b,
            m_b_pl_threshold;//< b, b_pl_threshold are absolute distances
    std::unique_ptr<ATOOLS::TwoDim_Table> p_N_xb;

    size_t m_nxbins, m_nbbins;

    double SampleImpactParameter()
    {
      // Sampling from b_abs/(b_abs^2 + R^2) with b_abs = b * R
      m_b = sqrt((ATOOLS::sqr(m_bmin) + 1.) * std::exp(2. * ATOOLS::ran->Get() *
                                                       m_norm_distribution) -
                 1.);
      double weight = m_norm_distribution * (ATOOLS::sqr(m_b) + 1.) / m_b * m_R;
      // multiply radius to convert b to physical distance
      m_b *= m_R;
      return weight;
    };
    virtual void FillTables();

  public:
    EPA_FF_Base(const ATOOLS::Flavour& beam, int dir);
    virtual ~EPA_FF_Base() = default;
    virtual double N(const double& x)
    {
      double wt = SampleImpactParameter();
      if (m_b < m_b_pl_threshold) return (*p_N_xb)(x, m_b) * wt;
      // use point-like approximation for large impact parameters
      double chi = x * m_mass * m_b;
      return 2 * m_Zsquared * m_b * x * m_mass2 *
             (ATOOLS::sqr(ATOOLS::SF.Kn(1, chi)) +
              2. * m_mass / ATOOLS::rpa->gen.Ecms() *
                      ATOOLS::sqr(ATOOLS::SF.Kn(0, chi))) *
             wt;
    };
    virtual double FF(const double& Q2) = 0;

    void SetRadius(const double& R) { m_R = R; }
    void SetPT2Max(const double& pt2max) { m_pt2max = pt2max; }
    void SetQ2Max(const double q2max) { m_q2max = q2max; }

    double Mass() const { return m_mass; }
    double Radius() const { return m_R; }
    double ImpactParameter() const { return m_b; }

    double Q2min(const double& x) const
    {
      return ATOOLS::Max(m_q2min, ATOOLS::sqr(m_mass * x) / (1 - x));
    }
    // Use maximal angle for the scattered electron in Q^2_max,
    // compare with hep-ph/9610406 and hep-ph/9310350 for the
    // inclusion of term ~(x m)^2
    double Q2max(const double& x) const
    {
      return (m_pt2max > 0.)
                     ? ATOOLS::Min(m_q2max,
                                   (m_pt2max + ATOOLS::sqr(m_mass * x)) /
                                           (1 - x))
                     : m_q2max;
    }

    void OutputToCSV(const std::string& type);
  };

  class N_xb_int : public ATOOLS::Function_Base {
  private:
    EPA_FF_Base* p_ff;
    double       m_b, m_x;

  public:
    explicit N_xb_int(EPA_FF_Base* ff) : p_ff(ff), m_b(0.), m_x(0.) {}
    inline void SetXB(const double& x, const double& b)
    {
      m_x = x;
      m_b = b;
    }
    double operator()(double y) override
    {
      //////////////////////////////////////////////////////////////////////////
      // Integration argument y here is bT*qT as mandated by the Bessel
      // function:
      // - argument of form factor Q^2 = qT^2+x^2m^2 with qT^2 = y^2/bT^2
      // - we assume that m_b is in units of 1/GeV, qT is in GeV, overall
      // results are in GeV.
      //////////////////////////////////////////////////////////////////////////
      double qT = y / m_b, qT2 = ATOOLS::sqr(qT);
      double Q2 = (qT2 + ATOOLS::sqr(m_x * p_ff->Mass())) / (1. - m_x);
      // integration variable is qT, so dy/dqT = b, hence Jacobian 1/b needed
      double res = qT2 / Q2 * p_ff->FF(Q2) / m_b;
      return res;
    }
  };

  class EPA_Point final : public EPA_FF_Base {
  public:
    EPA_Point(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double N(const double& x) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_PointApprox final : public EPA_FF_Base {
  public:
    EPA_PointApprox(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double N(const double& x) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Proton final : public EPA_FF_Base {
  private:
    double m_Q02, m_mu2;

  public:
    EPA_Proton(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double FF(const double& Q2) override
    {
      return std::pow(1. + Q2 / m_Q02, -2);
    }
    double N(const double& x) override;
  };

  class EPA_ProtonApprox final : public EPA_FF_Base {
  private:
    double m_mu2;

  public:
    EPA_ProtonApprox(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double FF(const double& Q2) override { return 1.; }
    double N(const double& x) override;
  };

  class EPA_Gauss final : public EPA_FF_Base {
  private:
    double m_Q02;

  public:
    EPA_Gauss(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override { return exp(-Q2 / m_Q02 / 2.); }
  };

  class EPA_HCS final : public EPA_FF_Base {
  public:
    EPA_HCS(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& QR) override
    {
      return 3 *
             (std::sin(QR) / std::pow(QR, 3) - std::cos(QR) / ATOOLS::sqr(QR));
    }
  };

  class EPA_Dipole final : public EPA_FF_Base {
  private:
    double m_Q02;

  public:
    explicit EPA_Dipole(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      return std::pow(1. + Q2 / m_Q02, -2);
    }
  };

  class EPA_DipoleApprox final : public EPA_FF_Base {
  private:
    double m_Q02;

  public:
    explicit EPA_DipoleApprox(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_WoodSaxon final : public EPA_FF_Base {
  private:
    double                                m_d, m_R_WS, m_rho0;
    std::unique_ptr<ATOOLS::OneDim_Table> p_FF_Q2;
    static constexpr double               m_q2_min = 1.e-12, m_q2_max = 1.e5;
    static constexpr int                  m_q2_n = 1e7;

    static double IntegrateWithAdaptiveRange(const std::function<double(double)>& integrand,
        double initial_rmax, double tolerance);
    double CalculateDensity();
    void   InitFFTable();

  public:
    EPA_WoodSaxon(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      if (Q2 < m_q2_min) return 1.;
      if (Q2 > m_q2_max) return 0.;
      return (*p_FF_Q2)(Q2);
    }
  };

  class EPA_WoodSaxonApprox final : public EPA_FF_Base {
  private:
    double m_R_WS, m_a;

  public:
    explicit EPA_WoodSaxonApprox(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      // from eq. 14 in hep-ph/9902259
      // according to Ilkka, a = r = 0.7 fm
      // Normalising the form factor/nuclear density can be performed
      // analytically and it yields \rho_0 = A / 4 / \pi * 3 / R^3
      const double qR = sqrt(Q2) * m_R_WS;
      return 3. / std::pow(qR, 3) * (std::sin(qR) - qR * std::cos(qR)) / (1. + ATOOLS::sqr(m_a) * Q2);
    }
  };

  class EPA_IonApprox final : public EPA_FF_Base {
  public:
    EPA_IonApprox(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override;
    double FF(const double& Q2) override { return 0.; }
  };

  class EPA_IonApproxIntegrated final : public EPA_FF_Base {
  public:
    EPA_IonApproxIntegrated(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double N(const double& x) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Test final : public EPA_FF_Base {
  public:
    explicit EPA_Test(const ATOOLS::Flavour& beam, int dir);
    void   FillTables() override {};
    double FF(const double& Q2) override { return 0.; }
  };
}// end of namespace BEAM

#endif
