#ifndef BEAM_Spectra_EPA_FF_H
#define BEAM_Spectra_EPA_FF_H

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Org/Message.H"
#include <fstream>
#include <iostream>

namespace ATOOLS {
  class Settings;
}

namespace BEAM {
  enum class EPA_ff_type {
    point        = 0,
    proton       = 1,
    protonApprox = 2,
    Gauss        = 10,
    hcs          = 11,
    dipole       = 12,
    dipoleApprox = 13,
    WoodSaxon    = 14,
    Test         = 99
  };

  class EPA_FF_Base {
  protected:
    ATOOLS::Flavour       m_beam;
    double                m_mass, m_mass2, m_R;
    double                m_q2min, m_q2max, m_pt2max;
    double                m_xmin, m_xmax, m_bmin, m_bmax;
    ATOOLS::OneDim_Table* p_Nred_x;
    ATOOLS::TwoDim_Table *p_N_xb, *p_Inv_xb;

    size_t m_nxbins, m_nbbins;
    bool   m_approx;

    virtual void Fill_Nxb_Table(ATOOLS::axis& xaxis, ATOOLS::axis& baxis);
    virtual void Fill_Invxb_Table();
    virtual void TestInvxb();

  public:
    EPA_FF_Base(const ATOOLS::Flavour& beam, int dir);
    virtual ~EPA_FF_Base()
    {
      if (p_Nred_x) delete p_Nred_x;
      if (p_N_xb) delete p_N_xb;
      if (p_Inv_xb) delete p_Inv_xb;
    }
    virtual void   FillTables(const size_t& nx, const size_t& nb);
    virtual double operator()(const double& x, const double& Q2) = 0;
    virtual double N(const double& x, const double& b)
    {
      double w((*p_N_xb)(x, b));
      //msg_Out() << METHOD << "("<<x<<", "<<b<<") = " << w << "\n";
      return (*p_N_xb)(x, b);
    };
    virtual double SelectB(double x);
    virtual double FF(const double& Q2) = 0;

    void         SetApprox(bool approx) { m_approx = approx; }
    virtual void SetRadius(const double& R) { m_R = R; }
    virtual void SetPT2Max(const double& pt2max) { m_pt2max = pt2max; }
    void         SetQ2Max(const double q2max) { m_q2max = q2max; }

    double Mass() const { return m_beam.Mass(); }
    double Radius() const { return m_R; }

    virtual ATOOLS::TwoDim_Table* GetN_xb() const { return p_N_xb; }
    virtual ATOOLS::OneDim_Table* GetNred_x() const { return p_Nred_x; }

    inline double Q2min(const double& x) const
    {
      return ATOOLS::Max(m_q2min, ATOOLS::sqr(m_mass * x) / (1 - x));
    }
    // Use maximal angle for the scattered electron in Q^2_max,
    // compare with hep-ph/9610406 and hep-ph/9310350 for the
    // inclusion of term ~(x m)^2
    inline double Q2max(const double& x) const
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
    double operator()(double y) override;
  };

  class EPA_Point final : public EPA_FF_Base {
  public:
    EPA_Point(const ATOOLS::Flavour& beam, int dir);
    void   FillTables(const size_t& nx, const size_t& nb) override;
    double operator()(const double& x, const double& Q2) override;
    double SelectB(double x) override;
    double N(const double& x, const double& b) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Proton final : public EPA_FF_Base {
  private:
    double m_Q02, m_mu2;

  public:
    EPA_Proton(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override
    {
      return std::pow(1. + Q2 / m_Q02, -2);
    }
  };

  class EPA_ProtonApprox final : public EPA_FF_Base {
  private:
    double m_mu2;

  public:
    EPA_ProtonApprox(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Gauss final : public EPA_FF_Base {
  private:
    double m_Q02, m_Zsquared;

  public:
    EPA_Gauss(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return exp(-Q2 / m_Q02 / 2.); }
  };

  class EPA_HCS final : public EPA_FF_Base {
  private:
    double m_Zsquared;

  public:
    EPA_HCS(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& QR) override
    {
      return 3 * (std::sin(QR) / std::pow(QR, 3) -
             std::cos(QR) / ATOOLS::sqr(QR));
    }
  };

  class EPA_Dipole final : public EPA_FF_Base {
  private:
    double m_Q02, m_Zsquared;

  public:
    explicit EPA_Dipole(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return std::pow(1. + Q2 / m_Q02, -2); }
  };

  class EPA_DipoleApprox final : public EPA_FF_Base {
  private:
    double m_Q02, m_Zsquared;

  public:
    explicit EPA_DipoleApprox(const ATOOLS::Flavour& beam, int dir);
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_WoodSaxon final : public EPA_FF_Base {
  private:
    class Rho_argument : public ATOOLS::Function_Base {
    private:
      double m_R, m_d;

    public:
      Rho_argument(const double& R, const double& d) : m_R(R), m_d(d) {}
      double operator()(double r) override
      {
        return 4. * M_PI * r * r / (1. + exp((r - m_R) / m_d));
      }
    };
    class WS_potential : public ATOOLS::Function_Base {
    private:
      double m_R, m_d, m_q;

    public:
      WS_potential(const double& R, const double& d) : m_R(R), m_d(d), m_q(1.)
      {}
      inline void SetQ(const double& q) { m_q = q; }
      double      operator()(double r) override
      {
        if (m_q * r < 1.e-6)
          return 4. * M_PI * r * r / (1. + exp((r - m_R) / m_d));
        return 4. * M_PI * sin(m_q * r) / m_q * r / (1. + exp((r - m_R) / m_d));
      }
    };
    class N_argument : public ATOOLS::Function_Base {
    private:
      EPA_WoodSaxon* p_ws;
      double         m_x;

    public:
      explicit N_argument(EPA_WoodSaxon* ws) : p_ws(ws), m_x(0.) {}
      inline void SetX(const double& x) { m_x = x; }
      double operator()(double q2) override { return (*p_ws)(m_x, q2) / q2; }
    };

    double                m_d, m_rho0;
    ATOOLS::OneDim_Table *p_FF_Q2, *p_N;

    double CalculateDensity();
    void   InitFFTable(const double& q2min, const double& q2max);
    void   InitNTable(const double& xmin, const double& xmax);

  public:
    EPA_WoodSaxon(const ATOOLS::Flavour& beam, int dir);
    ~EPA_WoodSaxon() override
    {
      if (p_N) {
        delete p_N;
        p_N = nullptr;
      }
      if (p_FF_Q2) {
        delete p_FF_Q2;
        p_FF_Q2 = nullptr;
      }
    }
    void   FillTables(const size_t& nx, const size_t& nb) override;
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return (*p_FF_Q2)(Q2); }
  };

  class EPA_Test final : public EPA_FF_Base {
  public:
    explicit EPA_Test(const ATOOLS::Flavour& beam, int dir);
    void   FillTables(const size_t& m_nxbins, const size_t& m_nbbins) override;
    double operator()(const double& x, const double& Q2) override;
    double FF(const double& Q2) override { return 0.; }
  };
}// end of namespace BEAM

#endif
