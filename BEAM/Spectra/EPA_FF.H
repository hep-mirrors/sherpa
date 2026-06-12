#ifndef BEAM_Spectra_EPA_FF_H
#define BEAM_Spectra_EPA_FF_H

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Lookup_Tables.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Special_Functions.H"
#include "ATOOLS/Phys/Flavour.H"
#include <ATOOLS/Org/Run_Parameter.H>
#include <cmath>
#include <functional>

namespace ATOOLS {
  class Settings;
}

namespace BEAM {
  enum class EPA_ff_type {
    point = 0,
    pointApprox = 1,
    proton = 2,
    protonApprox = 3,
    ionApproxInt = 4,
    Gauss = 10,
    hcs = 11,
    dipole = 12,
    dipoleApprox = 13,
    WoodSaxon = 14,
    WoodSaxonApprox = 15,
    ionApprox = 16,
    Test = 99,
    Undefined = -1
  };
  inline std::ostream& operator<<(std::ostream& os, EPA_ff_type type)
  {
    switch (type) {
    case EPA_ff_type::point:
      return os << "Lepton";
    case EPA_ff_type::pointApprox:
      return os << "Approx_Lepton";
    case EPA_ff_type::proton:
      return os << "Proton";
    case EPA_ff_type::protonApprox:
      return os << "Approx_Proton";
    case EPA_ff_type::ionApproxInt:
      return os << "Point-like_Integrated";
    case EPA_ff_type::Gauss:
      return os << "Gaussian";
    case EPA_ff_type::hcs:
      return os << "HCS";
    case EPA_ff_type::dipole:
      return os << "Dipole";
    case EPA_ff_type::dipoleApprox:
      return os << "Approx_Dipole";
    case EPA_ff_type::WoodSaxon:
      return os << "Woods-Saxon";
    case EPA_ff_type::WoodSaxonApprox:
      return os << "Approx_Woods-Saxon";
    case EPA_ff_type::ionApprox:
      return os << "Point-like";
    case EPA_ff_type::Test:
      return os << "Test";
    case EPA_ff_type::Undefined:
      return os << "Undefined";
    }

    return os << "Unknown EPA form factor " << type << "\n";
  }

  inline std::istream& operator>>(std::istream& is, EPA_ff_type& type)
  {
    std::string token;
    if (!(is >> token)) {
      return is;
    }

    if (token == "Lepton" || token == "0")
      type = EPA_ff_type::point;
    else if (token == "Approx_Lepton" || token == "1")
      type = EPA_ff_type::pointApprox;
    else if (token == "Proton" || token == "2")
      type = EPA_ff_type::proton;
    else if (token == "Approx_Proton" || token == "3")
      type = EPA_ff_type::protonApprox;
    else if (token == "Point-like_Integrated" || token == "4")
      type = EPA_ff_type::ionApproxInt;
    else if (token == "Gaussian" || token == "10")
      type = EPA_ff_type::Gauss;
    else if (token == "HCS" || token == "11")
      type = EPA_ff_type::hcs;
    else if (token == "Dipole" || token == "12")
      type = EPA_ff_type::dipole;
    else if (token == "Approx_Dipole" || token == "13")
      type = EPA_ff_type::dipoleApprox;
    else if (token == "Woods-Saxon" || token == "14")
      type = EPA_ff_type::WoodSaxon;
    else if (token == "Approx_Woods-Saxon" || token == "15")
      type = EPA_ff_type::WoodSaxonApprox;
    else if (token == "Point-like" || token == "16")
      type = EPA_ff_type::ionApprox;
    else if (token == "Test" || token == "99")
      type = EPA_ff_type::Test;
    else {
      type = EPA_ff_type::Undefined;
    }

    return is;
  }

  class EPA_FF_Base {
  protected:

    ATOOLS::Flavour m_beam;
    size_t m_A;
    double m_mass, m_mass2, m_R;
    double m_q2min, m_q2max, m_pt2max;
    double m_xmin, m_xmax;
    double m_bmin,
        m_bmax; //< bmin, bmax only relative to radius
    double m_chi_cut; //< upper cut on chi = x*m*b
    double m_Zsquared, m_b,
        m_b_pl_threshold; //< b, b_pl_threshold are absolute distances
    std::unique_ptr<ATOOLS::TwoDim_Table> p_N_xb;

    size_t m_nxbins, m_nbbins;

    // Importance-sample the impact parameter for this form factor's beam
    // parameters; see SampleB for the 1/x support cap. Stores b in m_b.
    double SampleImpactParameter(const double& x, const double& ran)
    {
      return SampleB(x, ran, m_bmin, m_bmax, m_chi_cut, m_mass, m_R, m_b);
    }
    virtual void FillTables();

    // ---- table caching (implemented in EPA_FF.C) -------------------------
    // The on-disk format and the MPI-safe load/save go through the generic
    // ATOOLS::Table_Cache facility; the logic below only builds the cache key
    // and decides where the file lives. The tables are stored *inside* the
    // result archive <RESULT_DIRECTORY>.zip (as entries EPA/<file>.tab), written
    // via the My_File zip-DB just like the matrix-element results -- this is why
    // BuildNxbTable opens that DB around the cache I/O. Caching is therefore
    // gated on GENERATE_RESULT_DIRECTORY exactly as the ME results are.

    // Content version baked into every cache key (see BaseCacheId). Bump
    // whenever the table-construction logic changes in a way no other key
    // parameter captures, to invalidate previously written caches.
    // v2: fixed the EPA_HCS form factor argument (q*R instead of Q^2).
    static constexpr int s_cache_version = 2;
    bool m_cache{true};    //< caching enabled (EPA:CacheTables) for this beam
    std::string m_respath; //< resolved RESULT_DIRECTORY root; cache lives inside
                           //< <m_respath>.zip; empty => caching disabled

    // Common part of the cache key (version + axes + beam scalars). Concrete
    // form factors prepend a type tag and append their own parameters so that
    // any change to a relevant input invalidates the cache exactly.
    std::string BaseCacheId(const std::string& tag) const;
    virtual std::string CacheId() const { return BaseCacheId("base"); }
    std::string CacheFileName(const std::string& tag,
                              const std::string& key) const;

    // True iff table caching is active for this beam: enabled via
    // EPA:CacheTables AND a result archive is being written. An empty m_respath
    // means GENERATE_RESULT_DIRECTORY is off or RESULT_DIRECTORY is empty, i.e.
    // there is no <RESULT_DIRECTORY>.zip to cache into.
    bool CachingActive() const { return m_cache && !m_respath.empty(); }

    // Load-or-(fill + save) a cached table via ATOOLS::Table_Cache. On a hit
    // every rank loads it (rank 0 reads and broadcasts the bytes); on a miss
    // every rank runs fill() (which populates target) and rank 0 writes the
    // result for the next run. Shared by the N(x,b) and form-factor tables.
    template <class Table>
    void LoadOrBuild(const std::string& tag, const std::string& key,
                     std::unique_ptr<Table>& target,
                     const std::function<void()>& fill);

    // Build p_N_xb: opens the result zip-DB (so this and any prerequisite
    // form-factor table land inside <RESULT_DIRECTORY>.zip), then PrepareFill()
    // builds that prerequisite (e.g. the Woods-Saxon form-factor table) before
    // FillTables() fills the N(x,b) grid.
    void BuildNxbTable();
    virtual void PrepareFill() {}

  public:

    EPA_FF_Base(const ATOOLS::Flavour& beam, int dir);
    virtual ~EPA_FF_Base() = default;

    // Pure importance-sampling of the (dimensionless) impact parameter for a
    // given x. The upper limit is constrained by the physical b_max and the
    // chi_max set by the user. The latter should be preferred as it does cut
    // including the correct scaling with x.
    static double SampleB(double x, double ran, double bmin, double bmax,
                          double chi_cut, double mass, double R, double& b_phys)
    {
      const double bmaxx =
          ATOOLS::Max(bmin, ATOOLS::Min(bmax, chi_cut / (x * mass * R)));
      const double norm =
          0.5 * std::log((ATOOLS::sqr(bmaxx) + 1.) / (ATOOLS::sqr(bmin) + 1.));
      // Sampling from b_abs/(b_abs^2 + R^2) with b_abs = b * R
      double b =
          std::sqrt((ATOOLS::sqr(bmin) + 1.) * std::exp(2. * ran * norm) - 1.);
      double weight = norm * (ATOOLS::sqr(b) + 1.) / b * R;
      b_phys = b * R; // convert to physical distance
      return weight;
    }

    virtual double N(const double& x, const double& ran)
    {
      double wt = SampleImpactParameter(x, ran);
      DEBUG_FUNC("x = " << x << ", ran = " << ran << ", b = " << m_b);
      if (m_b < m_b_pl_threshold) return (*p_N_xb)(x, m_b) * wt;
      // use point-like approximation for large impact parameters
      double chi = x * m_mass * m_b;
      return 2 * m_Zsquared * m_b * x * m_mass2 *
             (ATOOLS::sqr(ATOOLS::SF.Kn(1, chi)) +
              2. * m_mass / ATOOLS::rpa->gen.Ecms() *
                  ATOOLS::sqr(ATOOLS::SF.Kn(0, chi))) *
             wt;
    }
    virtual double FF(const double& Q2) = 0;

    void SetRadius(const double& R) { m_R = R; }
    void SetPT2Max(const double& pt2max) { m_pt2max = pt2max; }
    void SetQ2Max(const double q2max) { m_q2max = q2max; }

    double Mass2() const { return m_mass2; }
    double Radius() const { return m_R; }
    double ImpactParameter() const { return m_b; }
    double Xmin() const { return m_xmin; }
    double Xmax() const { return m_xmax; }

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
                               (m_pt2max + ATOOLS::sqr(m_mass * x)) / (1 - x))
                 : m_q2max;
    }

    void OutputToCSV(const std::string& type);
  };

  class N_xb_int: public ATOOLS::Function_Base {
  private:

    EPA_FF_Base* p_ff;
    double m_b, m_x;

  public:

    explicit N_xb_int(EPA_FF_Base* ff): p_ff(ff), m_b(0.), m_x(0.) {}
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
      double Q2 = (qT2 + m_x * m_x * p_ff->Mass2()) / (1. - m_x);
      // integration variable is qT, so dy/dqT = b, hence Jacobian 1/b needed
      double res = qT2 / Q2 * p_ff->FF(Q2) / m_b;
      return res;
    }
  };

  class EPA_Point final: public EPA_FF_Base {
  public:

    EPA_Point(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double N(const double& x, const double& ran) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_PointApprox final: public EPA_FF_Base {
  public:

    EPA_PointApprox(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double N(const double& x, const double& ran) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Proton final: public EPA_FF_Base {
  private:

    double m_Q02, m_mu2;

  public:

    EPA_Proton(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double FF(const double& Q2) override
    {
      return std::pow(1. + Q2 / m_Q02, -2);
    }
    double N(const double& x, const double& ran) override;
  };

  class EPA_ProtonApprox final: public EPA_FF_Base {
  private:

    double m_mu2;

  public:

    EPA_ProtonApprox(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double FF(const double& Q2) override { return 1.; }
    double N(const double& x, const double& ran) override;
  };

  class EPA_Gauss final: public EPA_FF_Base {
  private:

    double m_Q02;
    std::string CacheId() const override;

  public:

    EPA_Gauss(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override { return exp(-Q2 / m_Q02 / 2.); }
  };

  class EPA_HCS final: public EPA_FF_Base {
  private:

    std::string CacheId() const override;

  public:

    EPA_HCS(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& QR) override
    {
      return 3 *
             (std::sin(QR) / std::pow(QR, 3) - std::cos(QR) / ATOOLS::sqr(QR));
    }
  };

  class EPA_Dipole final: public EPA_FF_Base {
  private:

    double m_Q02;
    std::string CacheId() const override;

  public:

    explicit EPA_Dipole(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      return std::pow(1. + Q2 / m_Q02, -2);
    }
  };

  class EPA_DipoleApprox final: public EPA_FF_Base {
  private:

    std::string CacheId() const override;

  public:

    explicit EPA_DipoleApprox(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_WoodSaxon final: public EPA_FF_Base {
  private:

    double m_d, m_R_WS;
    size_t m_rnodes;     //< number of r-grid intervals for the Filon transform
    double m_rmaxfactor; //< r_max = R_WS + m_rmaxfactor * d
    std::unique_ptr<ATOOLS::OneDim_Table> p_FF_Q;
    static constexpr double m_q_min = 0., m_q_max = 1.e2;
    static constexpr int m_q_n = 1e5;

    void InitFFTable();   //< Filon sine transform of the Woods-Saxon density
    void BuildFFQTable(); //< load-or-(InitFFTable + save) for p_FF_Q
    // RWS / d / q-axis / r-grid part of the cache key (shared by the FF-table
    // key and the N(x,b) key, since N(x,b) depends on the form factor).
    std::string WSParamString() const;
    std::string CacheId() const override;
    void PrepareFill() override { BuildFFQTable(); }

  public:

    EPA_WoodSaxon(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      const double Q = sqrt(Q2);
      if (Q > m_q_max) return 0.;
      // p_FF_Q is prebuilt on an N(x,b) cache miss (PrepareFill ->
      // BuildFFQTable). It is null only after an N(x,b) cache hit (which skips
      // the form-factor table) or when FF() is first reached via the debug CSV
      // output; build it on demand so any caller gets the correct value instead
      // of dereferencing a null table. InitFFTable() is pure local computation
      // (no MPI/cache collectives), so building it from a single rank is safe.
      if (!p_FF_Q) InitFFTable();
      return (*p_FF_Q)(Q);
    }
  };

  class EPA_WoodSaxonApprox final: public EPA_FF_Base {
  private:

    double m_R_WS, m_a;
    std::string CacheId() const override;

  public:

    explicit EPA_WoodSaxonApprox(const ATOOLS::Flavour& beam, int dir);
    double FF(const double& Q2) override
    {
      // from eq. 14 in hep-ph/9902259
      // according to Ilkka, a = r = 0.7 fm
      // Normalising the form factor/nuclear density can be performed
      // analytically and it yields \rho_0 = A / 4 / \pi * 3 / R^3
      const double qR = sqrt(Q2) * m_R_WS;
      return 3. / std::pow(qR, 3) * (std::sin(qR) - qR * std::cos(qR)) /
             (1. + ATOOLS::sqr(m_a) * Q2);
    }
  };

  class EPA_IonApprox final: public EPA_FF_Base {
  public:

    EPA_IonApprox(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override;
    double FF(const double& Q2) override { return 0.; }
  };

  class EPA_IonApproxIntegrated final: public EPA_FF_Base {
  public:

    EPA_IonApproxIntegrated(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double N(const double& x, const double& ran) override;
    double FF(const double& Q2) override { return 1.; }
  };

  class EPA_Test final: public EPA_FF_Base {
  public:

    explicit EPA_Test(const ATOOLS::Flavour& beam, int dir);
    void FillTables() override {};
    double FF(const double& Q2) override { return 0.; }
  };
} // end of namespace BEAM

#endif
