#ifndef METOOLS_Current_Library_FF_P_X_H
#define METOOLS_Current_Library_FF_P_X_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"

namespace METOOLS {

  //////////////////////////////////////////////////////////////////////////
  //
  // Generic q^2-dependent INVARIANT form factor F(q^2) for a
  // heavy-to-light weak transition P(p) -> X(p_X), q = p - p_X.
  //
  // This is the workhorse behind ALL of VA_P_P / VA_P_S / VA_P_V /
  // VA_P_A / VA_P_T: each of those currents owns one FF_P_X instance
  // per invariant form factor (f+,f0 / f+^S,f0^S / V,A0,A1,A2 /
  // A,V0,V1,V2 / V_T,A0^T,A1^T,A2^T) and never implements any q^2
  // shape itself. Every shape listed in Sec."Generic form-factor
  // parametrizations" and Sec."z-expansion parametrizations" of
  // charm_semileptonic.tex is available here and is selectable, per
  // form factor, from the decay YAML.
  //
  // DELIBERATE DESIGN CHOICE, differing from FF_0_PP/FF_0_PPP: those
  // classes are CHANNEL-indexed (one class per hadronic final state,
  // model number selects a whole resonance tower). Here the object is
  // FORM-FACTOR-indexed instead - the "channel" enters only through
  // M, m and the numerical parameters, because for semileptonic decays
  // the Lorentz decomposition is fixed by J^P of X alone and every
  // channel reuses literally the same scalar functions. Mixing the two
  // schemes would have meant ~25 near-identical channel classes.
  //
  //////////////////////////////////////////////////////////////////////////

  enum class ffq2_shape {
    // --- 0-9: elementary shapes, charm_semileptonic.tex
    //     Sec."Generic form-factor parametrizations" ---
    constant      = 0,   // F(0)
    taylor        = 1,   // F(0)[1+c1 x+c2 x^2],       x=q^2/Lambda^2
    simple_pole   = 2,   // F(0)/(1-q^2/mpole^2)
    modified_pole = 3,   // F(0)/[(1-x)(1-alpha x)],   x=q^2/mpole^2
    BK_plus       = 4,   // Becirevic-Kaidalov f+ branch
    BK_zero       = 5,   // Becirevic-Kaidalov f0 branch
    BZ            = 6,   // Ball-Zwicky two-component
    rational      = 7,   // F(0)/(1-a x+b x^2)
    exponential   = 8,   // F(0)exp(a x+b x^2)
    multipole     = 9,   // sum_n R_n/(m_n^2-q^2)
    ISGW2         = 10,  // f+(qmax^2)[1+r^2/12 (qmax^2-q^2)]^-2
    // --- 11-15: z-expansions, Sec."z-expansion parametrizations" ---
    BGL           = 11,  // 1/(P phi) sum a_n z^n, charm outer function
    charm_series  = 12,  // F(0)-normalised charm z-series (r_n=a_n/a_0)
    BCL_plus      = 13,  // BCL with the n/N z^N threshold subtraction
    BCL_zero      = 14,  // BCL scalar branch, optional 1-q^2/m0+^2
    SSE           = 15,  // pole-factored shifted z-expansion (Eq.SSE)
    unknown       = 999
  };

  //////////////////////////////////////////////////////////////////////////
  // Shared small helpers for the D->X currents. Kept here (rather than
  // in Tools.H) because they are only ever needed by the P->X family.
  //////////////////////////////////////////////////////////////////////////

  /// Minkowski dot product of a complex and a real four-vector,
  /// written out component-wise on purpose: it must not silently pick
  /// up a conjugation from an operator* overload.
  inline Complex CDot(const ATOOLS::Vec4C & a,const ATOOLS::Vec4D & b) {
    return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
  }

  /// NOTE on the Levi-Civita contraction: eps^{mu nu rho sigma}
  /// a_nu b_rho c_sigma with a COMPLEX first argument is provided
  /// directly by ATOOLS::cross(const Vec4C&, const Vec4D&, const
  /// Vec4D&) - the legacy HADRONS++ currents use exactly that call,
  /// cross(conj(eps),p0,p1). An earlier version of this header carried
  /// a hand-rolled real/imaginary split (CCross) written because that
  /// overload had not been confirmed to exist. It does, so the helper
  /// is gone: one call site, one convention, nothing to keep in sync.
  ///
  /// The SIGN convention is likewise settled by the legacy code rather
  /// than assumed. HADRONS++ multiplies the Levi-Civita term of its
  /// charm currents by fac_costhl = -1 (and by +1 for bottom), where
  /// the axial terms carry an explicit i that the vector term does
  /// not. Factoring that global i out of the legacy current maps its
  /// charm convention onto EPS_SIGN = +1 here, which is the default.
  /// The bottom convention would be EPS_SIGN = -1. Since this decay
  /// table is charm only, the default is the validated one.

  class FF_P_X : public FormFactor_Base {
  protected:
    ffq2_shape m_shape;

    // Normalisation. m_F0 is F(q^2=0) for every shape that is
    // normalised that way; for ISGW2 it is F(q_max^2) and for
    // multipole/BGL it is an overall multiplicative factor (1 by
    // default) on top of the residues/coefficients.
    double m_F0;

    // Elementary-shape parameters.
    double m_mpole, m_mpole2;      // simple/modified pole, BZ, SSE, BCL
    double m_alpha, m_beta, m_r;   // modified pole / BK / BZ
    double m_c1, m_c2;             // Taylor
    double m_a, m_b, m_Lambda2;    // rational, exponential
    double m_rISGW2;               // ISGW2 slope [GeV^-1]

    // Multipole / residue expansion.
    std::vector<double> m_Rn, m_mn;
    bool                m_superconv; // impose sum_n R_n = 0 on the last

    // z-expansion block.
    std::vector<double> m_an;      // BGL a_n / charm r_n / BCL b_n / SSE a_n
    std::vector<double> m_blaschke;// subthreshold pole masses for P_F
    double m_mc;                   // charm mass in the outer function
    double m_m0plus;               // optional explicit 0+ pole for BCL_zero
    bool   m_useBlaschke;
    /// Cached <name>_T0 override. Cached, not re-read, because
    /// p_model dangles after the owning current's setup returns.
    double m_t0user;
    bool   m_hast0user;

    // Kinematics. M is always the nominal parent mass. m is the
    // NOMINAL daughter mass by default; with m_dynamicmass=true the
    // event-by-event sqrt(s_X) is used instead everywhere t_+/t_- and
    // q_max^2 appear (note item 3/4: consistent, but an additional
    // dynamical model - hence off by default).
    double m_M, m_m;
    bool   m_dynamicmass;

    // q^2 floor for the regularised 1/q^2 structures, see
    // VA_P_X_Base::RegularisedOverQ2.
    double m_q2min;

    void   ReadParameters();
    void   SetupZ(const double & m,double & tp,double & tm,double & t0) const;

    // z(q^2,t0) and the regularised ratio z(q^2,a)/(a-q^2), which has
    // the finite analytic limit 1/(4(t_+-a)) at q^2=a. Needed by the
    // charm outer function, where several factors are individually
    // 0/0 (note, text under Eq.(charmphi)).
    double Zvar(const double & q2,const double & t0,const double & tp) const;
    double ZOverDelta(const double & q2,const double & a,
                      const double & tp) const;
    double OuterPhiPlus(const double & q2,const double & t0,
                        const double & tp,const double & tm) const;
    double BlaschkeP(const double & q2,const double & tp) const;
    void   PrintSetup() const;

  public:
    FF_P_X(const FF_Parameters & params);
    ~FF_P_X() {}

    /// Real-valued core: F(q^2) for a given daughter mass mX.
    double Value(const double & q2,const double & mX) const;

    /// FormFactor_Base interface. moms are the FULL external momenta;
    /// m_pi[0] is the parent, m_pi[1] the daughter hadron, so that
    /// q = p_parent - p_daughter (NOT the sum, unlike FF_0_PP).
    Complex operator()(const ATOOLS::Vec4D_Vector& moms);

    /// Convenience for the currents, which already have q2/mX to hand.
    inline double operator()(const double & q2,const double & mX) const {
      return Value(q2,mX);
    }

    inline double        GetF0()  const          { return m_F0; }
    inline void          SetF0(const double & f) { m_F0 = f;    }
    inline ffq2_shape    Shape()  const          { return m_shape; }
    inline bool          DynamicMass() const     { return m_dynamicmass; }
    inline double        NominalMass() const     { return m_m; }
  };
}

#endif
