#ifndef METOOLS_Current_Library_FF_PP_Wave_H
#define METOOLS_Current_Library_FF_PP_Wave_H

#include "METOOLS/HadronCurrents/FormFactors/FormFactor_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"

namespace METOOLS {

  //////////////////////////////////////////////////////////////////////////
  //
  // Complex partial-wave amplitude A_L(s) for a two-pseudoscalar
  // system produced in a semileptonic decay, s = (p_1+p_2)^2.
  //
  // BUILT ON THE EXISTING PROPAGATOR MACHINERY, NOT BESIDE IT.
  // An earlier version of this class carried its own Breit-Wigner,
  // its own breakup momentum, its own Blatt-Weisskopf barriers and
  // its own Gounaris-Sakurai routine. All four already exist and are
  // better: Propagator.H supplies BreitWigner (fixed/running/GS),
  // Summed_Propagator and Multiplied_Propagator, and the running
  // widths come from Total_Width_Base via the global LineShapes
  // registry, i.e. from the resonance's ACTUAL decay channels rather
  // than from a single-channel approximation. Reimplementing them
  // meant the K*(892) here had a different width function from the
  // K*(892) everywhere else in the library. This class now only
  // assembles propagators; it computes no lineshape of its own.
  //
  // What remains genuinely local is the LASS background phase, which
  // is a property of elastic K pi scattering rather than of any one
  // resonance and has no counterpart in Propagator.H.
  //
  // Index convention as in FF_0_PP, NOT FF_P_X: m_pi[0], m_pi[1] are
  // the two pseudoscalars and s is built from their SUM.
  //
  //////////////////////////////////////////////////////////////////////////

  enum class pp_wave_model {
    /// Constant - pure phase space. A null test, nothing more.
    flat          = 0,
    /// One resonance, via BreitWigner or Complex_Pole (see _TYPE).
    resonance     = 1,
    /// Coherent sum of two resonances with a relative magnitude and
    /// phase, assembled with Summed_Propagator. The natural
    /// "kappa + K*_0(1430)" or "sigma + f_0(980)" admixture.
    two_resonance = 2,
    /// LASS: elastic effective-range background phase-shifted onto a
    /// resonance propagator. S wave only. This is what the BaBar and
    /// BESIII D -> K pi l nu amplitude analyses fit, so its parameters
    /// transfer from a publication without reinterpretation.
    LASS          = 3,
    unknown       = 999
  };

  class FF_PP_Wave : public FormFactor_Base {
  protected:
    pp_wave_model m_model;
    /// Orbital angular momentum, 0/1/2. Only a consistency label here:
    /// the L-dependence of the running width already lives in the
    /// Partial_Width_Base objects behind each propagator, so applying
    /// a barrier factor again would double count it.
    int    m_L;

    double m_m1, m_m2, m_sthr;

    /// Assembled propagator structure. Owned.
    Propagator_Base * p_prop;
    /// The resonance sitting under the LASS background, if used.
    Propagator_Base * p_lass;

    // LASS background: cot(delta_B) = 1/(a p) + r p / 2.
    double m_a, m_r;
    /// Optional s above which LASS is frozen. It is an ELASTIC
    /// parametrization with no validity above the K eta' threshold;
    /// off by default so that nothing is silently truncated.
    double m_cutoff;

    /// Breakup momentum. Local because the LASS background needs it;
    /// the propagators compute their own internally.
    double Breakup(const double & s) const;
    /// Build one propagator from the YAML keys <name>_KF<i>, _TYPE<i>,
    /// and for a complex pole _MPOLE<i>, _GPOLE<i>.
    Propagator_Base * BuildOne(const std::string & tag);
    Complex LASSAmp(const double & s);

    void ReadParameters();
    void PrintSetup() const;

  public:
    FF_PP_Wave(const FF_Parameters & params);
    ~FF_PP_Wave();

    Complex Value(const double & s);
    Complex operator()(const ATOOLS::Vec4D_Vector& moms);

    inline pp_wave_model Model() const { return m_model; }
    inline int           L()     const { return m_L; }
    inline double        Threshold() const { return m_sthr; }
    /// The dominant propagator, so that the owning current can read
    /// its pole mass and on-shell width to derive the decay coupling.
    inline Propagator_Base * Prop() const {
      return (p_prop!=NULL ? p_prop : p_lass);
    }
  };
}

#endif
