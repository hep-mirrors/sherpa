#ifndef METOOLS_HadronCurrents_Exchange_ME_H
#define METOOLS_HadronCurrents_Exchange_ME_H

#include "METOOLS/HadronCurrents/ME_Base.H"
#include "METOOLS/HadronCurrents/Current_Base.H"
#include "METOOLS/HadronCurrents/FormFactors/Propagator.H"

namespace METOOLS {

  // One amplitude, built as a sum over exchanged bosons each joining two
  // currents:
  //
  //   M = sum_V  g_V  D_V(q^2) [ J1_V . J2_V - (J1_V.q)(J2_V.q)/M_V^2 ]
  //
  // declared as
  //
  //   Currents:
  //     Exchange:
  //     - Boson: 22                                  (omit for a contact term)
  //       J1: {Type: VA_F_F, Indices: [1, 0], ...}
  //       J2: {Type: VA_P_P, Indices: [0, 1], ...}
  //
  // and used from BOTH sides of the code: a hadron decay channel reads the
  // same block as a process declaration in a run card. This class only ever
  // builds the amplitude and stores it per helicity; what happens next --
  // squaring and averaging for a cross section, or contracting with a
  // production density matrix for a decay -- belongs to the caller.
  //
  // q is DERIVED from J1's legs, with the sign that makes each one outgoing,
  // so nothing here needs telling whether an exchange is s-, t- or u-like,
  // nor whether it is looking at a decay or a scattering process.
  //
  // Every key on a current other than Type and Indices is passed through to
  // that current verbatim, which is how form-factor models and their fitted
  // parameters arrive.
  class Exchange_ME : public ME_Base {
  public:

    // One exchanged boson. It owns its OWN pair of currents, because the
    // couplings belong to the (boson, current) pair and not to the current:
    // a photon wants (v,a) = (Q,0) on the same line where a Z wants
    // (v_e, a_e), and one current object cannot be both.
    struct Exchange {
      Propagator_Base * p_prop;   // NULL for a contact term
      Current_Base * p_j1;
      Current_Base * p_j2;
      double            m_g;      // overall coupling, both vertices
      double            m_M2;     // 0 => no q^mu q^nu term
      Exchange(): p_prop(NULL), p_j1(NULL), p_j2(NULL), m_g(0.), m_M2(0.) {}
    };

  private:

    std::vector<Exchange> m_x;
    std::vector<int>      m_i1, m_i2;   // leg lists of J1 and J2
    std::vector<double>   m_qsign;      // -1 incoming, +1 outgoing
    std::vector<size_t>   m_nspin;      // spin states per leg, 1 for a scalar
    size_t                m_nhel;
    size_t                m_nin;        // 1 for a decay, 2 for scattering
    int                   m_chan;       // 1 = s, 2 = t, 3 = u, 0 = neither

    Current_Base * MakeCurrent(ATOOLS::Scoped_Settings s,
                               const std::vector<int> & idx);
    void AddExchange(ATOOLS::Scoped_Settings e);

    /// Reject a flavour-changing neutral current. See the definition.
    void CheckFlavourFlow(const bool contact,const ATOOLS::Flavour & bos,
                          const std::vector<int> & i1,
                          const std::vector<int> & i2) const;

  public:

    /// nin is how many external legs are incoming: 1 for a decay, where
    /// leg 0 is the decaying particle, 2 for a scattering process. It is
    /// not derivable from the flavours, so the caller states it; the decay
    /// default keeps the HADRONS++ path unchanged.
    Exchange_ME(const ATOOLS::Flavour_Vector& flavs,
                const std::vector<int>& indices,
                const std::string& name,
                const size_t nin=1);
    ~Exchange_ME();

    void SetModelParameters(ATOOLS::Scoped_Settings& s);
    void Calculate(const ATOOLS::Vec4D_Vector& moms, bool anti=false);

    // The momentum the exchange carries, from J1's legs alone.
    ATOOLS::Vec4D Q(const ATOOLS::Vec4D_Vector& moms) const;

    // Which Mandelstam that is, for a 2 -> 2: 1 = s, 2 = t, 3 = u, 0 for
    // anything else. Classified by WHICH TWO LEGS J1 PAIRS, never by whether
    // it happens to contain leg 2 or 3 -- the latter is only right while J1
    // also holds leg 0, which stops being true as soon as the initial state
    // is reordered (Sherpa sorts heavier first once the masses differ).
    inline int Channel() const { return m_chan; }

    inline const std::vector<Exchange> & Exchanges() const { return m_x; }
    // Non-const, so a caller can retune a built exchange -- in particular
    // switch one off by setting its coupling to zero, which is how EXTRA_XS
    // isolates pure QED for its analytic cross check. Zeroing rather than
    // erasing, because the list owns its currents and propagators.
    inline std::vector<Exchange> & Exchanges() { return m_x; }
    /// Whether external leg i is incoming. See m_nin.
    inline bool IsIncomingLeg(const int i) const { return i<int(m_nin); }

    inline const std::vector<int> & J1Indices() const { return m_i1; }
    inline const std::vector<int> & J2Indices() const { return m_i2; }
  };
}

#endif
