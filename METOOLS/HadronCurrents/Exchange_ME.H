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
  //
  // EACH EXCHANGE IS ONE DIAGRAM, with its own q. That stops being a
  // formality as soon as a real photon is attached to a current: in
  // e+e- -> mu+mu- gamma the initial-state diagram has the photon in J1 and
  // q^2 = (p0+p1-k)^2, the final-state one has it in J2 and q^2 = s. Both
  // are needed, so q cannot belong to the matrix element as a whole.
  //
  // A photon is declared in one of two ways:
  //
  //   explicitly, by listing its leg in the Indices of the current it is
  //   emitted from. The diagram is then exactly as written -- one entry per
  //   placement -- which is how to ask for, say, initial-state radiation
  //   alone.
  //
  //   or with a block-level key
  //
  //     Currents:
  //       Photon: 4                       (or a list, [4, 5])
  //       Exchange: [ ...the photon-less diagrams... ]
  //
  //   which expands every entry into one diagram per charged line, the
  //   photon's leg appended to that current's Indices. Complete by
  //   construction.
  //
  // The coupling at the emission vertex is the model's alpha_QED, the same
  // one the photon exchange defaults to, so a pure-QED amplitude is a
  // consistent e^3. Photon_Alpha: on the block overrides it -- e.g. alpha(0)
  // for a real photon while the exchange runs at the hard scale.
  //
  // Nothing in here knows what a photon DOES. Emitting one is the current's
  // job; this class only decides which current gets the leg and derives q
  // from J1's legs exactly as before, which gives p0+p1-k automatically.
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
      // This diagram's legs. Per exchange, NOT per matrix element: a photon
      // in J1 or J2 changes q. Photon legs included.
      std::vector<int>    m_i1, m_i2;
      std::vector<double> m_qsign;  // per J1 leg: -1 incoming, +1 outgoing
      Exchange(): p_prop(NULL), p_j1(NULL), p_j2(NULL), m_g(0.), m_M2(0.) {}
    };

  private:

    // One Exchange: entry with its legs resolved -- after any Photon:
    // expansion, before any current is built. Every structural check runs
    // on these, so a malformed block is refused before a current can throw
    // something less specific.
    struct Diagram {
      ATOOLS::Scoped_Settings m_e;
      std::vector<int>        m_i1, m_i2;
      int                     m_kf;     // 0 for a contact term
      Diagram(ATOOLS::Scoped_Settings e): m_e(e), m_kf(0) {}
    };

    std::vector<Exchange> m_x;
    std::vector<int>      m_line1, m_line2; // J1/J2 legs, photons removed
    std::vector<int>      m_photons;        // external photon legs
    double                m_photonalpha;    // alpha at an emission vertex
    std::vector<size_t>   m_nspin;      // spin states per leg, 1 for a scalar
    size_t                m_nhel;
    size_t                m_nin;        // 1 for a decay, 2 for scattering
    int                   m_chan;       // 1 = s, 2 = t, 3 = u, 0 = neither

    Current_Base * MakeCurrent(ATOOLS::Scoped_Settings s,
                               const std::vector<int> & idx);
    void AddExchange(const Diagram & d);

    bool IsPhotonLeg(const int i) const;
    std::vector<int> StripPhotons(const std::vector<int> & idx) const;
    /// Whether any leg on this line carries electric charge.
    bool IsChargedLine(const std::vector<int> & line) const;
    /// sum over the line's legs of +Q outgoing, -Q incoming.
    double ChargeFlow(const std::vector<int> & line) const;

    /// Read the Exchange: list, expanding by Photon: if given.
    std::vector<Diagram> ReadDiagrams(ATOOLS::Scoped_Settings & s);
    /// Every check that can be made on the legs alone. See the definition.
    void CheckDiagrams(const std::vector<Diagram> & ds);
    /// Whether the photon placements are gauge invariant. See the definition.
    void CheckPhotonCoverage(const std::vector<Diagram> & ds) const;

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

    // The momentum one exchange carries, from its J1's legs alone.
    ATOOLS::Vec4D Q(const Exchange & x,const ATOOLS::Vec4D_Vector& moms) const;

    // Which Mandelstam the exchange carries, for a 2 -> 2: 1 = s, 2 = t,
    // 3 = u, 0 for anything else -- including every process with a photon,
    // which is never a 2 -> 2. Classified by WHICH TWO LEGS J1 PAIRS, never
    // by whether it happens to contain leg 2 or 3 -- the latter is only
    // right while J1 also holds leg 0, which stops being true as soon as the
    // initial state is reordered (Sherpa sorts heavier first once the
    // masses differ).
    inline int Channel() const { return m_chan; }

    // Electroweak order in the sense of a run card's Order: block (powers of
    // the coupling in the amplitude): 2 for the exchange -- G_F counts as
    // g^2/M_W^2 for a contact term -- plus one per photon.
    inline int OrderEW() const { return 2+int(m_photons.size()); }

    inline const std::vector<Exchange> & Exchanges() const { return m_x; }
    // Non-const, so a caller can retune a built exchange -- in particular
    // switch one off by setting its coupling to zero, which is how EXTRA_XS
    // isolates pure QED for its analytic cross check. Zeroing rather than
    // erasing, because the list owns its currents and propagators.
    inline std::vector<Exchange> & Exchanges() { return m_x; }
    /// Whether external leg i is incoming. See m_nin.
    inline bool IsIncomingLeg(const int i) const { return i<int(m_nin); }

    // The two lines every diagram joins: J1's and J2's legs with any photon
    // removed. Identical for every exchange -- CheckDiagrams enforces it --
    // which is what lets a caller assign colour once.
    inline const std::vector<int> & Line1() const { return m_line1; }
    inline const std::vector<int> & Line2() const { return m_line2; }
    inline const std::vector<int> & Photons() const { return m_photons; }
    inline double PhotonAlpha() const { return m_photonalpha; }

    /// Test scaffold for XS_CURRENTS_GAUGE: every current replaces eps* by
    /// k, so a gauge-invariant amplitude must vanish. Off in production.
    void SetPhotonGaugeTest(const bool on);
  };
}

#endif
