#ifndef METOOLS_HadronCurrents_Propagator_H
#define METOOLS_HadronCurrents_Propagator_H


#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"
#include <map>
#include <vector>

namespace METOOLS {
  enum class resonance_type {
    fixed    = 0,
    running  = 1,
    GS       = 2,
    bespoke  = 99
  };

  class Total_Width_Base;

  class Propagator_Base {
  protected:
    Total_Width_Base * p_width;
    resonance_type     m_type;
    double             m_M, m_M2, m_Gamma;

    // Width entering the propagator denominator: the running width from the
    // lineshape, or the fixed on-shell width if there is none to run with.
    const double Width(const double & s) const;
  public:
    // M and Gamma default to the particle table entry behind 'width'; pass
    // them explicitly to override it.  Form factors need that override: their
    // amplitudes and phases come out of a fit that determined the resonance
    // parameters at the same time, so the two sets have to be used together
    // and must not be mixed with the PDG values in the particle table.
    Propagator_Base(Total_Width_Base * width,
		    const resonance_type & type=resonance_type::running,
		    const double & M=-1.,const double & Gamma=-1.);
    virtual ~Propagator_Base() {}
    virtual const Complex operator()(const double & s) = 0;
    virtual const Complex Normalised(const double & s) = 0;
    virtual const double Normalised2(const double & s) = 0;
  };

  class Unity : public Propagator_Base {
  public:
    Unity(Total_Width_Base * width=NULL,
	  const resonance_type & type=resonance_type::running) :
      Propagator_Base(width,type) {}
    const Complex operator()(const double & s)  { return Complex(1.,0.); }
    const Complex Normalised(const double & s)  { return Complex(1.,0.); }
    const double  Normalised2(const double & s) { return 1.; }
  };

  class BreitWigner : public Propagator_Base {
  public:
    BreitWigner(Total_Width_Base * width,
		const resonance_type & type=resonance_type::running) :
      Propagator_Base(width,type) {}
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  class WeightedBreitWigner : public Propagator_Base {
  public:
    WeightedBreitWigner(Total_Width_Base * width,
			const resonance_type & type=resonance_type::running,
			const double & M=-1.,const double & Gamma=-1.) :
      Propagator_Base(width,type,M,Gamma) {}
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // m^2/(m^2-s-i m Gamma).  BreitWigner always runs its width off the
  // lineshape, which is not an option for a resonance that has none, so this
  // is the plain fixed-width form taking its parameters directly.
  class FixedBreitWigner : public Propagator_Base {
  public:
    FixedBreitWigner(Total_Width_Base * width=NULL,
		     const double & M=-1.,const double & Gamma=-1.) :
      Propagator_Base(width,resonance_type::fixed,M,Gamma) {}
    // A fixed width does not switch off below threshold, so evaluating this
    // literally at s = 0 leaves 1/(1 - i Gamma/M) instead of 1 
    const Complex operator()(const double & s)
    { if (std::abs(s)<1.e-12) return Complex(1.,0.);
      return m_M2/Complex(m_M2-s,-m_M*m_Gamma); }
    const Complex Normalised(const double & s)  { return (*this)(s); }
    const double  Normalised2(const double & s) { return norm((*this)(s)); }
  };

  // arXiv:hep-ph/0512180 Eq. (5): 1/(s/m^2 - 1 + i Gamma/m) = -m^2/(m^2-s-i
  // m Gamma). NOT normalised to 1 at s=0 (it evaluates to about -1 there for
  // a narrow resonance) and must be used RAW - do not wrap this in
  // Summed_Propagator(normalise=true), which divides by the sum of weights
  // and silently forces F(0)=1, corrupting this formula's actual magnitude
  // and phase. Found and fixed after a direct 3-way comparison against the
  // tau-decays-in-sherpa fork's own eq-10 implementation showed this
  // codebase's Summed_Propagator-wrapped version badly overshoot the BELLE
  // tail (2-3 GeV) that the fork's raw, unnormalised formula tracks well -
  // see [[hep-ph-0512180-eq10-ef-terms]] memory note.
  class HepPh0512180_BreitWigner : public Propagator_Base {
  public:
    HepPh0512180_BreitWigner(const double & M,const double & Gamma) :
      Propagator_Base(NULL,resonance_type::fixed,M,Gamma) {}
    const Complex operator()(const double & s)
    { return -m_M2/Complex(m_M2-s,-m_M*m_Gamma); }
    const Complex Normalised(const double & s)  { return (*this)(s); }
    const double  Normalised2(const double & s) { return norm((*this)(s)); }
  };

  class GounarisSakurai : public Propagator_Base {
  private:
    double  m_mpi, m_mpi2;
    Complex m_ppiM, m_d, m_hM2, m_dhM2;
    // Everything below the 2 pi threshold is reached by analytic
    // continuation, so these have to be complex.
    const Complex ppi(const double & s) const;
    const Complex h(const double & s) const;
    const Complex dh(const double & s) const;
    const Complex f(const double & s) const;
    // The pi pi width the GS parametrisation is derived with.
    const Complex GammaGS(const double & s) const;
  public:
    GounarisSakurai(Total_Width_Base * width,
		    const resonance_type & type=resonance_type::running,
		    const double & M=-1.,const double & Gamma=-1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // The conformal polynomial C_p of arXiv:1907.01556 eq. (2.14), in the
  // variable z(s) = (sqrt(s_inel-s_1) - sqrt(s_inel-s)) /
  //                 (sqrt(s_inel-s_1) + sqrt(s_inel-s)),
  //
  //     C_p(s) = Sum_{i=1..p} c_i ( z(s)^i - z(0)^i ),
  //
  // so C_p(0) = 0 by construction.  Its role is the large-q^2 behaviour: with
  // the anomaly sum rule imposed on top of it (see the constraint solve in
  // Meson_FormFactors.C) the isoscalar a(q^2) is driven to ZERO at large q^2
  // instead of tending to a constant.  That is what stops the three-pion
  // cross section flattening out above ~2 GeV, where the data falls steeply.
  // Ported from tau-decays-in-sherpa (METOOLS-hadron-currents).
  class Conformal_Polynomial : public Propagator_Base {
  private:
    std::vector<double> m_c;        // c_1 ... c_p; index 0 holds c_1
    double m_sinel, m_s1, m_outer, m_z0;

    Complex Z(const double & s) const;
  public:
    Conformal_Polynomial(const std::vector<double> & c,
                         const double & s_inel=1.,const double & s1=-1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);

    // The limit q^2 -> infinity, where z -> -1.  This is the quantity the
    // anomaly sum rule constrains.
    double AtInfinity() const;

    inline size_t Degree() const { return m_c.size(); }
    inline double Z0() const { return m_z0; }
    inline const std::vector<double> & Coefficients() const { return m_c; }
    void SetCoefficient(const size_t & i,const double & c);
  };

  // c0 + Sum_i w_i ( P_i(s) - P_i(0) ), a once-subtracted sum.
  //
  // This is NOT Summed_Propagator with different weights, and the difference
  // is structural rather than a normalisation choice:
  //   at s = 0 it is exactly c0, so a physical constant - the Wess-Zumino-
  //     Witten anomaly for the three-pion isoscalar - survives to the origin
  //     instead of being divided away;
  //   at large s every P_i dies and it tends to the CONSTANT c0 - Sum_i w_i
  //     P_i(0), which is what holds a form factor up between resonances.  A
  //     tower normalised to 1 at the origin necessarily falls as 1/s instead.
  // Ported from tau-decays-in-sherpa (METOOLS-hadron-currents branch),
  // METOOLS/HadronCurrents/FormFactors/Propagator.{H,C}, which implements
  // arXiv:1907.01556 eq. (2.14).
  class Subtracted_Propagator : public Propagator_Base {
  private:
    struct Term { Propagator_Base * p_prop; Complex m_weight, m_atzero; };
    std::vector<Term> m_terms;
    Complex m_c0;
  public:
    Subtracted_Propagator(const Complex & c0=Complex(1.,0.));
    ~Subtracted_Propagator();
    // Takes ownership of what it is given, so copying it would double-delete.
    Subtracted_Propagator(const Subtracted_Propagator &) = delete;
    Subtracted_Propagator & operator=(const Subtracted_Propagator &) = delete;

    // Takes ownership of prop.
    void  Add(Propagator_Base * prop,const Complex & weight=1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    inline const Complex & Subtraction() const { return m_c0; }
  };

  /*
    The terms of a summed or multiplied propagator, in the order they were
    added. (A std::map keyed by the pointers summed them in heap-address
    order, which changes from run to run, so F differed in the last bits
    between identical runs.)
  */
  typedef std::vector<std::pair<Propagator_Base *, Complex> > Prop_List;

  class Summed_Propagator : public Propagator_Base {
  private:
    Prop_List m_props;
    Complex m_norm;
    bool    m_normalise;
  public:
    // Sums are divided by the sum of their weights by default.  Pass
    // normalise=false for sums that are already normalised by construction,
    // e.g. a "1 + mixing" factor whose mixing term vanishes at s = 0.
    Summed_Propagator(const bool & normalise=true);
    ~Summed_Propagator();
    // Takes ownership of what it is given, so copying it would double-delete.
    Summed_Propagator(const Summed_Propagator &) = delete;
    Summed_Propagator & operator=(const Summed_Propagator &) = delete;

    void  Add(Propagator_Base * prop,const Complex & weight=1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    inline Prop_List & GetAll() { return m_props; }
  };

  class Multiplied_Propagator : public Propagator_Base {
  private:
    Prop_List m_props;
    Complex m_norm;
  public:
    Multiplied_Propagator(Propagator_Base * prop=NULL);
    ~Multiplied_Propagator();
    // Takes ownership of what it is given, so copying it would double-delete.
    Multiplied_Propagator(const Multiplied_Propagator &) = delete;
    Multiplied_Propagator & operator=(const Multiplied_Propagator &) = delete;

    void  Add(Propagator_Base * prop,const Complex & weight=1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    inline Prop_List & GetAll() { return m_props; }
  };
}

#endif
