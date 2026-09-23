#ifndef METOOLS_HadronCurrents_Propagator_H
#define METOOLS_HadronCurrents_Propagator_H


#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"
#include <map>
#include <string>

namespace METOOLS {
  enum class resonance_type {
    fixed        = 0,
    running      = 1,
    GS           = 2,
    // Pole of the T-matrix on the second Riemann sheet, sqrt(s_pole)
    // = M - i Gamma/2, with NO running width and no real-axis mass.
    // For a state as broad as the f_0(500)/sigma this is not a
    // refinement of the Breit-Wigner, it is the only well-defined
    // option: the sigma pole sits ~275 MeV off the real axis, the
    // "width" exceeds the "mass", and different BW parametrizations
    // of the same data return masses spread over hundreds of MeV.
    // The T-matrix pole is the parametrization-independent statement
    // and is what the PDG actually quotes for these states.
    complex_pole = 3,
    // Two-channel Flatte. The right description for a resonance
    // sitting ON a threshold - f_0(980) and a_0(980) both straddle
    // K Kbar - where a Breit-Wigner is wrong for a structural reason:
    // below the second threshold that channel cannot contribute to the
    // width, but it DOES contribute to the real part of the
    // self-energy, and it is that subthreshold real part which
    // produces the characteristic cusp. A running-width BW has no
    // mechanism for it, and a complex pole has no threshold at all.
    flatte       = 4,
    bespoke      = 99
  };

  class Total_Width_Base;
  
  class Propagator_Base {
  protected:
    Total_Width_Base * p_width;
    resonance_type     m_type;
    double             m_M, m_M2;
  public:
    Propagator_Base(Total_Width_Base * width,
		    const resonance_type & type=resonance_type::running);
    virtual const Complex operator()(const double & s) = 0;
    virtual const Complex Normalised(const double & s) = 0;
    virtual const double Normalised2(const double & s) = 0;
    // Public accessors for diagnostic dumps (request #1) - pole mass
    // and on-shell (i.e. Gamma evaluated at s=M^2) width. Defined out
    // of line in Propagator.C, NOT inline here: Total_Width_Base is
    // only forward-declared in this header (its full definition lives
    // in Resonance_Base.H, included by Propagator.C but not here), so
    // calling its methods requires the complete type to be visible.
    inline double Mass() const { return m_M; }
    // virtual so that propagators without a Total_Width_Base (the
    // complex pole) can still report a meaningful width to the
    // diagnostic dump instead of a silent zero.
    virtual double OnShellWidth() const;
    // Constant numerator this propagator carries on top of the BARE
    // pole 1/(M^2-s-i...). BreitWigner returns M^2/(...), i.e. it is
    // normalised to 1 at s=0, and Complex_Pole returns s_pole/(...).
    // A coupling derived from a decay width assumes the bare pole, so
    // anything combining the two must divide by this. Getting it wrong
    // is a pure M^4 error in the rate and is easy to miss: it left the
    // D -> K pi widths at 0.60 of PDG and D0 -> pi- pi0 at 0.349,
    // which is m_rho^4 = 0.3515 to within half a percent.
    virtual Complex PoleNumerator() const { return Complex(m_M2,0.); }
    ATOOLS::Flavour Flav() const;
  };

  class BreitWigner : public Propagator_Base {
  private:
    // Gounaris-Sakurai specific quantities (only used if m_type==GS).
    // m_mdau2 is the mass^2 of the dominant two-body daughter (e.g. the
    // charged pion for rho->pi pi); it enters the analytic k(s), h(s)
    // functions of the GS parametrization (Gounaris & Sakurai, PRL 21
    // (1968) 244; conventions as used, e.g., in Fujikawa et al.
    // [Belle], Phys.Rev.D78 (2008) 072006, and summarised in
    // arXiv:1509.09140 Eq.(2.2)).
    // All pole quantities (m_GammaPole, m_kPole, m_hPole, m_hPrimePole,
    // m_d) are fixed once at construction time, since they only depend
    // on the (fixed) resonance mass and on-pole width.
    double m_mdau2;
    double m_GammaPole, m_kPole, m_hPole, m_hPrimePole, m_d;

    double k(const double & s) const;
    double h(const double & s) const;
    void   InitGS();
    const Complex ValueGS(const double & s);
  public:
    // mdau<0 (default) means "use the charged pion mass" - the only
    // case currently needed (rho -> pi pi).  Pass an explicit mdau if
    // a GS propagator is ever needed for a different two-body channel.
    BreitWigner(Total_Width_Base * width,
		const resonance_type & type=resonance_type::running,
		const double & mdau=-1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // Breit-Wigner variant matching the literal denominator convention
  // printed in the RChL 3pi current, e.g. Eq.(6) of arXiv:1203.3955:
  // "s - M^2 - i*M*Gamma(s)" (no M^2 numerator - the papers' formulas
  // supply their own numerator factors externally, e.g. the "3*s1" in
  // front of the propagator in Eq.6). This is deliberately a DIFFERENT
  // sign/scaling convention from BreitWigner above (which returns
  // M^2/(M^2-s-i*sqrt(s)*Gamma(s))): the KKpi current (Eq.13-15 of the
  // same paper) instead literally prints "M^2-s-iM*Gamma(s))", i.e. the
  // opposite overall sign from the 3pi current's Eq.6. Rather than
  // "fixing" what may or may not be an inconsistency in the source
  // paper, each channel's code uses whichever convention is literally
  // printed for that channel's own equations - see FF_0_PPP.C's
  // comments at the point of use for which is which.
  class RChL_BW : public Propagator_Base {
  public:
    RChL_BW(Total_Width_Base * width) :
      Propagator_Base(width,resonance_type::running) {}
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // Complex-pole propagator,
  //     D(s) = s_pole / (s_pole - s),   s_pole = (M - i Gamma/2)^2,
  // normalised to D(0) = 1 so that it matches the convention of
  // BreitWigner above (M^2/(M^2-s-i sqrt(s) Gamma(s)) -> 1 at s=0) and
  // can be dropped into a Summed_Propagator alongside one.
  //
  // Takes no Total_Width_Base: there is no running width to evaluate,
  // which is the whole point. The flavour is carried separately, only
  // so that DumpPropagatorStructure can label it.
  class Complex_Pole : public Propagator_Base {
  private:
    Complex         m_spole;
    double          m_Mpole, m_Gpole;
    ATOOLS::Flavour m_flav;
  public:
    Complex_Pole(const double & Mpole,const double & Gpole,
                 const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_none));
    ATOOLS::Flavour Flav() const { return m_flav; }
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    virtual double OnShellWidth() const { return m_Gpole; }
    virtual Complex PoleNumerator() const { return m_spole; }
    inline Complex PolePosition() const { return m_spole; }
  };

  // Flatte, two channels:
  //
  //   D(s) = M^2 / ( M^2 - s - i [ g_1 rho_1(s) + g_2 rho_2(s) ] ),
  //
  // with the two-body phase-space factors
  //   rho_ab(s) = sqrt( (1-(m_a+m_b)^2/s)(1-(m_a-m_b)^2/s) )
  // continued ANALYTICALLY below threshold, where rho becomes
  // i|rho| and the term moves from the imaginary to the real part of
  // the denominator. The couplings g_i carry GeV^2.
  //
  // The M^2 numerator matches the BreitWigner convention in this file
  // (D -> 1 at s=0 up to the subthreshold rho's), so PoleNumerator()
  // returns M^2 and callers deriving a coupling from a width divide it
  // out exactly as they do for a Breit-Wigner.
  class Flatte : public Propagator_Base {
  private:
    double m_g1, m_g2;
    double m_ma1, m_mb1, m_ma2, m_mb2;
    ATOOLS::Flavour m_flav;
    Complex Rho(const double & s,const double & ma,
                const double & mb) const;
  public:
    Flatte(const double & M,const double & g1,const double & g2,
           const double & ma1,const double & mb1,
           const double & ma2,const double & mb2,
           const ATOOLS::Flavour & flav=ATOOLS::Flavour(kf_none));
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    // Width at the pole from the OPEN channel only, so that a caller
    // deriving a decay coupling gets something meaningful.
    virtual double OnShellWidth() const;
    virtual Complex PoleNumerator() const { return Complex(m_M2,0.); }
    ATOOLS::Flavour Flav() const { return m_flav; }
  };

  class Summed_Propagator : public Propagator_Base {
  private:
    std::map<Propagator_Base *, Complex> m_props;
    Complex m_norm;
  public:
    Summed_Propagator(Propagator_Base * prop=NULL);
    ~Summed_Propagator();
    
    void  Add(Propagator_Base * prop,const Complex & weight=1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    inline std::map<Propagator_Base *, Complex> & GetAll() { return m_props; }
    // Ambiguous for a sum - each component carries its own. Returns 1
    // so that a caller dividing by it is left with the components'
    // own normalisations, which is the only defensible default.
    virtual Complex PoleNumerator() const { return Complex(1.,0.); }
  };

  // Diagnostic dump (request #1): print the constituent resonances of
  // a (possibly Summed_) propagator structure - mass, on-shell width,
  // and complex weight - plus the model identifier and a channel
  // label, via msg_Out(). Called once per constructed form factor
  // instance; see the Construct() methods in FF_0_PP.C/FF_0_PPP.C for
  // call sites. Safe to call with props=NULL (e.g. a channel/model
  // combination with no dedicated propagator structure at all).
  void DumpPropagatorStructure(const std::string & label,
				const int & ffmodel_id,
				Propagator_Base * props);

  class Multiplied_Propagator : public Propagator_Base {
  private:
    std::map<Propagator_Base *, Complex> m_props;
    Complex m_norm;
  public:
    Multiplied_Propagator(Propagator_Base * prop=NULL);
    ~Multiplied_Propagator();
    
    void  Add(Propagator_Base * prop,const Complex & weight=1.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
    inline std::map<Propagator_Base *, Complex> & GetAll() { return m_props; }
    // Ambiguous for a sum - each component carries its own. Returns 1
    // so that a caller dividing by it is left with the components'
    // own normalisations, which is the only defensible default.
    virtual Complex PoleNumerator() const { return Complex(1.,0.); }
  };
}

#endif
