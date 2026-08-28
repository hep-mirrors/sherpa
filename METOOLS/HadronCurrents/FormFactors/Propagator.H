#ifndef METOOLS_HadronCurrents_Propagator_H
#define METOOLS_HadronCurrents_Propagator_H


#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"
#include <map>
#include <string>

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
    double OnShellWidth() const;
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

  // FM95 Eq.(34): the NORMALIZED, CONSTANT-width Breit-Wigner that
  // hep-ph/9503474 specifies for the K1 envelopes,
  //     BW(s) = (-m^2 + i m G) / (s - m^2 + i m G),
  // identical to Tools::BreitWignerFix(s,m^2,m*G).  This is deliberately NOT
  // the running-width BreitWigner above: FM95 says "normalized Breit-Wigner
  // propagators with constant widths", and a running Gamma(s) grows with s,
  // which inflates the denominator and suppresses the high-mass tail that the
  // K1(1400) is supposed to supply.  Mass and width are passed explicitly
  // rather than taken from a Total_Width_Base, because the paper's Eq.(35)
  // values are parameters of its own fit, not the flavour's PDG pole.
  class FM95_Fixed_BW : public Propagator_Base {
  private:
    double m_G;
  public:
    FM95_Fixed_BW(const double & mass,const double & width);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  class RChL_BW : public Propagator_Base {
  public:
    RChL_BW(Total_Width_Base * width) :
      Propagator_Base(width,resonance_type::running) {}
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // ------------------------------------------------------------------
  // Propagators of the Novosibirsk four-pion parameterisation, as
  // TAUOLA implements them (tauola-factory/Binp/binp.f).
  //
  // These are NOT the BreitWigner above and NOT the paper's Eq.(18).
  // Two differences matter:
  //   * they are NORMALISED to be dimensionless - binp.f divides by
  //     M^2 (and the rho additionally by 1 + (Gamma/M) dm(0)) - with
  //     the dimensions restored by the current's overall 1/M_rho^4/Q;
  //   * the rho carries a Gounaris-Sakurai dispersive real part dm(s).
  //     That is the same physics as BreitWigner's resonance_type::GS,
  //     but binp.f groups and normalises it differently, so reusing
  //     ValueGS() here would not reproduce the published widths.
  //
  // They return D, not 1/D, matching binp.f's z_drho/z_dsigma/z_da1/
  // z_domega - callers divide.  Which resonance is selected by the
  // novo4pi_resonance tag rather than by four near-identical classes.
  // ------------------------------------------------------------------
  enum class novo4pi_resonance { rho, sigma, a1, omega };

  class Novo4Pi_Propagator : public Propagator_Base {
  private:
    novo4pi_resonance m_res;
    double         m_G;        // pole width
    double         m_scaleA;   // a1 form-factor scale (this is Lambda^2)
    Complex        m_z;        // sigma/rho mixing, needed by the a1 width
    double         m_mdau;     // daughter mass, for the sigma phase space
    double         m_rhonorm;  // cached 1/(1 + (Gamma/M) dm(0))
    double         m_g1;       // cached width normalisation at the pole
    bool           m_keepshape;
  public:
    // width may be NULL: these propagators use the parameterisation's
    // own analytic/tabulated widths, not a registered line shape.
    // width is normally NULL: these propagators use the
    // parameterisation's own analytic/tabulated widths.  Passing a
    // registered line shape instead makes the a1 use
    //   D = s/M^2 - 1 + i Gamma(s)/M
    // i.e. Bondar's form with his tabulated g_a1 replaced by Sherpa's
    // running width.  Provided so the two can be compared; see the
    // note in VA_0_Novosibirsk4Pi.H on why the tabulated one is the
    // default.
    Novo4Pi_Propagator(const novo4pi_resonance & res,
                    const double & M,const double & G,
                    const double & mdau=0.13957,
                    const Complex & z=Complex(1.269,0.591),
                    const double & scaleA=1.2,
                    Total_Width_Base * width=NULL,
                    // rho only: keep Bondar's Gounaris-Sakurai real part
                    // and his normalisation, swapping ONLY the width
                    // function for the registered one.  Exists to
                    // isolate which of the three differences between
                    // his propagator and a plain registry one matters.
                    const bool & keepshape=false);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
  };

  // Two-channel Flatte, for a resonance sitting ON a threshold, where a
  // Breit-Wigner - fixed OR running width - cannot work: the channel phase
  // space must go IMAGINARY below threshold, which is the whole point.  The
  // motivating case is K1(1270), whose pole is essentially on the K rho
  // threshold (1.269 vs a 1.27 pole).
  //
  // TRAP, and it is a real one: the Flatte's BARE pole is NOT the PDG
  // Breit-Wigner mass.  For K1(1270) the bare pole is 1.254 while PDG quotes
  // 1.272, and 1.272 sits ABOVE the K rho threshold (~1.264), which INVERTS
  // the sign of the dispersive term that makes the Flatte work.  They are
  // different quantities; never substitute the registry mass here.
  class Two_Channel_Flatte : public Propagator_Base {
  private:
    Complex m_g1, m_g2;
    double m_m11, m_m12, m_m21, m_m22;
    double m_gamma0;

    Complex PhaseSpace(const double & s,const double & m1,
                       const double & m2) const;

  public:
    Two_Channel_Flatte(Total_Width_Base * width,const double & mass,
                       const Complex & g1,const double & m11,
                       const double & m12,
                       const Complex & g2,const double & m21,
                       const double & m22,const double & gamma0=0.);
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double  Normalised2(const double & s);
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
  };
}

#endif
