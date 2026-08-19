#ifndef METOOLS_HadronCurrents_Propagator_H
#define METOOLS_HadronCurrents_Propagator_H


#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MyComplex.H"

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

  class Summed_Propagator : public Propagator_Base {
  private:
    std::map<Propagator_Base *, Complex> m_props;
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
    inline std::map<Propagator_Base *, Complex> & GetAll() { return m_props; }
  };

  class Multiplied_Propagator : public Propagator_Base {
  private:
    std::map<Propagator_Base *, Complex> m_props;
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
    inline std::map<Propagator_Base *, Complex> & GetAll() { return m_props; }
  };
}

#endif
