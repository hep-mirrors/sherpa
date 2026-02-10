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
    double             m_M, m_M2;
  public:
    Propagator_Base(Total_Width_Base * width,
		    const resonance_type & type=resonance_type::running);
    virtual const Complex operator()(const double & s) = 0;
    virtual const Complex Normalised(const double & s) = 0;
    virtual const double Normalised2(const double & s) = 0;
  };

  class BreitWigner : public Propagator_Base {
  public:
    BreitWigner(Total_Width_Base * width,
		const resonance_type & type=resonance_type::running) :
      Propagator_Base(width,type) {}
    const Complex operator()(const double & s);
    const Complex Normalised(const double & s);
    const double Normalised2(const double & s);
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
    const double Normalised2(const double & s);
    inline std::map<Propagator_Base *, Complex> & GetAll() { return m_props; }
  };
}

#endif
