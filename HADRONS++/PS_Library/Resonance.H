#ifndef HADRONS_PS_Library_Resonance_H
#define HADRONS_PS_Library_Resonance_H

#include "ATOOLS/Org/Getter_Function.H"
#include "ATOOLS/Phys/Flavour.H"
#include <list>
#include <string>
#include <vector>

namespace HADRONS {
  enum class resonance_type {
    fixed   = 0,
    running = 1,
    GS      = 2,
    bespoke = 10
  };

  struct Resonance_Parameters {
    ATOOLS::Flavour              m_inflav;
    std::vector<ATOOLS::Flavour> m_outflavs;
    resonance_type               m_type;
    double m_OSmass, m_OSwidth, m_phase;
    
    Resonance_Parameters(const ATOOLS::Flavour & inflav,
			 const std::vector<ATOOLS::Flavour> & outflavs,
			 const resonance_type & type=resonance_type::running,
			 const double & OSmass=-1.,const double & OSwidth=-1.,
			 const double & phase=0.
			 ) :
    m_inflav(inflav), m_outflavs(outflavs), m_type(type),
    m_OSmass(OSmass), m_OSwidth(OSwidth), m_phase(phase) {}
  };
  
  typedef std::list<Resonance_Parameters> ResonanceParameters;
  
  
  class Resonance_Base {
  protected:
    ATOOLS::Flavour m_inflav;
    resonance_type  m_type;
    double          m_OSmass, m_OSmass2, m_OSwidth, m_OSwidth2, m_phase;
    double          m_threshold, m_threshold2, m_lambda;
    int             m_exponent;
    std::string     m_name;
    std::vector<double> m_decmasses2;

    virtual double CalculateWidth(const double & s) = 0;
  public:
    Resonance_Base(const Resonance_Parameters & params);    
    virtual const ATOOLS::Flavour & Flav()           { return m_inflav; }
    virtual const double Mass(const double & s=0.)   { return m_OSmass; }
    virtual const double Mass2(const double & s=0.)  { return m_OSmass2; }
    virtual const double Width(const double & s=0.)  { return CalculateWidth(s); }
    virtual const double Width2(const double & s=0.) { return ATOOLS::sqr(Width(s)); }
    virtual const double Phase()                     { return m_phase; }
    
    inline double Lambda(const double & M2,const double & m12,const double & m22) {
      return sqrt( ATOOLS::Max(0., ATOOLS::sqr(M2-m12-m22)-4.*m12*m22)/(4.*M2) );
    }
    virtual Complex BreitWigner(const double & s,const double & M2,const double & MW) {
      return M2/Complex( M2-s, -MW );
    }
    virtual Complex AltBreitWigner(const double & s,const double & M2,const double & MW) {
      return s/Complex( M2-s, -MW );
    }
    virtual Complex BreitWigner(const double & s) = 0;
  };

  class FixedWidth_Resonance : public Resonance_Base {
    double m_MGamma;
    double CalculateWidth(const double & s) { return m_OSwidth; }
  public:
    FixedWidth_Resonance(const Resonance_Parameters & params) :
      Resonance_Base(params), m_MGamma(m_OSmass*m_OSwidth) {}
    Complex BreitWigner(const double & s) { return Resonance_Base::BreitWigner(s,m_OSmass2,m_MGamma); }
  };
  
  class RunningWidth_Resonance : public Resonance_Base {
    double CalculateWidth(const double & s);
  public:
    RunningWidth_Resonance(const Resonance_Parameters & params);
    Complex BreitWigner(const double & s); 
  };


  class GS_Resonance : public Resonance_Base {
    double m_m_pi, m_d, m_h0, m_dh0;
    double f(const double & s);
    double h(const double & s);
    double CalculateWidth(const double & s);
  public:
    GS_Resonance(const Resonance_Parameters & params);
    Complex BreitWigner(const double & s); 
  };

}
#endif
