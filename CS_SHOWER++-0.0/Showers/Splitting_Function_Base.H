#ifndef Splitting_Function_Base_h
#define Splitting_Function_Base_h

#include "Flavour.H"
#include "Vector.H"
#include "Random.H"

namespace CS_SHOWER {
  struct cstp {
    enum code {
      FF   = 11,
      FI   = 12,
      IF   = 21,
      II   = 22,
      none = 0
    };
  };

  class Splitting_Function_Base {
  protected:
    const static double s_Nc;
    const static double s_CF;
    const static double s_CA;
    const static double s_TR;
    const static double s_Tf;
    ATOOLS::Flavour m_flavs[3]; // a -> b c
    // Domain of last z-integral and its value
    double m_zmin,m_zmax,m_lastint;
    cstp::code m_type;
  public:
    Splitting_Function_Base() : m_type(cstp::none) {}
    virtual ~Splitting_Function_Base() {}

    virtual double operator() (const double,const double,
			       const double,const double)      = 0;
    virtual double OverIntegrated(const double,const double,
				  const double,const double)   = 0;
    virtual double Overestimated(const double,const double)    = 0;
    virtual double RejectionWeight(const double,const double,
				   const double,const double)  = 0;
    virtual double Z()                                         = 0;         
    virtual double Phi(double z)            const { return 2.*M_PI*ATOOLS::ran.Get(); }
    virtual double Last()                   const { return m_lastint; } 
    virtual void   SetLast(const double last)     { m_lastint = last; }
    virtual cstp::code        GetType()     const { return m_type; } 
    virtual ATOOLS::Flavour & GetFlavourA()       { return m_flavs[0]; }
    virtual ATOOLS::Flavour & GetFlavourB()       { return m_flavs[1]; }
    virtual ATOOLS::Flavour & GetFlavourC()       { return m_flavs[2]; }

    friend std::ostream& operator<<(std::ostream &,const Splitting_Function_Base &);
  };
}

#endif
