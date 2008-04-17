#ifndef Splitting_Kernel_H
#define Splitting_Kernel_H

#include "Flavour.H"

namespace PDF {

  class Splitting_Kernel {
  private:
    
    ATOOLS::Flavour m_fla, m_flb, m_flc;

  public:

    // constructor
    Splitting_Kernel(const ATOOLS::Flavour &a,
		     const ATOOLS::Flavour &b,
		     const ATOOLS::Flavour &c);

    // destructor
    virtual ~Splitting_Kernel();

    // member functions
    virtual double operator()(const double &z) const = 0;

    virtual double Value(const double &z) const;
    virtual double Integral(const double &zmin,const double &zmax) const = 0;

    // inline functions
    inline ATOOLS::Flavour GetA() const { return m_fla; }
    inline ATOOLS::Flavour GetB() const { return m_flb; }
    inline ATOOLS::Flavour GetC() const { return m_flc; }
    
  };// end of class Splitting_Kernel

}// end of namespace PDF

#endif
