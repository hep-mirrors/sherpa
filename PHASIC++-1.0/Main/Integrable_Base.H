#ifndef Integrable_Base_H
#define Integrable_Base_H

#include "Flavour.H"
#include "Selector.H"
#include "Vector.H"

namespace PHASIC {

  class Integrable_Base {
  private:

    double m_overflow;

  public:
    // constructor
    Integrable_Base();

    // member functions
    virtual int Nin()  = 0;
    virtual int Nout() = 0;
    virtual int Nvec() = 0;
    virtual int Size() = 0;

    virtual std::string      Name()  = 0;
    virtual ATOOLS::Flavour *Flavs() = 0;

    virtual double Differential(ATOOLS::Vec4D *) = 0;
    virtual double Differential2()               = 0;

    virtual double KFactor(double)        = 0;
    virtual double Scale(ATOOLS::Vec4D *) = 0;
    virtual double ISRThreshold()         = 0;
    virtual double Max()                  = 0;

    virtual ATOOLS::Selector_Base *Selector() = 0;
    virtual Integrable_Base       *Selected() = 0;

    virtual void SwapInOrder()    = 0;
    virtual void RestoreInOrder() = 0;
    virtual void AddPoint(double) = 0;

    virtual void DeSelect()     = 0;
    virtual void SelectOne()    = 0;
    virtual void SetMax(double) = 0;

    virtual void SetMomenta(ATOOLS::Vec4D *) = 0;
    virtual void SetHelicitySampling(bool);
    
    // inline functions
    inline void   SetOverflow(double overflow) { m_overflow=overflow; }
    inline double Overflow()                   { return m_overflow;   }
    
    virtual long int Points() = 0;
    virtual double   Total() = 0;
    virtual double   TotalError() = 0;

  };// end of class Integrable_Base

}// end of namespace PHASIC

#endif


