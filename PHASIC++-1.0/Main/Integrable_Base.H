#ifndef Integrable_Base_H
#define Integrable_Base_H

#include "Flavour.H"
#include "Cut_Data.H"
#include "Selector.H"
#include "Vector.H"

namespace PHASIC {
  class Integrable_Base;

  class Integrable_Base {
  public:
    Integrable_Base();
    virtual int                        Nin()                           = 0;
    virtual int                        Nout()                          = 0;
    virtual int                        Nvec()                          = 0;
    virtual std::string                Name()                          = 0;
    virtual ATOOLS::Flavour          * Flavs()                         = 0;

    virtual double                     Differential(ATOOLS::Vec4D *)   = 0;
    virtual double                     Differential2()                 = 0;
    virtual double                     KFactor(double)                 = 0;
    virtual double                     Scale(ATOOLS::Vec4D *)          = 0;
    virtual ATOOLS::Cut_Data         * Cuts()                          = 0;
    virtual ATOOLS::Selector_Base    * Selector()                      = 0;
    virtual void                       UpdateCuts(double,double)       = 0;
    virtual double                     ISRThreshold()                  = 0;
    virtual double                     Max()                           = 0;
    virtual void                       AddPoint(double)                = 0;
    virtual void                       DeSelect()                      = 0;
    virtual void                       SelectOne()                     = 0;
    virtual void                       SetMax(double)                  = 0;
    virtual void                       SetMomenta(ATOOLS::Vec4D *)     = 0;
    virtual Integrable_Base          * Selected()                      = 0;
    virtual void                       SetHelicitySampling(bool) {}
    virtual void                       SwapInOrder()                   = 0;
    virtual void                       RestoreInOrder()                = 0;
    void SetOverflow(double);
    double Overflow();
  private:
    double m_overflow;
  };

  inline void Integrable_Base::SetOverflow(double overflow) { m_overflow=overflow; }
  inline double Integrable_Base::Overflow() { return m_overflow; }

}

#endif


