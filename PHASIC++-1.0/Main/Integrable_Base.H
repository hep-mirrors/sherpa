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
    virtual int                        Nin()                           = 0;
    virtual int                        Nout()                          = 0;
    virtual int                        Nvec()                          = 0;
    virtual std::string                Name()                          = 0;
    virtual APHYTOOLS::Flavour       * Flavs()                         = 0;

    virtual double                     Differential(AMATOOLS::Vec4D *) = 0;
    virtual double                     Differential2()                 = 0;
    virtual double                     KFactor(double)                 = 0;
    virtual double                     Scale(AMATOOLS::Vec4D *)        = 0;
    virtual APHYTOOLS::Cut_Data      * Cuts()                          = 0;
    virtual APHYTOOLS::Selector_Base * Selector()                      = 0;
    virtual void                       UpdateCuts(double,double)       = 0;
    virtual double                     ISRThreshold()                  = 0;
    virtual double                     Max()                           = 0;
    virtual void                       AddPoint(double)                = 0;
    virtual void                       DeSelect()                      = 0;
    virtual void                       SelectOne()                     = 0;
    virtual void                       SetMax(double)                  = 0;
    virtual void                       SetMomenta(AMATOOLS::Vec4D *)   = 0;
    virtual Integrable_Base          * Selected()                      = 0;
  };
}

#endif


