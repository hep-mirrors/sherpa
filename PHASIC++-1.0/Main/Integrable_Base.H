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
    inline virtual int                        Nin()                           = 0;
    inline virtual int                        Nout()                          = 0;
    inline virtual int                        Nvec()                          = 0;
    inline virtual std::string                Name()                          = 0;
    inline virtual APHYTOOLS::Flavour       * Flavs()                         = 0;

    inline virtual double                     Differential(AMATOOLS::Vec4D *) = 0;
    inline virtual double                     Differential2()                 = 0;
    inline virtual double                     KFactor(double)                 = 0;
    inline virtual double                     Scale(AMATOOLS::Vec4D *)        = 0;
    inline virtual APHYTOOLS::Cut_Data      * Cuts()                          = 0;
    inline virtual APHYTOOLS::Selector_Base * Selector()                      = 0;
    inline virtual void                       UpdateCuts(double,double)       = 0;
    inline virtual double                     ISRThreshold()                  = 0;
    inline virtual double                     Max()                           = 0;
    inline virtual void                       AddPoint(double)                = 0;
    inline virtual void                       DeSelect()                      = 0;
    inline virtual void                       SelectOne()                     = 0;
    inline virtual void                       SetMax(double)                  = 0;
    inline virtual void                       SetMomenta(AMATOOLS::Vec4D *)   = 0;
    inline virtual Integrable_Base          * Selected()                      = 0;
  };
}

#endif


