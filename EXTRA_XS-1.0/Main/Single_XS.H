#ifndef Single_XS_H
#define Single_XS_H

#include "XS_Base.H"
#include "Standard_Selector.H"

namespace EXTRAXS {
  class Single_XS : public XS_Base {
  protected:
    double               thres;
    virtual void GenerateName();
  public:
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
    Single_XS(int,int,APHYTOOLS::Flavour *);
    virtual void SetISR(ISR::ISR_Handler *);
    virtual void SetBeam(BEAM::Beam_Handler *);
    virtual double Differential(AMATOOLS::vec4d *);
    virtual double Differential2();
    virtual double DSigma(double,double,double);
    virtual double DSigma2();
    virtual double operator() (double,double,double);
    virtual XS_Base * Selected() { return this; }
  };
}



#endif



