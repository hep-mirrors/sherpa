#ifndef ISR_Base_h
#define ISR_Base_h

#include "PDF_Base.H"
#include "Flavour.H"
#include "Message.H"

namespace ISR {
  class ISR_Base { 
  protected:
    double weight;
    std::string  type;
  public:
    virtual ~ISR_Base()          {};
    virtual bool   CalculateWeight(double,double) {
      AORGTOOLS::msg.Error()<<" Virtual Method ISR_Base::CalculateWeight(). "<<std::endl; 
      return false;
    }
    virtual double Weight(APHYTOOLS::Flavour) =0;
    virtual bool            On()              { return 0; }
    virtual std::string     Type()            { return type; }
    virtual double          Exponent()        { return 0.; }
    virtual double          Xmax()            { return 1.; }
    virtual PDF::PDF_Base * PDF()             { return 0; }
  };
}
#endif
