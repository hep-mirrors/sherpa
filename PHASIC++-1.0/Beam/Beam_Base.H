#ifndef Beam_Base_h
#define Beam_Base_h

#include "Flavour.H"
#include "Pol_Info.H"
#include "Message.H"

namespace BEAM {
  class Beam_Base { 
  protected:
    double weight;
    /*!
      The type of the Beam_Base, either "None" or "LaserBackScattering" or "Beamstrahlung".
    */
    std::string  type;
  public:
    virtual ~Beam_Base()          {};
    virtual bool   CalculateWeight(double,double) {
      AORGTOOLS::msg.Error()<<" Virtual Method Beam_Base::CalculateWeight(). "<<std::endl; 
    }
    virtual double Weight()                   { return weight; }
    virtual double Weight(APHYTOOLS::Flavour) { return weight; }
    virtual bool            On()              { return 0; }
    virtual std::string     Type()            { return type; }
    virtual double          Exponent()        { return 0.; }
    virtual double          Xmax()            { return 1.; }
    virtual double          Peak()            { return 1.; }
  };
}
#endif
