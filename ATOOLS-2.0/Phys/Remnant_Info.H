#ifndef Remnant_Info_H
#define Remnant_Info_H

namespace ATOOLS {

  class Flavour;

  class Remnant_Info {
  public:

    virtual double MinimalEnergy(const ATOOLS::Flavour &flavour) = 0;

  };// end of class Remnant_Info

}// end of namespace ATOOLS

#endif
