#ifndef Intact_h
#define Intact_h

#include "ISR_Base.H"

namespace PDF {
  class Intact : public ISR_Base {
  public:
    Intact(APHYTOOLS::Flavour);     
    bool       CalculateWeight(double,double);
    double     Weight(APHYTOOLS::Flavour);
    bool       On()           { return 0; }
    PDF_Base * PDF()          { return 0; }
  };
  /*!
    \file 
    \brief contains the class PDF::Intact.
  */
  /*!
    \class PDF::Intact
    \brief For the handling of no initial state radiation.

    This class is used as a kind of dummy wrapper for particles experiencing no
    initial state radiation. Therefore the internal structure is pretty simple,
    and there is no PDF_Base related to it, quite in contrast to the Structure_Function.

    Its type is "(None)".
  */
  /*!
    \fn Intact::Intact(Flavour)
    Sets the flavour of the bunch and the type, fixes the weight to 1.
  */
}
#endif



