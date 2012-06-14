#ifndef PDF_Main_Intact_H
#define PDF_Main_Intact_H

#include "PDF/Main/ISR_Base.H"

namespace PDF {

  class Intact : public ISR_Base {
  public:

    // constructor
    Intact(ATOOLS::Flavour);     

    // member functions
    bool   CalculateWeight(double x,double z,double kp2,double Q2);
    double Weight(ATOOLS::Flavour flavour);

  };// end of class Intact

  /*!
    \file 
    \brief contains the class PDF::Intact.
  */
  /*!
    \class PDF::Intact
    \brief For the handling of no initial state radiation.

    This class is used as a kind of "dummy" wrapper for particles experiencing no
    initial state radiation. Therefore the internal structure is pretty simple,
    and there is no PDF_Base related to it, quite in contrast to the Structure_Function.

    Its type is "(None)".
  */
  /*!
    \fn Intact::Intact(Flavour)
    Sets the flavour of the bunch and the type, fixes the weight to 1.
  */
  /*!
    \fn bool Intact::CalculateWeight(const double,const double)
    No weight has to be calculated here.
  */
  /*!
    \fn double Intact::Weight(ATOOLS::Flavour)
    The trivial weight is returned.
  */

}// end of namespace PDF

#endif



