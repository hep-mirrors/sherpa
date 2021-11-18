#ifndef BEAM_Main_Monochromatic_H
#define BEAM_Main_Monochromatic_H

#include "BEAM/Main/Beam_Base.H"


namespace BEAM {
  class Monochromatic : public Beam_Base {
  public:
    Monochromatic(const ATOOLS::Flavour,const double,const double,const int);
    Beam_Base * Copy();     
    ATOOLS::Flavour Remnant();
    bool        CalculateWeight(const double,const double);
    double      Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none));
    bool        On() { return 0; }
  };
  /*!
    \file 
    \brief contains the class BEAM::Monochromatic.
  */
  /*!
    \class Monochromatic
    \brief This class is for no beam spectra, i.e. monochromatic beams.
    In the unpolarized case, the weight will be exactly one, otherwise it should be
    given by
    \f[
    {\cal W} = \frac{\lambda_{\rm in}+\lambda_{\rm out}}{2}\,,
    \f]
    where both \f$\lambda_{\rm in,out}$\f are the polarization degrees of the incoming
    and outgoing particle.
    \todo This is still to be implemented !
  */
  /*!
    \fn Monochromatic::Monochromatic(const ATOOLS::Flavour,const double,const double,const int)
    This is the constructor of the class. Vectors are constructed through a call to the 
    constructor of the Beam_Base.
  */
  /*!
    \fn Beam_Base * Monochromatic::Copy()
    Returns an exact copy of this class.
  */  
  /*!
    \fn 
    bool Monochromatic::CalculateWeight(const double,const double)
    No weight is calculated so far.
  */  
  /*!
    \fn double Monochromatic::Weight(ATOOLS::Flavour=ATOOLS::Flavour(kf_none))
    The weight is returned. It has been initialised in the constructor.
  */  
}
#endif
