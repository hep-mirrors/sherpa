#ifndef Monochromatic_h
#define Monochromatic_h

#include "Beam_Base.H"


namespace BEAM {
  class Monochromatic : public Beam_Base {
  public:
    Monochromatic(const APHYTOOLS::Flavour,const double,const double,const int);
    Beam_Base * Copy();     
    bool        CalculateWeight(const double,const double);
    double      Weight(APHYTOOLS::Flavour);
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
    \todo
    This is still to be implemented !
  */
  /*!
    \fn Monochromatic::Monochromatic(const APHYTOOLS::Flavour,const double,const double,const int)
    This is the constructor of the class. here, apart from setting internal variables,
    also the lab vector is constructed. It reads
    \f[
    \vec p_{\rm in} = (E,0,0,\pm \sqrt{E^2-m^2})^T\,,
    \f]
    where the sign of the \f z\f-component of the vector is given by the sign of the energy
    that is handed over. Of course, internally, the energies are all positive.
  */
}
#endif
