#ifndef Beam_Base_h
#define Beam_Base_h

#include "Flavour.H"
#include "Vector.H"
#include "Message.H"

namespace BEAM {
  class Beam_Base { 
  protected:
    std::string        m_type;
    APHYTOOLS::Flavour m_beam;
    double             m_energy,m_polarization;
    AMATOOLS::Vec4D    m_lab,m_vecout;
    int                m_dir;
    double             m_weight,m_x,m_Q2;
  public:
    Beam_Base(std::string _type,const APHYTOOLS::Flavour _beam,
	      const double _energy,const double _polarization,const int _dir) :
      m_type(_type), m_beam(_beam), m_energy(_energy),m_polarization(_polarization),
      m_dir(_dir),m_x(1.),m_Q2(0.),m_weight(1.)
    {
      double disc      =  1.-AMATOOLS::sqr(m_beam.PSMass()/m_energy);
      if (disc<0) {
	AORGTOOLS::msg.Error()<<"Error in Beam_Base :"<<m_type<<std::endl
			      <<"   Mismatch of energy and mass of beam particle : "
			      <<m_beam<<" / "<<m_energy<<std::endl
			      <<"   Will lead to termination of program."<<std::endl;
	abort();
      }
      m_lab    = AMATOOLS::Vec4D(m_energy,0.,0.,_dir*m_energy*sqrt(disc));
      m_vecout = AMATOOLS::Vec4D(m_energy,0.,0.,_dir*m_energy*sqrt(disc));
    }

    virtual ~Beam_Base() {};
    virtual Beam_Base * Copy()                                     = 0;
    virtual bool        CalculateWeight(const double,const double) = 0;

    virtual APHYTOOLS::Flavour Flav()                     { return m_beam; }
    virtual AMATOOLS::Vec4D    OutMomentum()              { return m_vecout; }
    virtual AMATOOLS::Vec4D    Momentum()                 { return m_lab; }
    virtual double             Energy()                   { return m_energy; }
    virtual double             Polarization()             { return m_polarization; }
    virtual double             Weight()                   { return m_weight; }
    virtual double             Weight(APHYTOOLS::Flavour) { return m_weight; }
    virtual bool               On()                       { return 0; }
    virtual std::string        Type()                     { return m_type; }
    virtual double             Exponent()                 { return 0.; }
    virtual double             Xmax()                     { return 1.; }
    virtual double             Peak()                     { return 1.; }
  };

  /*!
    \file 
    \brief contains the class BEAM::Beam_Base.
  */
  /*!
    \class Beam_Base
    \brief A mainly virtual class for different beam spectra.
    
    This class is the base class for all possible beam spectra that are implemented and 
    for interfaces to new beam spectra provided from outside the program. So far, two 
    spectra have been implemented, one for monochromatic beams BEAM::Monochromatic,
    and one for Laser back scattering of an electron BEAM::Laser_Backscattering.

    In the corresponding constructors, the four momenta of the incoming beams in the lab 
    frame are constructed. They are important especially for asymmetric collisions, like,
    e.g., at HERA.

    \todo
    A number of beam spectra are still missing : 
    Gaussian, Beamstrahlung and Weiszaecker-Williams off both electron and protons. 
  */
  /*!
    \var std::string Beam_Base::m_type
    The type of the Beam_Base, either "Monochromatic", "LaserBackScattering",
    "Gaussian", or "Beamstrahlung". The latter two have not been implemented yet.
    The handling of different beam types is closely related to the class Beam_Type. 
    \sa AORGTOOLS::Beam_Type
  */
  /*!
    \var APHYTOOLS::Flavour Beam_Base::m_beam
    The flavour of the incoming beam particle.
  */
  /*!
    \var double Beam_Base::m_energy
    Energy of the incoming beam particle.
  */
  /*!
    \var double Beam_Base::m_polarization
    The polarization degree of the incoming particle. It should
    range between -1 and 1 and it should - preferably - be a helicity.
  */
  /*!
    \var AMATOOLS::Vec4D Beam_Base::m_lab
    The incoming momentum of the beam particle before the effect of the spectrum kicks in.
  */
  /*!
    \var double Beam_Base::m_weight
    The internal weight for the spectrum. It is calculated through the method
    CalculateWeight(const double,const double).
  */
  /*!
    \fn Beam_Base * Beam_Base::Copy()
    This method produces an exact copy of the beam base in question. All internal parameters -
    excluding the actual weight, m_weight, are copied.
  */
  /*!
    \fn bool Beam_Base::CalculateWeight(const double,const double)
    This is the method that will calculate the weight for a given beam spectrum in
    dependence on the energy fraction of the outgoing particle w.r.t. the incoming one
    and in dependence on the scale of the process. The scale might, for instance, become
    important once we consider Weizsaecker-Williams type spectra.
  */
  /*!
    \fn double Beam_Base::Weight()
    Returns the weight, m_weight.
  */
  /*!
    \fn double Beam_Base::Weight(APHYTOOLS::Flavour)
    Returns the weight after the flavours are checked.
  */
  /*!
    \fn bool Beam_Base::On()
    Returns .false. if no beam-handling is necessary, i.e. for monochromatic beams and
    .true. in all other cases.
  */
  /*!
    \fn double Beam_Base::Exponent() 
    This is a helper for the phase space integration. It returns a characteristic value -
    usually an exponent - for simple channel structures. 
  */
  /*!
    \fn double Beam_Base::Xmax()
    The maximally accessible energy fraction an outgoing particle can have.
  */
  /*!
    \fn double Beam_Base::Peak()
    The peak position for spectra. in most cases this will be 0 or 1; in contrast, for
    laser back scattering it has a different value. This is also needed for the construction
    of integration channels.
  */
}
#endif
