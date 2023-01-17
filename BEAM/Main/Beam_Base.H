#ifndef BEAM_Main_Beam_Base_H
#define BEAM_Main_Beam_Base_H

#include "BEAM/Main/Beam_Parameters.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"

namespace BEAM {

  class Beam_Base {
  protected:
    beamspectrum::code m_type;
    ATOOLS::Flavour    m_beam, m_bunch;
    ATOOLS::Vec4D      m_lab, m_vecout, m_position;
    int                m_dir;
    double             m_energy, m_polarisation;
    double             m_x, m_Q2, m_weight;
    bool               m_on;
    
  public:
    Beam_Base(beamspectrum::code type,const ATOOLS::Flavour& _beam,
	      const double _energy,const double _polarisation,
	      const int _dir, int mode = 0);
    virtual ~Beam_Base() = default;

    virtual Beam_Base * Copy()                                     = 0;
    virtual bool        CalculateWeight(const double,const double) = 0;
    virtual void FixPosition() {}
    virtual inline const beamspectrum::code & Type()    const { return m_type; }
    virtual inline const ATOOLS::Flavour    & Beam()    const { return m_beam; }
    virtual inline const ATOOLS::Flavour    & Bunch()   const { return m_bunch; }
    virtual inline const ATOOLS::Flavour    & Remnant() const { return m_beam; }

    virtual inline const ATOOLS::Vec4D & OutMomentum() const { return m_vecout; }
    virtual inline const ATOOLS::Vec4D & InMomentum()  const { return m_lab; }
    virtual inline const ATOOLS::Vec4D & Position()    const { return m_position; }

    virtual inline const double & Energy()       const { return m_energy; }
    virtual inline const double & Polarisation() const { return m_polarisation; }
    virtual inline double Weight(ATOOLS::Flavour fl=kf_none) { return m_weight; }
    virtual inline double Exponent()      const { return 0.; }
    virtual inline double X()             const { return m_x; }
    virtual inline double Xmin()          const { return 0.; }
    virtual inline double Xmax()          const { return 1.; }
    virtual inline double Peak()          const { return 1.; }

    virtual inline bool On()              const { return m_on; }
    virtual inline bool PolarisationOn() const { return (m_polarisation!=0.); }

    virtual inline void SetInMomentum(const ATOOLS::Vec4D &in) { m_lab = in; }
    virtual inline void SetOutMomentum(const ATOOLS::Vec4D &out) { m_vecout = out; }

    inline void SetX(const double &x) { m_x=x; }
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
    \todo
    A number of beam spectra are still missing :
    Gaussian, Beamstrahlung and Weiszaecker-Williams off both electron and protons.
  */
  /*!
    \var std::string Beam_Base::m_type
    The type of the Beam_Base, either "Monochromatic", "LaserBackScattering",
    "Gaussian", or "Beamstrahlung". The latter two have not been implemented yet.
    The handling of different beam types is closely related to the class Beam_Type.
    \sa ATOOLS::Beam_Type
  */
  /*!
    \var ATOOLS::Flavour Beam_Base::m_beam
    The flavour of the incoming beam particle.
  */
  /*!
    \var double Beam_Base::m_energy
    Energy of the incoming beam particle.
  */
  /*!
    \var double Beam_Base::m_polarisation
    The polarization degree of the incoming particle. It should
    range between -1 and 1 and it should - preferably - be a helicity.
  */
  /*!
    \var ATOOLS::Vec4D Beam_Base::m_lab
    The momentum of the incoming beam particle, i.e. before the spectrum becomes effective.
  */
  /*!
    \var ATOOLS::Vec4D Beam_Base::m_vecout
    The momentum of the outgoing bunch particle, i.e. after the spectrum becomes effective.
  */
  /*!
    \var int m_dir
    A direction flag for the construction of the beam momentum.
  */
  /*!
    \var double Beam_Base::m_weight
    The internal weight for the spectrum. It is calculated through the method
    CalculateWeight(const double,const double).
  */
  /*!
    \var double Beam_Base::m_x
    The energy portion the outgoing bunch particle carries away from the beam particle.
    This is in the c.m. frame of the two incoming beams.
  */
  /*!
    \var double Beam_Base::m_Q2
    The scale of the transition of the beams into bunches. This can - in principle - be set
    individually, for instance for Weizsaecker-Williams photons. This option has not been
    fully implemented yet because all beam spectra included so far depend on the energy fraction
    only.
  */
  /*!
    \fn Beam_Base::Beam_Base(std::string,const ATOOLS::Flavour,const double,const double,const int)
    \brief The constructor, initializes the full structure of this beam spectrum.

    In this constructor the incoming vectors are constructed from the direction flag, the energy of
    the beam in the lab frame and the mass of the beam particle. By default, also the outgonig momentum
    is set as the incoming one.
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
    important once Weizsaecker-Williams type spectra are considered.
  */
  /*!
    \fn double Beam_Base::Weight(ATOOLS::Flavour fl=Flavour(kf_none))
    Returns the weight. A flavour checking could be included in the future.
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
