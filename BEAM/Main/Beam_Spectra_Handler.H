#ifndef BEAM_Main_Beam_Spectra_Handler_H
#define BEAM_Main_Beam_Spectra_Handler_H

#include "BEAM/Main/Beam_Parameters.H"
#include "BEAM/Main/Kinematics_Base.H"
#include "BEAM/Main/Weight_Base.H"
#include <array>

namespace BEAM {
  class Beam_Spectra_Handler {
  private:
    std::array<Beam_Base*, 2> p_BeamBase;
    Kinematics_Base  * p_kinematics;
    Weight_Base      * p_weight;
    Beam_Parameters    m_parameters;
    beammode           m_beammode;
    collidermode       m_collidermode;
    bool               m_on;
    int                m_mode,m_polarisation;
    std::string        m_name,m_type;

    bool InitTheBeams();
    bool InitTheKinematics();
    bool InitTheWeight();
  public:
    Beam_Spectra_Handler();
    ~Beam_Spectra_Handler();

    void FixPositions();
    inline bool                       On()           const { return m_on; }
    inline const beammode           & BeamMode()     const { return m_beammode; }
    inline const collidermode       & ColliderMode() const { return m_collidermode; }
    inline Beam_Base                * GetBeam(size_t i)    { return p_BeamBase[i]; }

    inline double Upper1() const { return p_BeamBase[0]->Xmax(); }
    inline double Upper2() const { return p_BeamBase[1]->Xmax(); }
    inline double Peak()   const { return p_BeamBase[0]->Peak()*p_BeamBase[1]->Peak(); }

    inline bool MakeBeams(ATOOLS::Vec4D_Vector& moms) { return (*p_kinematics)(moms); }
    inline void AssignKeys(ATOOLS::Integration_Info *const info) const {
      p_kinematics->AssignKeys(info);
      p_weight->AssignKeys(info);
    }
    inline void SetLimits() const { p_kinematics->SetLimits(); }
    inline void BoostInCMS(ATOOLS::Vec4D & p) const { p_kinematics->BoostInCMS(p); }
    inline void BoostInLab(ATOOLS::Vec4D & p) const { p_kinematics->BoostInLab(p); }
    inline void BoostBackLab(ATOOLS::Vec4D & p) { p_kinematics->BoostBackLab(p); }
    inline void BoostInCMS(ATOOLS::Vec4D_Vector& p) const {
      p_kinematics->BoostInCMS(p);
    }
    inline void BoostInLab(ATOOLS::Vec4D_Vector& p) const {
      p_kinematics->BoostInLab(p);
    }

    inline void SetType(const std::string & type) { m_type = type; }
    inline void SetSprimeMin(const double & smin) const { p_kinematics->SetSprimeMin(smin); }
    inline void SetSprimeMax(const double & smax) const { p_kinematics->SetSprimeMax(smax); }
    inline double SprimeMin() const { return p_kinematics->SprimeMin(); }
    inline double SprimeMax() const { return p_kinematics->SprimeMax(); }
    inline double Sprime()    const { return p_kinematics->Sprime(); }
    inline double Y()         const { return p_kinematics->Y(); }
    inline double Exponent(const size_t & i) const {
      return p_kinematics->Exponent(i);
    }
    inline bool CalculateWeight(const double & scale) const {
      return p_weight->Calculate(scale);
    }
    inline double Weight() const {
      return (*p_weight)();
    }

    // possible deprecated methods
    int    Polarisation() const { return m_polarisation; }
    bool   CheckConsistency(const ATOOLS::Flavour *,const ATOOLS::Flavour *);
    bool   CheckConsistency(const ATOOLS::Flavour *);
    void   BoostFixedTarget();
    double SqLam(double s,double s1,double s2);
  };

  /*!
    \namespace BEAM
    The namespace BEAM houses all classes that are employed to generate
    beam spectra. In the framework of both the SHERPA package and of the
    program AMEGIC the following nomenclature is assumed :
    - There are incoming beams at a certain energy, the nominal energy of the
      beam in the collider,
    - these beams then result in bunches of interacting particles which have
      an energy distribution, and, maybe, a \f$k_\perp\f$ distribution of
      transverse momenta w.r.t.the beam axis.
    - The bunch particles then can have a substructure, i.e. they might consist of
      partons, which then interact.

    As an illustrative example, consider the case of DIS of an electron on a photon.
    The incoming electron beam emits bunches of photons that might or might not
    resolved into partons during the interaction with the proton. In the BEAM
    namespace, the energy (and, maybe, virtuality) distribution of the photons
    is handled.
  */
  /*!
    \file
    \brief Contains the Beam_Spectra_Handler
  */
  /*!
    \class Beam_Spectra_Handler
    \brief Handler of all Beam related issues.
    This class manages the set-up of the incoming (bunch-)particles from the beams according to
    the strategy defined through the global Settings instance. Before coming
    into full effect during integration or event generation, this handler
    initalises a suitable Beam treatment (Beam_Bases) for both
    beams and uses them to generate corresponding weights, i.e. energy distributions. At the moment,
    all outgonig bunch particles are still collinear to the incoming beams, but this is going
    to change in the (near) future.
    \todo Allow for non-collinear bunches
  */
  /*!
    \var std::array<Beam_Base *, 2> Beam_Spectra_Handler::p_BeamBase
    Pointers to the two beam bases, one for each beam. They are initialized (through the method
    SpecifiySpectra) and owned by the Beam_Spectra_Handler which in turn has to delete them after
    the run. At the moment two types are supported, namely
    - monochromatic beams
    - photon beams stemming from electron through Laser backscattering.
    More types will be implemented soon, exmples include
    - Weizsaecker-Williams photons
    - electrons after Beamstrahlung
    - beams with a Gaussian energy-spread.
    \sa Beam_Base.
    \todo include more types of BeamBases
  */
  /*!
    \var int Beam_Spectra_Handler::m_mode
    Indicates the mode of the beam handling, or, better phrased, whether either or both of
    the beams are monochromatic:
    - 0 both beams are monochromatic, no beam handling needed.
    - 1 beam one is not monochromatic, beam two is.
    - 2 beam one is monochromatic, beam two is not.
    - 3 both beams are not monochromatic.
  */
  /*!
    var std::string Beam_Spectra_Handler::m_name
    Name of the Beam_Spectra_Handler.
  */
  /*!
    var std::string Beam_Spectra_Handler::m_type
    Type of the Beam_Spectra_Handler, it consists of the types of the BeamBases.
  */
  /*!
    \fn Beam_Spectra_Handler::Beam_Spectra_Handler()
    The explicit constructor managing the initialisation of the beams through calls to
    SpecifySpectra() and InitKinematics().
    Having succeeded in the initialization of the two Beam_Bases, m_mode is also
    determined. It is foreseen that at this point also the beam geometries - if
    necessary - will be fixed. Having read in the incoming energies also the c.m. energy
    of the collision is set in the Run_Paraemeters.

    \todo Enable beam geometries, especially for pile-up events etc.
  */
  /*!
    \fn bool Beam_Spectra_Handler::CheckConsistency(const ATOOLS::Flavour *,const ATOOLS::Flavour *)
    This checks whether the two sets of flavours match the flaovurs inherent to the
    two BeamBases. If this is the case, true is returned.
  */
  /*!
    \fn bool Beam_Spectra_Handler::CheckConsistency(const ATOOLS::Flavour *)
    This checks whether the flavours match the bunches of the two BeamBases. If this is the case,
    true is returned.
  */
  /*!
    \fn bool Beam_Spectra_Handler::MakeBeams(ATOOLS::Vec4D *,double,double);
    Depending on the \f$s'\f$-value handed over as arguments, two matching vectors for the
    outgoing bunches in their c.m. frame (out) are constructed. Then the energy fractions in the
    c.m. system (in) of the incoming beams are determined with help of the other argument, the
    rapidity \f$y\f$ according to
    \f[
    \hat E^{(in)}_{1,2} = \exp\left(\pm y\right)
    \f]
    and the boost linking the two frames, CMBoost is initialized. This boost is then used
    to bring the c.m. vectors into the correct frame, i.e. the c.m. frame
    of the beams, i.e.
     \f[
     p^{(out)}_{1,2} \Longrightarrow p^{(in)}_{1,2}\,.
     \f]
  */
  /*!
    \fn bool Beam_Spectra_Handler::CalculateWeight(double)
    This method calculates the two beam densities, i.e their weights, according to the spectra
    depending on the relative energy fractions and - eventually - on a scale, which is passed as an
    argument. The weight calculation proceeds via calls to the specific methods
    Beam_Base::CalculateWeight(double,double),
    where the energy fractions and the scale are passed as arguments.
  */
  /*!
    \fn double Beam_Spectra_Handler::Weight()
    The weight corresponding to CalculateWeight(double), basically the product of the two
    individual weights. Unlike in the ISR handling this does not depend on the
    particle flavours for the implemented spectra, i.e. EPA and Pomeron.
  */
}




#endif
