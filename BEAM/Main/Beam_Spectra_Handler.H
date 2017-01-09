#ifndef BEAM_Main_Beam_Spectra_Handler_H
#define BEAM_Main_Beam_Spectra_Handler_H

#include "BEAM/Main/Beam_Base.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Info_Key.H"

namespace ATOOLS { class Default_Reader; }

namespace BEAM {
  class Beam_Spectra_Handler {
  protected:
    Beam_Base          ** p_BeamBase;
    ATOOLS::Poincare      m_CMSBoost;
    ATOOLS::Poincare      m_LABBoost;

    bool                  m_asymmetric;
    int                   m_mode, m_polarisation; 
    double                m_mass12, m_mass22, m_x1, m_x2;
    double                m_exponent[2];
    double                m_splimits[3],m_ylimits[2];
    ATOOLS::Vec4D         m_fiXVECs[2];
    std::string           m_name,m_type;
    ATOOLS::Info_Key      m_spkey, m_ykey, m_xkey;
  private:
    bool InitKinematics(ATOOLS::Default_Reader *);
    bool SpecifySpectra(ATOOLS::Default_Reader *);
    bool InitializeLaserBackscattering(ATOOLS::Default_Reader *,int);
    bool InitializeSpectrumReader(ATOOLS::Default_Reader *,int);
    bool InitializeMonochromatic(ATOOLS::Default_Reader *,int);
    bool InitializeEPA(ATOOLS::Default_Reader *,int);
    void InitializeFlav(kf_code);
  public:
    Beam_Spectra_Handler(ATOOLS::Default_Reader *);
    ~Beam_Spectra_Handler();
    
    bool   Init();
    void   Output();
    bool   CheckConsistency(ATOOLS::Flavour *,ATOOLS::Flavour *);
    bool   CheckConsistency(ATOOLS::Flavour *);
    bool   MakeBeams(ATOOLS::Vec4D *);
    bool   CalculateWeight(double);
    double Weight(ATOOLS::Flavour * flin = NULL);
    double Weight(int * pol_types, double * dofs);

    void   BoostInCMS(ATOOLS::Vec4D *,int);
    void   BoostInLab(ATOOLS::Vec4D *,int);

    void   AssignKeys(ATOOLS::Integration_Info *const info);
    void   SetLimits();

    void   SetName(std::string _name)     { m_name         = _name; }
    void   SetSprimeMin(double _spl);
    void   SetSprimeMax(double _spl);
    void   SetYMin(double _yl)            { m_ylimits[0]   = _yl; }
    void   SetYMax(double _yl)            { m_ylimits[1]   = _yl; }
    std::string     Name()                { return m_name; }
    std::string     Type()                { return m_type; }
    double        * SprimeRange()         { return m_splimits; }
    double        * YRange()              { return m_ylimits; }
    double          SprimeMin()           { return m_splimits[0]; }
    double          SprimeMax()           { return m_splimits[1]; }
    double          YMin()                { return m_ylimits[0]; }
    double          YMax()                { return m_ylimits[1]; }
    double          Upper1()              { return p_BeamBase[0]->Xmax(); }
    double          Upper2()              { return p_BeamBase[1]->Xmax(); }
    double          Peak()                { return p_BeamBase[0]->Peak() *
					           p_BeamBase[1]->Peak(); }
    double          Exponent(const int i) { return m_exponent[i]; }
    int             On()                  { return m_mode; }
    int             Polarisation()        { return m_polarisation; }
    Beam_Base *     GetBeam(const int i)  { return p_BeamBase[i]; }
    
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
    the strategy named in the file that is inherently passed to the Beam_Spectra_Handler through 
    the pointer to the Data_Read object. Before coming into full efect during integration or 
    event generation, this handler initalises a suitable Beam treatment (Beam_Bases) for both 
    beams and uses them to generate corresponding weights, i.e. energy distributions. At the moment, 
    all outgonig bunch particles are still collinear to the incoming beams, but this is going 
    to change in the (near) future.
    \todo Allow for non-collinear bunches
  */  
  /*!
    \var Beam_Base ** Beam_Spectra_Handler::p_BeamBase
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
    \var ATOOLS::Poincare Beam_Spectra_Handler::m_CMSBoost
    A boost from the c.m. system of the incoming beams to the c.m. system of the
    outgoing bunches.
  */
  /*!
    \var ATOOLS::Poincare Beam_Spectra_Handler::m_LABBoost
    A boost from the lab system of the incoming beams to their c.m. system.
  */
  /*!
    \var bool Beam_Spectra_Handler::m_asymmetric
    This flag indicates wether the incoming beams are in an asymmetric situation, i.e.
    if the boost from the c.m. system of the incoming beams to their lab system is 
    mandatory.
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
    \var double Beam_Spectra_Handler::m_exponent[2]
    Characteristic exponents used for the integration. This is mainly a service for better
    performance of the phase space integration. At the moment,
    \verbatim
    m_exponent[0] = 0.5
    \endverbatim
    and
    \verbatim
    m_exponent[1] = 0.98 * ( p_BeamBase[0]->Exponent() + p_BeamBase[1]->Exponent()) 
    \endverbatim
    is given by characteristic exponents of the BeamBases.
  */
  /*!
    \var double Beam_Spectra_Handler::m_splimits[3]
    \f$s'\f$-limits and characteristics:
    m_splimits[0,1] = \f$s'_{\rm min, max}\f$
    and 
    m_splimits[2] = \f$s_{\rm beams}\f$.
  */ 
  /*!
    \var Beam_Spectra_Handler::m_ylimits[2]
    The rapidity region covered. It is per default set to the range \f$y \in [-10,10]\f$. In fact 
    this range should be calculated from the range of the BeamBases.
    \todo Rapidity range from BeamBases.
  */
  /*!
    \var double Beam_Spectra_Handler::m_mass12
    Square of the mass of the first incoming particle.
  */
  /*!
    \var double Beam_Spectra_Handler::m_mass22
    Square of the mass of the second incoming particle.
  */
  /*!
    \var double Beam_Spectra_Handler::m_x1
    The energy fractions \f$x_{1}\f$ the outgoing bunch 1 have w.r.t. the corresponding 
    incoming beams.
  */
  /*!
    \var double Beam_Spectra_Handler::m_x2
    The energy fractions \f$x_{2}\f$ the outgoing bunch 2 have w.r.t. the corresponding 
    incoming beams.
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
    \var ATOOLS::Vec4D Beam_Spectra_Handler::m_fiXVECs[2]
    The c.m. momenta of the two incoming beams.
  */
  /*!
    \fn Beam_Spectra_Handler::Beam_Spectra_Handler(ATOOLS::Data_Read *)
    The explicit constructor managing the initialisation of the beams through calls to
    SpecifySpectra(ATOOLS::Data_Read *) and InitKinematics(ATOOLS::Data_Read *).
    Having succeeded in the initialization of the two Beam_Bases, m_mode is also
    determined. It is foreseen that at this point also the beam geometries - if
    necessary - will be fixed. Having read in the incoming energies also the c.m. energy
    of the collision is set in the Run_Paraemeters.

    \todo Enable beam geometries, especially for pile-up events etc.
  */
  /*!
    \fn bool Beam_Spectra_Handler::SpecifySpectra(ATOOLS::Data_Read *);
    SpecifySpectra reads in the specific spectra of the beams and the corresponding generator 
    and initializes the spectra through appropriate InitializeXXX methods. At the moment, two 
    such methods are available, namely InitializeLaserBackscattering(ATOOLS::Data_Read *,int)
    and InitializeMonochromatic(ATOOLS::Data_Read *,int). For every new spectrum, such
    a method has to be provided as well and SpecifySpectra has to be modified accordingly.
    The idea is that in the long run, beam generators might be linked from outside by
    experimenters aiming at higher precision/better fit to beam data.
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitializeLaserBackscattering(ATOOLS::Data_Read *,int)
    Characteristic parameters for laser back scattering are read in and the spectrum is
    instantiated through a call to 
    Laser_Backscattering::Laser_Backscattering(const ATOOLS::Flavour,const double,const double,
                                               const double,const double,const int,
					       const int,const int,bool &)
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitializeMonochromatic(ATOOLS::Data_Read *,int);
    Characteristic parameters for monochromatic beams are read in and the spectrum is
    instantiated through a call to 
    Monochromatic::Monochromatic(const ATOOLS::Flavour,const double,const double,bool &)
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitKinematics(ATOOLS::Data_Read *);
    Here the two lab momenta of the BeamBases are combined to determnie the c.m. energy of
    the process in their c.m. frame. Furthermore, the two corresponding c.m. momenta of the
    incoming beams are constructed and it is decided whether additional boosts to the lab frame
    are mandatory. Also, the characteristic exponents are filled. 
  */
  /*!
    \fn bool Beam_Spectra_Handler::CheckConsistency(ATOOLS::Flavour *,ATOOLS::Flavour *)
    This checks whether the two sets of flavours match the flaovurs inherent to the
    two BeamBases. If this is the case, true is returned.
  */
  /*!
    \fn bool Beam_Spectra_Handler::CheckConsistency(ATOOLS::Flavour *)
    This checks whether the flavours match the bunches of the two BeamBases. If this is the case, 
    true is returned.
  */
  /*!
    \fn void Beam_Spectra_Handler::Output()
    Some control output of the type and name of this Handler and the setup of the BeamBases.
  */
  /*!
    \fn void Beam_Spectra_Handler::BoostInCMS(ATOOLS::Vec4D *,int)
    Boosts the n vectors into the CMS frame of the beam particles.
  */
  /*!
    \fn void Beam_Spectra_Handler::BoostInLab(ATOOLS::Vec4D *,int)
    Boosts the n vectors into the lab frame.
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
    \fn double Beam_Spectra_Handler::Weight(ATOOLS::Flavour * flin = NULL)
    The weight corresponding to CalculateWeight(double), basically the product of the two
    individual weights. This clearly depends on the flavours of the bunch particles here
    and might therefore allow to have, e.g., a spectrum of the outgoing electrons in
    laser backscattering. Such an option, however, has not been implemented yet. 
  */
}



 
#endif



