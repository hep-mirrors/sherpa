#ifndef Beam_Spectra_Handler_h
#define Beam_Spectra_Handler_h

#include "Beam_Base.H"
#include "Vector.H"
#include "Poincare.H"
#include "Data_Read.H"

namespace BEAM {
  class Beam_Spectra_Handler {
  protected:
    Beam_Base          ** p_BeamBase;
    AMATOOLS::Poincare    m_CMSBoost;
    AMATOOLS::Poincare    m_LABBoost;

    bool                  m_asymmetric;
    int                   m_mode; 
    double                m_mass12, m_mass22, m_x1, m_x2;
    double                m_exponent[2];
    double                m_smin,m_smax,m_splimits[3],m_ylimits[2];
    AMATOOLS::Vec4D       m_fiXVECs[2];
    std::string           m_name,m_type;

  private:
    bool InitKinematics(AORGTOOLS::Data_Read *);
    bool SpecifySpectra(AORGTOOLS::Data_Read *);
    bool InitializeLaserBackscattering(AORGTOOLS::Data_Read *,int);
    bool InitializeMonochromatic(AORGTOOLS::Data_Read *,int);
  public:
    Beam_Spectra_Handler(AORGTOOLS::Data_Read *);
    Beam_Spectra_Handler(int *,APHYTOOLS::Flavour *,double *,APHYTOOLS::Flavour *,double *,double *);
    ~Beam_Spectra_Handler();
    
    void   Output();
    bool   CheckConsistency(APHYTOOLS::Flavour *,APHYTOOLS::Flavour *);
    bool   CheckConsistency(APHYTOOLS::Flavour *);
    bool   MakeBeams(AMATOOLS::Vec4D *,double,double);
    bool   CalculateWeight(double);
    double Weight(APHYTOOLS::Flavour *);
    double Weight();

    void   BoostInCMS(AMATOOLS::Vec4D *,int);
    void   BoostInLab(AMATOOLS::Vec4D *,int);

    void   SetName(std::string _name)     { m_name         = _name; }
    void   SetSprimeMin(double _spl)      { m_splimits[0]  = _spl; }
    void   SetSprimeMax(double _spl)      { m_splimits[1]  = _spl; }
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
    \var Beam_Base * Beam_Spectra_Handler::BeamBase[2]
    The two beam bases, one for each beam. They are initialized (through the method SpecifiySpectra)
    and owned by the Beam_Spectra_Handler which in turn has to delete them after the run. At the 
    moment two types are supported, namely 
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
    \var AMATOOLS::Poincare Beam_Spectra_Handler::m_m_CMSBoost
    A boost from the c.m. system of the incoming beams to the c.m. system of the
    outgoing bunches.
  */
  /*!
    \var AMATOOLS::Poincare Beam_Spectra_Handler::m_LABBoost
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
    \verbatim
    m_splimits[0,1] = \endverbatim
    \f$s'_{\rm min, max}\f$
    and 
    \verbatim
    m_splimits[2] = \endverbatim
    \f$s_{\rm beams}\f$.
  */ 
  /*!
    \var Beam_Spectra_Handler::m_ylimits[2]
    The rapidity region covered. It is per default set to the range \f$y \in [-10,10]\f$. In fact 
    this range should be calculated from the range of the BeamBases.
    \todo Rapidity range from BeamBases.
  */
  /*!
    \var AMATOOLS::Vec4D Beam_Spectra_Handler::m_fiXVECs[2]
    The c.m. momenta of the two incoming beams.
  */
  /*!
    \fn bool Beam_Spectra_Handler::Beam_Spectra_Handler(AORGTOOLS::Data_Read *)
    The explicit constructor managing the initialisation of the beams through calls to
    SpecifySpectra(AORGTOOLS::Data_Read *) and InitKinematics(AORGTOOLS::Data_Read *).

    Having succeeded in the initialization of the two Beam_Bases, m_mode is also
    determined.
  */
  /*!
    \fn bool Beam_Spectra_Handler::SpecifySpectra(AORGTOOLS::Data_Read *);
    SpecifySpectra reads in the specific spectra of the beams and the generator and initializes
    the spectra through the corresponidng Initialize methods. At the moment, two such methods
    are available, namely InitializeLaserBackscattering(AORGTOOLS::Data_Read *,int)
    and InitializeMonochromatic(AORGTOOLS::Data_Read *,int). For every new spectrum, such
    a method has to be provided as well and SpecifySpectra has to be modified accordingly.
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitKinematics(AORGTOOLS::Data_Read *);
    Here the two lab momenta of the BeamBases are combined to determnie the c.m. energy of
    the process in their c.m. frame. Furthermore, the two corresponding c.m. momenta of the
    incoming beams are constructed and it is decided whether additional boosts to the lab frame
    are mandatory. Also, the characteristic exponents are filled. 
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitializeLaserBackscattering(AORGTOOLS::Data_Read *,int)
    Characteristic parameters for laser back scattering are read in and the spectrum is
    instantiated through a call to 
    Laser_Backscattering::Laser_Backscattering(const APHYTOOLS::Flavour,const double,const double,
                                               const double,const double,const int,
					       const int,const int,bool &)
  */
  /*!
    \fn bool Beam_Spectra_Handler::InitializeMonochromatic(AORGTOOLS::Data_Read *,int);
    Characteristic parameters for monochromatic beams are read in and the spectrum is
    instantiated through a call to 
    Monochromatic::Monochromatic(const APHYTOOLS::Flavour,const double,const double,bool &)
  */
  /*!
    \fn bool Beam_Spectra_Handler::CheckConsistency(APHYTOOLS::Flavour *,APHYTOOLS::Flavour *)
    This checks whether the incoming and outgoing flavours match the flaovurs needed by the
    two BeamBases.
  */
  /*!
    \fn bool Beam_Spectra_Handler::MakeBeams(AMATOOLS::Vec4D *,double,double);
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
    \fn double Beam_Spectra_Handler::Weight(APHYTOOLS::Flavour * fl)
    The weight corresponding to CalculateWeight(double), basically the product of the two
    individual weights.
  */
}



 
#endif



