#ifndef ISR_Handler_h
#define ISR_Handler_h

#include "Vector.H"
#include "Poincare.H"
#include "ISR_Base.H"

namespace PDF {
  class ISR_Handler {
  protected:
    ISR_Base           ** p_ISRBase;
    ATOOLS::Poincare    m_CMSBoost;
    int                   m_mode; 
    double                m_mass12, m_mass22, m_x1, m_x2;
    double                m_exponent[2];
    double                m_splimits[3],m_ylimits[2];
    ATOOLS::Vec4D       m_fiXVECs[2];
    std::string           m_name, m_type;
  public:
    ISR_Handler(int *,ATOOLS::Flavour *,ATOOLS::Flavour *,double *);
    ISR_Handler(ISR_Base **,double *);
    ~ISR_Handler();
    void   Init(double *);
    bool   CheckConsistency(ATOOLS::Flavour *,ATOOLS::Flavour *);
    bool   CheckConsistency(ATOOLS::Flavour *);
    void   SetPartonMasses(ATOOLS::Flavour * _fl); 
    bool   MakeISR(ATOOLS::Vec4D *,double,double);
    bool   CalculateWeight(double);
    bool   CalculateWeight2(double);
    double Weight(ATOOLS::Flavour *);
    double Weight2(ATOOLS::Flavour *);

    void   BoostInCMS(ATOOLS::Vec4D *,int);
    void   BoostInLab(ATOOLS::Vec4D *,int);

    void   SetName(std::string _name) { m_name = _name; }
    void   SetSprimeMin(double _spl)  { m_splimits[0]  = _spl; }
    void   SetSprimeMax(double _spl)  { m_splimits[1]  = _spl; }
    void   SetPole(double _spl)       { m_splimits[2] = _spl; }
    void   SetYMin(double _yl)        { m_ylimits[0]  = _yl; }
    void   SetYMax(double _yl)        { m_ylimits[1]  = _yl; }

    std::string          Name()            { return m_name; }
    std::string          Type()            { return m_type; }
    double               Exponent(int i)   { return m_exponent[i]; }
    double             * SprimeRange()     { return m_splimits; }
    double             * YRange()          { return m_ylimits; }
    double               SprimeMin()       { return m_splimits[0]; }
    double               SprimeMax()       { return m_splimits[1]; }
    double               Pole()            { return m_splimits[2]; }
    double               YMin()            { return m_ylimits[0]; }
    double               YMax()            { return m_ylimits[1]; }
    
    // Don't use this !!! It is for the ISRChannels only !!!
    double               Upper1()          { return p_ISRBase[0]->Xmax(); }
    double               Upper2()          { return p_ISRBase[1]->Xmax(); }
    double               X1()              { return m_x1; }
    double               X2()              { return m_x2; }
  
    PDF_Base           * PDF(int beam)     { return p_ISRBase[beam]->PDF(); }
    ATOOLS::Flavour      Flav(int beam)    { return p_ISRBase[beam]->Flav(); }
    int                  On()              { return m_mode; }
  };

  /*!
    \namespace PDF
    The namespace PDF houses all classes that are employed to generate
    parton spectra. In the framework of both the SHERPA package and of the
    program AMEGIC the following nomenclature is assumed :
    - There are incoming beams at a certain energy, the nominal energy of the 
      beam in the collider, which then result in bunches of interacting 
      particles which have an energy distribution, and, maybe, a \f$k_\perp\f$ 
      distribution of transverse momenta w.r.t.the beam axis.
    - The bunch particles then can have a substructure, i.e. they might consist of
      partons, which then interact in a hard subprocess.

    As an illustrative example, consider the case of DIS of an electron on a photon.
    The incoming electron beam emits bunches of photons that might or might not
    resolved into partons during the interaction with the proton. In the PDF namespace,
    the parton distribution funcitons of both the photon and the proton are handled.
  */
  /*!
    \class ISR_Handler
    \brief Manager of all Initial State Radiation that can be identifeid as parton
           distributions.
    This class manages all initial state radiation (ISR) according to the parton 
    distribution functions (PDFs) that are handed over. The ISR_Handler is initialized from 
    the SHERPA package or from Amegic. Before coming into full effect during integration
    or event generation, it initalises suitable ISR treatment through ISR_Bases that will
    contain the PDFs for each of the bunches.
  */
  /*!
    \var ISR_Base ** ISR_Handler::p_ISRBase
    Pointers to the two ISR bases, one for each bunch.

    \sa ISR_Base
  */
  /*!
     \var int ISR_Handler::m_mode
     The m_mode flag indicates what kind of ISR treatment is to be considered:
     - 0 no ISR for both bunchs 
     - 1 only bunch 1 experiences ISR
     - 2 only bunch 2 experiences ISR
     - 3 both bunches experience ISR.     
  */  
  /*!
    \var ATOOLS::Poincare ISR_Handler::m_CMSBoost
    A boost from the c.m. system of the incoming bunches to the c.m. system of the
    outgoing partons, which form the initial state of the hard interaction.
  */
  /*!
    \var double ISR_Handler::m_exponent[2]
    Characteristic exponents used for the integration.
  */
  /*!
    \var double ISR_Handler::m_splimits[3]
    \f$s'\f$-limits and characteristics:
    m_splimits[0,1] = \f$s'_{\rm min, max}\f$
    and 
    m_splimits[2] = \f$s_{\rm beams}\f$.
  */ 
  /*!
    \var ISR_Handler::m_ylimits[2]
    The rapidity region covered. It is per default set to the range \f$y \in [-10,10]\f$. In fact 
    this range should be calculated from the range of the BeamBases.
    \todo Rapidity range from BeamBases.
  */
  /*!
    \var double ISR_Handler::m_mass12
    Square of the mass of the first incoming particle.
  */
  /*!
    \var double ISR_Handler::m_mass22
    Square of the mass of the second incoming particle.
  */
  /*!
    \var double ISR_Handler::m_x1
    The energy fractions \f$x_{1}\f$ the outgoing bunch 1 have w.r.t. the corresponding 
    incoming beams.
  */
  /*!
    \var double ISR_Handler::m_x2
    The energy fractions \f$x_{2}\f$ the outgoing bunch 2 have w.r.t. the corresponding 
    incoming beams.
  */
  /*!
    var std::string ISR_Handler::m_name
    Name of the ISR_Handler.
  */
  /*!
    var std::string ISR_Handler::m_type
    Type of the ISR_Handler, it consists of the types of the BeamBases.
  */
  /*!
    \var ATOOLS::Vec4D ISR_Handler::m_fiXVECs[2]
    The c.m. momenta of the two incoming bunches.
  */
  /*!
    \fn bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *,ATOOLS::Flavour *)
    This checks whether the two sets of flavours match the flavours of the incoming
    and outgonig particles of the two p_ISRBases. If this is the case, true is returned.
    This method is largely similar to the corresponding one in the BEAM::Beam_Spectra_Handler.
  */
  /*!
    \fn bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *)
    This checks whether the flavours are allowed to be used as outgonig flavours in the
    two p_ISRBases. If this is the case, true is returned. This method is largely similar to 
    the corresponding one in the BEAM::Beam_Spectra_Handler.
  */  
  /*!
    \fn void ISR_Handler::SetPartonMasses(ATOOLS::Flavour * _fl)
    This sets the masses squared and the vectors such that they fit to the masses of the
    flavours.
  */
  /*!
    \fn bool ISR_Handler::MakeISR(ATOOLS::Vec4D *,double,double);
    Depending on the \f$s'\f$-value handed over as arguments, two matching vectors for the
    outgoing partons in their c.m. frame (out) are constructed. Then the energy fractions in the
    c.m. system (in) of the incoming bunches are determined with help of the other argument, the
    rapidity \f$y\f$ according to
    \f[
    \hat E^{(in)}_{1,2} = \exp\left(\pm y\right) 
    \f]
    and the boost linking the two frames, CMBoost is initialized. This boost is then used
    to bring the c.m. vectors into the correct frame, i.e. the c.m. frame 
    of the bunches, i.e.
    \f[
    p^{(out)}_{1,2} \Longrightarrow p^{(in)}_{1,2}\,.
    \f]
    This method is largely similar to the corresponding one in the
    BEAM::Beam_Spectra_Handler.
  */
  /*!
    \fn bool ISR_Handler::CalculateWeight(double)
    This method calculates the parton densities according to the two pdf's depending
    on the external x-values and the scale, which is passed as an argument. It should
    be noted that when calling the calculation of pdf weights, usually all weights for
    all partons inside a bunch particle are evaluated. In general this translates
    into eleven weights for a proton, one for each of the five lightest quark
    flavours, five for their antiflavours, and one for the gluon.
  */
  /*!
    \fn bool ISR_Handler::CalculateWeight2(double)
    This method calculates the parton densities according to the two pdf's depending
    on the external x-values and the scale, which is passed as an argument. The main point
    here is that the x-values are interchanged to allow for parton 1 coming from bunch 2
    and vice versa.
  */  
  /*!
    \fn double ISR_Handler::Weight(ATOOLS::Flavour * fl)
    This is the product of the actual values of the pdf's for the specific initial state 
    partons. It corresponds to the "untwisted" weight.
  */
  /*!
    \fn double ISR_Handler::Weight2(ATOOLS::Flavour * fl)
    This is the product of the actual values of the pdf's for the specific initial state 
    partons. It corresponds to the "twisted" weight, i.e. for parton 1 stemming from
    bunch 2 and vice versa.
  */
} 
#endif



