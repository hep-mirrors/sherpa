#ifndef ISR_Handler_H
#define ISR_Handler_H

#include "Vector.H"
#include "Poincare.H"
#include "ISR_Base.H"
#include "Info_Key.H"

namespace PDF {

  class ISR_Handler {
  protected:

    ISR_Base **p_isrbase;

    ATOOLS::Integration_Info *p_info;
    ATOOLS::Info_Key          m_spkey, m_ykey, m_xkey, m_zkey[2], m_kpkey[2];

    std::string m_name, m_type;

    int m_mode, m_kmrmode; 

    double m_mass2[2], m_exponent[2], m_x[2], m_mu2[2];
    double m_splimits[3],m_ylimits[2], m_kplimits[3], m_zlimits[2];
    double m_fixed_smin, m_fixed_smax, m_weight, m_flux;

    ATOOLS::Poincare m_cmsboost, m_kmrboost, m_kmrrot;
    ATOOLS::Vec4D    m_fixvecs[2], m_kp[2], p_cms[2];

    void Init(const double *splimits,const double *kplimits);

    const ATOOLS::Flavour  DiQuark(const ATOOLS::Flavour &,const ATOOLS::Flavour &);
    bool CheckRemnantKinematics(const ATOOLS::Flavour &,double,int);

  public:

    // constructor
    ISR_Handler(ISR_Base **isrbase,const double *splimits,const double *kplimits);

    // destructor
    ~ISR_Handler();

    // member functions
    bool CheckConsistency(ATOOLS::Flavour *bunches,ATOOLS::Flavour *partons);
    bool CheckConsistency(ATOOLS::Flavour *partons);
    void SetPartonMasses(ATOOLS::Flavour *fl); 

    bool MakeISR(ATOOLS::Vec4D *const p,const size_t n);

    bool CalculateWeight(const double scale);
    bool CalculateWeight2(const double scale);

    double Weight(ATOOLS::Flavour *flavour);
    double Weight2(ATOOLS::Flavour *flavour);

    void BoostInCMS(ATOOLS::Vec4D *p,const size_t n);
    void BoostInLab(ATOOLS::Vec4D *p,const size_t n);

    void SetSprimeMin(const double spmin);
    void SetSprimeMax(const double spmax);

    void SetFixedSprimeMin(const double spmin);
    void SetFixedSprimeMax(const double spmax);

    void AssignKeys(ATOOLS::Integration_Info *const info);
    void SetLimits();

    // inline functions 
    inline void SetPole(const double pole)      { m_splimits[2] = pole; }
    inline void SetYMin(const double ymin)      { m_ylimits[0]  = ymin; }
    inline void SetYMax(const double ymax)      { m_ylimits[1]  = ymax; }
    inline void SetName(const std::string name) { m_name = name;        }

    inline std::string Name() { return m_name; }
    inline std::string Type() { return m_type; }

    inline PDF_Base        *PDF(const size_t beam)  { return p_isrbase[beam]->PDF();     }
    inline ATOOLS::Flavour  Flav(const size_t beam) { return p_isrbase[beam]->Flavour(); }

    inline int On()    { return m_mode;    }
    inline int KMROn() { return m_kmrmode; }

    inline double *SprimeRange()         { return m_splimits;    }
    inline double *YRange()              { return m_ylimits;     }
    inline double *KMRZRange()           { return m_zlimits;     }
    inline double *KMRKpRange()          { return m_kplimits;    }
    inline double  Exponent(const int i) { return m_exponent[i]; }

    inline double SprimeMin() { return m_splimits[0]; }
    inline double SprimeMax() { return m_splimits[1]; }
    inline double Pole()      { return m_splimits[2]; }
    inline double YMin()      { return m_ylimits[0];  }
    inline double YMax()      { return m_ylimits[1];  }
    
    inline double Upper1() { return p_isrbase[0]->XMax(); }
    inline double Upper2() { return p_isrbase[1]->XMax(); }

    inline double X1() { return m_x[0]; }
    inline double X2() { return m_x[1]; }

    inline double Flux() { return m_flux; }
  
  };// end of class ISR_Handler

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
    \var ISR_Base ** ISR_Handler::p_isrbase
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
    and outgonig particles of the two p_isrbases. If this is the case, true is returned.
    This method is largely similar to the corresponding one in the BEAM::Beam_Spectra_Handler.
  */
  /*!
    \fn bool ISR_Handler::CheckConsistency(ATOOLS::Flavour *)
    This checks whether the flavours are allowed to be used as outgonig flavours in the
    two p_isrbases. If this is the case, true is returned. This method is largely similar to 
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

}// end of namespace PDF 

#endif



