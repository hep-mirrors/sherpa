#ifndef SHERPA_Initialization_Initialization_Handler_H
#define SHERPA_Initialization_Initialization_Handler_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Exception.H"

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Tools/Analysis_Interface.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/Definitions.H"

#include <map>
#include <string>
#include <iostream> 

namespace ATOOLS {
  class Variations;
}

namespace SHERPA {

  class Matrix_Element_Handler;
  class Hard_Decay_Handler;
  class Beam_Remnant_Handler;
  class Colour_Reconnection_Handler;
  class Fragmentation_Handler;
  class Decay_Handler_Base;
  class MI_Handler;
  class Soft_Photon_Handler;
  class Lund_Interface;
  class Event_Reader_Base;
  class Analysis_Interface;
  class Soft_Collision_Handler;
  class Filter;
  
  class Initialization_Handler: public ATOOLS::Terminator_Object {
  private :
    std::string m_evtform;
    std::set<std::string> m_pdflibs;

    double          m_bunch_splimits[2];
    ATOOLS::Flavour m_bunch_particles[2];

    eventtype::code m_mode;
    bool            m_savestatus;

    MODEL::Model_Base           * p_model;
    BEAM::Beam_Spectra_Handler  * p_beamspectra;
    PDF::ISR_Handler_Map          m_isrhandlers;
    REMNANTS::Remnant_Handler   * p_remnants;
    Matrix_Element_Handler      * p_mehandler;
    Hard_Decay_Handler          * p_harddecays;
    Shower_Handler_Map            m_showerhandlers;
    Beam_Remnant_Handler        * p_beamremnants;
    Colour_Reconnection_Handler * p_reconnections;
    Fragmentation_Handler       * p_fragmentation;
    Soft_Collision_Handler      * p_softcollisions;
    Decay_Handler_Base          * p_hdhandler;
    MI_Handler                  * p_mihandler; 
    Soft_Photon_Handler         * p_softphotons;
    Event_Reader_Base           * p_evtreader;
    Analysis_Vector               m_analyses;
    Output_Vector                 m_outputs;
    ATOOLS::Variations          * p_variations;
    Filter                      * p_filter;
    
    void ShowParameterSyntax();
    void CheckVersion();
    bool CompareVersions(const size_t& a1,const size_t& b1,const size_t& c1,
                         const size_t& a2,const size_t& b2,const size_t& c2);
    void LoadLibraries();
  private :
    bool InitializeTheIO();
    bool InitializeTheModel();
    bool InitializeTheBeams();
    bool InitializeThePDFs();
    bool InitializeTheRemnants();
    bool InitializeTheMatrixElements();
    bool InitializeTheUnderlyingEvents();
    bool InitializeTheHardDecays();
    bool InitializeTheShowers();
    bool InitializeTheBeamRemnants();
    bool InitializeTheColourReconnections();
    bool InitializeTheFragmentation();
    bool InitializeTheSoftCollisions();
    bool InitializeTheHadronDecays();
    bool InitializeTheSoftPhotons();
    bool InitializeTheFilter();
    bool InitializeTheAnalyses();
    bool InitializeTheReweighting(ATOOLS::Variations_Mode);

    void RegisterDefaults();
    bool CheckBeamISRConsistency();
    void SetGlobalVariables();
    void PrepareTerminate();

  public :

    Initialization_Handler();
    ~Initialization_Handler();

    bool InitializeTheFramework(int =-1);
    bool CalculateTheHardProcesses();

    Matrix_Element_Handler       * GetMatrixElementHandler()  const { return p_mehandler; }
    Decay_Handler_Base           * GetHDHandler() { return p_hdhandler; }
    inline PDF::ISR_Handler      * GetISRHandler(const PDF::isr::id id) { return m_isrhandlers[id]; }
    Analysis_Vector              * GetAnalyses()                    { return &m_analyses;       }
    Output_Vector                * GetOutputs()                     { return &m_outputs;        }
    MODEL::Model_Base            * GetModel()                 const { return p_model;           }
    PDF::ISR_Handler_Map         * GetISRHandlers()                 { return &m_isrhandlers;    }
    BEAM::Beam_Spectra_Handler   * GetBeamSpectraHandler()    const { return p_beamspectra;     }
    Hard_Decay_Handler           * GetHardDecayHandler()      const { return p_harddecays;      }
    const Shower_Handler_Map     & GetShowerHandlers()        const { return m_showerhandlers;  }
    Beam_Remnant_Handler         * GetBeamRemnantHandler()    const { return p_beamremnants;    }
    Colour_Reconnection_Handler  * GetColourReconnectionHandler() const { return p_reconnections;  }
    Fragmentation_Handler        * GetFragmentationHandler()  const { return p_fragmentation;   }
    Soft_Collision_Handler       * GetSoftCollisionHandler()  const { return p_softcollisions;  }
    ATOOLS::Flavour              * GetBunches()                     { return m_bunch_particles; }
    MI_Handler                   * GetMIHandler()             const { return p_mihandler;       }
    Soft_Photon_Handler          * GetSoftPhotonHandler()     const { return p_softphotons;     }
    Event_Reader_Base            * GetEventReader()           const { return p_evtreader;       }
    ATOOLS::Variations           * GetVariations()            const { return p_variations;      }
    Filter                       * GetFilter()                const { return p_filter;          }
    
    eventtype::code Mode() const { return m_mode; };

  };

  /*!
    \file 
    \brief Contains the class SHERPA::Initialization_Handler
  */                                  

  /*!
    \class Initialization_Handler
    \brief It is in charge of initialising and deleting all physics modules used during a run.
    
    This class takes care of the proper initialization of all physics modules to be initialized
    during the set up of the event generator. In particular, these contain
    - the model which contains all physics parameters (like coupling constants, 
      mixing matrices, etc.)
    - beam spectra (which is particularly important for physics processes involving for instance
      laser back scattering of photons off electrons)
    - the pdf's through the ISR_Handler,
    - the matrix elements, if needed they are also evaluated to yield total cross sections,
    - the parton shower in the initial and final state,
    - the handling of the beam remnants,
    - the fragmentation which at the moment is dealt with through an interface to the Lund string
      of Pythia and the subsequent hadron decays by the same program.

    These modules are later on taken by the main steering class, Sherpa, and put into 
    corresponding Event_Phase_Handlers. After termination of the run, they are all deleted.
  */
  /*!
    \var double Initialization_Handler::m_bunch_splimits[2]
    The boundaries for the s' after ISR. These paraermters are read in in the Initialization_Handler
    and passed to the ISR_Handler. They are given in units of the total energy squared under
    consideration.
  */
  /*!
    \var ATOOLS::Flavour Initialization_Handler::m_bunch_particles[2]
    The bunch particles, i.e. the particles that are resolved into their constituents by the pdf's.
  */
  /*!
    \var MODEL::Model_Base * Initialization_Handler::p_model
    A pointer to the model that is initialized by the Initialization_Handler.
  */
  /*!
    \var BEAM::Beam_Spectra_Handler * Initialization_Handler::p_beamspectra
    A pointer to the beam spectra that are initialized by the Initialization_Handler.
  */
  /*!
    \var PDF::ISR_Handler * Initialization_Handler::p_isrhandler
    A pointer to the ISR handling (the pdf's) that is initialized by the Initialization_Handler.
  */
  /*!
    \var Matrix_Element_Handler * Initialization_Handler::p_mehandler
    A pointer to the abstract wrapper around the matrix elements that are initialized
    by the Initialization_Handler.
  */
  /*!
    \var Shower_Handler * Initialization_Handler::p_showerhandler
    A pointer to the abstract wrapper around the shower module. It is initialized by the 
    Initialization_Handler. 
  */
  /*!
    \var Beam_Remnant_Handler * Initialization_Handler::p_beamremnants
    A pointer to the beam remnant handling used thoughout the run. 
  */
  /*!
    \var Fragmentation_Handler * Initialization_Handler::p_fragmentation
    A pointer to the wrapper around the fragmentation code that will be used during the run.
  */
  /*!
    \var Hadron_Decay_Handler * Initialization_Handler::p_hadrondecays
    A pointer to the wrapper around the hadron decay code that will be used during the run.
  */
  /*!
    \var HDHandlersMap Initialization_Handler::m_hdhandlers
    This is map with pointers to the corresponding hadron decay handlers
    <TABLE>
    <TR> <TD>1st column</TD> <TD>2nd column</TD> </TR>
    <TR> <TD>Lund</TD> <TD>pointer on Pythia</TD> </TR>
    <TR> <TD>Hadrons</TD> <TD>pointer on HADRONS</TD> </TR>
    </TABLE>
   */
  /*!
    \fn Initialization_Handler::~Initialization_Handler()
    The destructor deletes all physdics modules that have been used during the run.
  */
  /*!
    \fn bool Initialization_Handler::InitializeTheFramework()
    This method first reads in all particles, before subsequently the model, the beams and the pdf's
    are set up. After this has been done, consistency of the bunch particles and the beam 
    handling is checked for. If successul, the matrix elements, showers, beam remanants, 
    fragmentation and hadron decays are instantiated. If everything worked out, true is returned
    and some first corss section evaluation will start.
  */
  /*!
    \fn bool Initialization_Handler::CalculateTheHardProcesses()
    Returns true if the calculation of total cross sections by the matrix elements
    of choice was successful, otherwise, flase will be returned.
  */

  /*!
    \fn Initialization_Handler::InitializeTheHadronDecays()
    This method reads the fragmentation settings, especially the desired decay model.
     
    If the model is \c Hadrons it calls the HADRONS constructor for initialisation
    and the Hadron_Decay_Handler::EraseTreated method
    to remove the particles that HADRONS can cope with from the unstable hadron list. If then
    there are still some particles in the list, i.e. those that HADRONS can't deal with, it needs to
    initialise Pythia as well.

    If, on the other hand, the model is \c Lund, it calls the Pythia for initialisation.
    */
  

}


#endif
