#ifndef SHERPA_Initialization_Initialization_Handler_H
#define SHERPA_Initialization_Initialization_Handler_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Exception.H"

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"
#include "SHERPA/PerturbativePhysics/MI_Handler.H"
#include "SHERPA/PerturbativePhysics/Shower_Handler.H"
#include "SHERPA/PerturbativePhysics/Matrix_Element_Handler.H"
#include "SHERPA/Tools/Analysis_Interface.H"
#include "SHERPA/Tools/Output_Base.H"
#include "SHERPA/Tools/Definitions.H"
#include "YFS/Main/YFS_Handler.H"

#include <map>
#include <string>
#include <iostream>

namespace ATOOLS {
  class Variations;
  class Fragmentation_Base;
  class Settings;
}

namespace SHERPA {

  class Matrix_Element_Handler;
  class Hard_Decay_Handler;
  class Beam_Remnant_Handler;
  class Colour_Reconnection_Handler;
  class Decay_Handler_Base;
  class Soft_Photon_Handler;
  class Event_Reader_Base;
  class Analysis_Interface;
  class Filter;

  class Initialization_Handler: public ATOOLS::Terminator_Object {
  private :
    std::string m_evtform;
    std::set<std::string> m_pdflibs;
    std::map<PDF::isr::id, std::array<std::string, 2> >  m_defsets;
    std::map<PDF::isr::id, std::vector<size_t> > m_bunchtags;

    ATOOLS::Flavour m_bunch_particles[2];

    eventtype::code m_mode;
    bool            m_savestatus;

    MODEL::Model_Base            * p_model;
    BEAM::Beam_Spectra_Handler   * p_beamspectra;
    PDF::ISR_Handler_Map           m_isrhandlers;
    REMNANTS::Remnant_Handler_Map  m_remnanthandlers;
    Matrix_Element_Handler       * p_mehandler;
    Hard_Decay_Handler           * p_harddecays;
    Shower_Handler_Map             m_showerhandlers;
    Beam_Remnant_Handler         * p_beamremnants;
    Colour_Reconnection_Handler  * p_reconnections;
    ATOOLS::Fragmentation_Base   * p_fragmentation;
    Soft_Collision_Handler_Map     m_schandlers;
    Decay_Handler_Base           * p_hdhandler;
    MI_Handler_Map                 m_mihandlers;
    Soft_Photon_Handler          * p_softphotons;
    Event_Reader_Base            * p_evtreader;
    Analysis_Vector                m_analyses;
    Output_Vector                  m_outputs;
    ATOOLS::Variations           * p_variations;
    Filter                       * p_filter;
    YFS::YFS_Handler            * p_yfshandler;
    
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
    bool InitializeTheYFS();

    void RegisterDefaults();
    bool CheckBeamISRConsistency();
    void LoadPDFLibraries(ATOOLS::Settings& settings);
    void InitISRHandler(const PDF::isr::id & pid,ATOOLS::Settings& settings);
    void DefineBunchFlavours(ATOOLS::Settings& settings);
    void SetGlobalVariables();
    void PrepareTerminate();

  public :

    Initialization_Handler();
    ~Initialization_Handler();

    bool InitializeTheFramework(int =-1);
    bool CalculateTheHardProcesses();
    inline const bool AddBeamRescatters() const {
      return m_isrhandlers.find(PDF::isr::bunch_rescatter)!=m_isrhandlers.end();
    }
    Matrix_Element_Handler           * GetMatrixElementHandler()      const { return p_mehandler; }
    Decay_Handler_Base               * GetHDHandler()                       { return p_hdhandler; }
    Analysis_Vector                  * GetAnalyses()                        { return &m_analyses; }
    Output_Vector                    * GetOutputs()                         { return &m_outputs; }
    MODEL::Model_Base                * GetModel()                     const { return p_model; }
    PDF::ISR_Handler_Map             * GetISRHandlers()                     { return &m_isrhandlers; }
    YFS::YFS_Handler                 * GetYFSHandler()                      { return p_yfshandler; }
    BEAM::Beam_Spectra_Handler       * GetBeamSpectraHandler()        const { return p_beamspectra; }
    Hard_Decay_Handler               * GetHardDecayHandler()          const { return p_harddecays; }
    const Shower_Handler_Map         & GetShowerHandlers()            const { return m_showerhandlers; }
    Beam_Remnant_Handler             * GetBeamRemnantHandler()        const { return p_beamremnants; }
    Colour_Reconnection_Handler      * GetColourReconnectionHandler() const { return p_reconnections; }
    ATOOLS::Fragmentation_Base       * GetFragmentation()             const { return p_fragmentation; }
    Soft_Collision_Handler_Map       * GetSoftCollisionHandlers()           { return &m_schandlers; }
    ATOOLS::Flavour                  * GetBunches()                         { return m_bunch_particles; }
    MI_Handler_Map                   * GetMIHandlers()                      { return &m_mihandlers; }
    REMNANTS::Remnant_Handler_Map    & GetRemnantHandlers()                 { return m_remnanthandlers; }
    Soft_Photon_Handler              * GetSoftPhotonHandler()         const { return p_softphotons; }
    Event_Reader_Base                * GetEventReader()               const { return p_evtreader; }
    ATOOLS::Variations               * GetVariations()                const { return p_variations; }
    Filter                           * GetFilter()                    const { return p_filter; }
    inline PDF::ISR_Handler       * GetISRHandler(const PDF::isr::id id)  { return m_isrhandlers[id]; }
    inline Soft_Collision_Handler * GetSCHandler(const PDF::isr::id id)   { return m_schandlers[id]; }
    inline MI_Handler             * GetMIHandler(const PDF::isr::id id)   { return m_mihandlers[id]; }

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
    - the fragmentation
    
    These modules are later on taken by the main steering class, Sherpa, and put into
    corresponding Event_Phase_Handlers. After termination of the run, they are all deleted.
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
    \var Initialization_Handler::m_isrhandlers
    The ISR handling (the pdf's) that is initialized by the Initialization_Handler.
  */
  /*!
    \var Matrix_Element_Handler * Initialization_Handler::p_mehandler
    A pointer to the abstract wrapper around the matrix elements that are initialized
    by the Initialization_Handler.
  */
  /*!
    \var Initialization_Handler::m_showerhandlers
    The abstract wrappers around the shower module. It is initialized by the
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
    \var Initialization_Handler::p_hdhandler
    A pointer to the wrapper around the hadron decay code that will be used during the run.
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
    */


}


#endif
