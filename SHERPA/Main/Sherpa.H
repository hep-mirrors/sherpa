#ifndef SHERPA_Main_Sherpa_H
#define SHERPA_Main_Sherpa_H

#include "ATOOLS/Org/Exception.H"
#include "SHERPA/Main/Filter.H"

namespace ATOOLS {
  class Blob_List;
  class Cluster_Amplitude;
}

namespace HepMC{
  class GenEvent;
}

namespace SHERPA {
  class Initialization_Handler;
  class Event_Handler;
  class Input_Output_Handler;
  class HepMC2_Interface;

  class Sherpa: public ATOOLS::Terminator_Object {
  private:
    std::string   m_path;
    long int      m_trials, m_debugstep, m_debuginterval, m_displayinterval;
    double        m_evt_starttime;
    int           m_showtrials;
    bool          m_filter;

    Initialization_Handler * p_inithandler;
    Event_Handler          * p_eventhandler;
    HepMC2_Interface       * p_hepmc2;
    Filter                 * p_filter;

    void DrawLogo(const int mode);  

    void PrepareTerminate();

  public:

    Sherpa();
    ~Sherpa();

    bool InitializeTheRun(int argc,char * argv[]);
    bool InitializeTheEventHandler();

    bool GenerateOneEvent(bool reset=true);
    void FillHepMCEvent(HepMC::GenEvent& event);

    double TotalXS();
    double TotalErr();

    std::string PDFInfo();

    bool SummarizeRun();

    long int NumberOfEvents() const;

    const ATOOLS::Blob_List &GetBlobList() const;

    double GetMEWeight(const ATOOLS::Cluster_Amplitude &ampl,const int mode=0) const;
    
    inline Initialization_Handler * GetInitHandler() const 
    { return p_inithandler; }
    inline Event_Handler * GetEventHandler() const 
    { return p_eventhandler; }
    
  };

  /*!
    \namespace SHERPA
    \brief The namespace for the organisation overhead of the event generator.
  */
  /*!
    \class Sherpa
    \brief Overall steering of the package.

    At the moment, event generation with Sherpa proceeds in two major phases, an initialization
    and, following that, corresponding calls to single event generation. For this, both weighted and
    unweighted events are available.

    The initialization is done through the SHERPA::Initialization_Handler. This object 
    - defines the physical model for the run, 
    - sets up the beams and pdf's, 
    - initializes and evaluates matrix elements for the hard processes, 
    - sets up the parton showers, 
    - the fragmentation including the beam remnant handling

    The single events are organized through the Event_Handler. This object owns a list of
    Event_Phase_Handlers (e.g. Signal_Process, JetEvolution, etc.) that take care of the 
    evolution of the event. This proceeds trough a list of Blob's (objects with n incoming 
    and m outgoing particles), which is passed through the list of Event_Phase_Handler's. Each of 
    them checks, if any "active" blob is of a type the specific Event_Phase_Handler can deal with. 
    If so, the Blob_List and the Blob will be modified accordingly to allow for 
    more modifications in other phases of the event.

    It should be noted that the specific Event_Phase_Handlers and their base class as well as 
    the Event_Handler can be found in the subdirectory SingleEvents, whereas specific handlers 
    that are filled into the Event_Phase_Handlers (such as Matrix_Element_Handler or Shower_Handler) 
    can be found in the subdirectories PerturbativePhysics or SoftPhysics.
  */
  /*!
    \var std::string Sherpa::m_path
    The path for the input files and, eventually, the output.
  */
  /*!
    \var long int Sherpa::m_errors
    The number of errors during event geenration. For each event, up to Sherpa::m_trials
    trials are done. If none leads to a successful event, this counter is incremented.
  */
  /*!
    \var long int Sherpa::m_trials
    The number of trials per single event. At the moment, m_trials = 100.
  */
  /*!
    \var Initialization_Handler * Sherpa::p_inithandler
    A pointer to the Initialization_Handler that manages the initialization of all parts of the 
    package and triggers any caluclation/evaluation (like, e.g. of matrix elements) that has to 
    be done before the simluation of single events. Eventually, the Initialization_Handler
    also deletes all modules after run termination.
  */
  /*!
    \var Event_Handler * Sherpa::p_eventhandler
    A pointer to the handler of single event production.
  */
  /*!
    \var Input_Output_Handler * Sherpa::p_iohandler
    A pointer to the handler of the global output. This module is not very evolved so far, 
    actual input by our clients is needed.
    \todo List all publications that should be cited when running Sherpa.
  */
  /*!
    \var Analysis_Handler * Sherpa::p_analysis
    A pointer to the internal analysis modules.
  */
  /*!
    \fn void Sherpa::DrawLogo()  
    Has to be filled with a nice logo.
  */
  /*!
    \fn Sherpa::Sherpa()
    Constructor of the class, sets the number of allowed trials, m_trials, per single event.
  */
  /*!
    \fn Sherpa::~Sherpa()
    Deletes the following objects: Analysis_Handler, Event_Handler, Output_Handler.
  */
  /*!
    \fn bool Sherpa::InitializeTheRun(std::string _path = std::string("./"))
    Initializes the full run through instantiation of the Initialization_Handler and corresponding
    calls to its methods. If the initialization of the framework was successful, cross sections are
    evaluated, if this worked out, true is returned. If any of the tasks above fails, false is 
    returned.
  */
  /*!
    \fn bool Sherpa::InitializeTheEventHandler()
    The event handler is instatiated, individual event phases are initialized and added to the event 
    handler's list of event phases to be dealt with. Per default, true is returned.
  */
  /*!
    \fn bool Sherpa::GenerateOneEvent()
    In this method, up to m_trials (default 100) trials are made to geenrate one single event. If 
    this succeeds, true is returned, otherwise the result is false. 
  */
  /*!
    \fn bool Sherpa::SummarizeRun()
    Print out event generation statistics and ask the Event_Handler instance to finish.
    Always returns true.
  */
}
#endif

