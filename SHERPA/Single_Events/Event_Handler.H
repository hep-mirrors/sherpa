#ifndef SHERPA_Single_Events_Event_Handler_H
#define SHERPA_Single_Events_Event_Handler_H

#include <vector>
#include <iostream>

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/Tools/Definitions.H"
#include "ATOOLS/Math/Uncertain.H"
#include "ATOOLS/Phys/Variations.H"
#include "ATOOLS/Phys/Weights.H"

namespace SHERPA {
  class Filter;
  
  typedef std::vector<Event_Phase_Handler *> Phase_List;
  typedef Phase_List::iterator Phase_Iterator;

  class Event_Handler {
  private:

    long int m_lastparticlecounter, m_lastblobcounter;
    double m_n, m_addn, m_maxweight;
    double m_mn;
    std::map<std::string,double> m_maxweights;
    int m_checkweight;
    size_t m_lastrss;
    int m_decayer;

    Phase_List        * p_phases;
    ATOOLS::Blob_List   m_blobs, m_sblobs;
    ATOOLS::Blob      * p_signal;
    Filter            * p_filter;

    const ATOOLS::Variations  * p_variations;
    ATOOLS::Weights_Map m_wgtmapsum, m_wgtmapsumsqr;
    ATOOLS::Weights_Map m_mwgtmapsum, m_mwgtmapsumsqr;
    
    void WriteRNGStatus(const std::string &file,
			const std::string &message) const;

    bool GenerateStandardPerturbativeEvent(eventtype::code & mode);
    bool GenerateMinimumBiasEvent(eventtype::code & mode);
    bool GenerateHadronDecayEvent(eventtype::code & mode);
    void InitialiseSeedBlob(ATOOLS::btp::code type,
			    ATOOLS::blob_status::code status);
    bool AnalyseEvent();
    bool WeightsAreGood(const ATOOLS::Weights_Map&);
    int  IterateEventPhases(eventtype::code & mode);

    void MPISync();

    // This entails costly operations and should not be done too often (e.g.
    // only after a successful run).
    void MPISyncXSAndErrMaps();

    // These are only valid after calling MPISyncXSAndErrMaps().
    ATOOLS::Weights_Map TotalXSMPI();
    ATOOLS::Weights_Map TotalErrMPI();

  public:

    Event_Handler();
    ~Event_Handler();

    void AddEventPhase(Event_Phase_Handler *);
    void EmptyEventPhases();

    void SetVariations(const ATOOLS::Variations * p) { p_variations = p; };

    void PrintGenericEventStructure();

    bool GenerateEvent(eventtype::code mode=eventtype::StandardPerturbative);
    void Reset();
    void ResetNonPerturbativePhases();
    void Finish();

    ATOOLS::Blob_List * GetBlobs() { return &m_blobs; }

    void SetFilter(Filter * filter) { p_filter = filter; }
    
    void PerformMemoryMonitoring();

    ATOOLS::Weights_Map TotalXS();
    ATOOLS::Weights_Map TotalErr();
    ATOOLS::Uncertain<double> TotalNominalXS();
    ATOOLS::Uncertain<double> TotalNominalXSMPI();

    std::string CurrentProcess() const;

  };


  /*!
    \file 
    \brief Contains the class SHERPA::Event_Handler
  */

  /*!
    \class Event_Handler
    \brief For the generation of single events.
    
    This class organizes and steers the generation of single events. The idea is the following: 
    In Sherpa, events are represented by blobs which contain a number of incoming and outgoing 
    particles that actually connect the blobs. The Event_Phase_Handler contained in the
    Event_Handler test - one by one - the list of blobs on whether they can treat any of them
    by inspecting the individual blobs type, its status and the like. By treating the blob
    in question, the blob and - potentially - the blob list gets modified. This is repeated
    until no Event_Phase_Handler can deal with any blob of the list. To clearly separate
    this treatment for the time being the Event_Phase_Handlers are divided into two types,
    (either "Perturbative" or "Hadronization") and the procedure outlined above is performed
    for the blobs of the perturbative kind first, before the phases of type "Hadronization" 
    take over. The different Event_Phases are filled in by the Sherpa class through the method 
    AddEventPhase (Event_Phase_Handler *), remember, this is the actual steering class.
  */ 
  /*!
    \var Phase_List * Event_Handler::p_phases
    The list of Event_Phase_Handlers. Through a typedef this list is realized by a vector from the
    standard library, the corresponding iterator is defined to be Phase_Iterator.
  */ 
  /*!
    \var ATOOLS::Blob_List Event_Handler::m_blobs
    The internal list of blobs constituting the event. This list is emptied before any new event
    is generated.
  */ 
  /*!
    \fn Event_Handler::Event_Handler()
    Instantiates the phase list.
  */ 
  /*!
    \fn Event_Handler::~Event_Handler()
    Deletes the phase list.
  */ 
  /*!
    \fn void Event_Handler::AddEventPhase(Event_Phase_Handler *)
    Checks if the particular Event_Phase has already been included. If so, it won't be added
    another time. in the way, Sherpa organizes event generation this would have no sense
    anyhow. 
  */ 
  /*!
    \fn void Event_Handler::EmptyEventPhases()
    Deletes all individual Event_Phase_Handlers and empties the list of event phases.
  */ 
  /*!
    \fn void Event_Handler::SetVariations()
    Set the Variations object for the entire event generation.
  */ 
  /*!
    \fn bool Event_Handler::GenerateEvent()
    Deletes all blobs that might have been produced but not deleted before. Then, a first
    empty blob is created of type "Signal Process". This blob is put into the blob list
    and the list is then dealt with by the event phases. In fact, the first blob will
    be filled by the SignalProcess phase which contains the matrix element generator.
    It is foreseen that pile-up is steered by generating a fluctuating number of
    further blobs of type, say "Minimum bias" that will then be filled by an anticipated
    PileUp Event_Phase_Handler. However, having filled the signal blob this method iterates 
    through all event phases of the type "Perturbative" until no blob inside the blob list 
    can be dealt with by any of them. If this is the case, the same procdure is
    repeated with the phases of type "Hadronization". 
  */ 
  /*!
    \fn void Event_Handler::PrintGenericEventStructure()
    Prints out all Event_Phase_Handlers used for event generation. This might give some
    idea of what to expect.
  */ 
}
#endif

