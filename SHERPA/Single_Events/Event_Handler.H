#ifndef SHERPA_Single_Events_Event_Handler_H
#define SHERPA_Single_Events_Event_Handler_H

#include <vector>
#include <iostream>

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {

  typedef std::vector<Event_Phase_Handler *> Phase_List;
  typedef Phase_List::iterator Phase_Iterator;

  class Event_Handler {
  private:

    long int m_lastparticlecounter, m_lastblobcounter;
    double m_n, m_addn, m_sum, m_sumsqr;

    Phase_List        * p_phases;
    ATOOLS::Blob_List   m_blobs, m_sblobs;
    ATOOLS::Blob      * p_signal;

  public:

    Event_Handler();
    ~Event_Handler();

    void                  AddEventPhase(Event_Phase_Handler *);
    Event_Phase_Handler * GetEventPhase(const size_t);
    size_t                NumberOfEventPhases();
    void                  EmptyEventPhases();

    void PrintGenericEventStructure();

    bool GenerateEvent(int=0);
    void Reset();
    void Finish();

    ATOOLS::Blob_List * GetBlobs() { return &m_blobs; }

    inline double TotalXS() const { return m_sum/m_n; }
    inline double TotalVar() const
    { return (m_sumsqr-m_sum*m_sum/m_n)/(m_n-1); }
    inline double TotalErr() const {
      if (m_n==1) return TotalXS();
      else if (ATOOLS::IsEqual
	       (m_sumsqr*m_n,m_sum*m_sum,1.0e-6)) return 0.0;
      else return sqrt(TotalVar()/m_n);
    }

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
    \fn void Event_Handler::CleanUpEvent()
    Cleans up in all event phases - this might be important if they change some of their
    internal parameters during an event that have to be resetted - , deletes all internal blobs
    and empties the blob list.
  */ 
  /*!
    \fn void Event_Handler::PrintGenericEventStructure()
    Prints out all Event_Phase_Handlers used for event generation. This might give some
    idea of what to expect.
  */ 
  /*!
    \fn void Event_Handler::PrintEvent(int)
    Prints out the blobs. A second option should be added printing out some particle-based
    event record. 
    \todo Enable to print out partons instead of blobs.
  */ 
  /*!
    \fn void Event_Handler::PrintBlobs()
    Prints the actual blob list.
  */ 
}
#endif

