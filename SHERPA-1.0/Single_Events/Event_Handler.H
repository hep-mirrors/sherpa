#ifndef Event_Handler_h
#define Event_Handler_h

#include <vector>
#include <iostream>

#include "Event_Phase_Handler.H"
#include "Sample_Analysis.H"

namespace SHERPA {
  typedef std::vector<Event_Phase_Handler *> Phase_List;
  typedef Phase_List::iterator Phase_Iterator;

  class Event_Handler {
    Phase_List             * p_phases;
    Sample_Analysis        * p_analysis;
    APHYTOOLS::Blob_List     m_blobs;
    APHYTOOLS::Parton_List   m_partons;

    void PrintBlobs();
    void PrintPartons();

    void CleanUpEvent();
  public:
    Event_Handler();
    ~Event_Handler();

    void AddEventPhase(Event_Phase_Handler *);
    void EmptyEventPhases();
    void PrintGenericEventStructure();

    bool GenerateEvent();
    void PerformAnalysis();
    void PrintEvent(int);

    void SetAnalysis(Sample_Analysis *);
  };
}
#endif

