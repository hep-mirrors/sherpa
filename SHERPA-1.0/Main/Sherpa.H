#ifndef Sherpa_h
#define Sherpa_h


#include "Initialization_Handler.H"
#include "Event_Handler.H"
#include "Output_Handler.H"
#include "Sample_Analysis.H"

namespace SHERPA {
  class Sherpa {
  private:
    std::string   m_path;
    long int      m_errors,m_trials;

    Initialization_Handler * p_inithandler;
    Event_Handler          * p_eventhandler;
    Output_Handler         * p_output;
    Sample_Analysis        * p_analysis;

    void DrawLogo();  

  public:
    Sherpa();
    ~Sherpa();

    bool     InitializeTheRun(std::string);
    bool     InitializeTheEventHandler();
    bool     GenerateOneEvent();
    bool     SummarizeRun();
    long int NumberOfErrors() { return m_errors; } 
  };
}
#endif

