#ifndef Analysis_Phase_H
#define Analysis_Phase_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Sample_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class Analysis_Phase : public Event_Phase_Handler {
    int                m_status;
    Analysis_Handler * p_analysis;
    std::string        m_outputpath;
  public:
    Analysis_Phase();
    Analysis_Phase(Analysis_Handler *,const  std::string & =std::string(""));
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string & =std::string("./Results"));

    void SetAnalysisHandler(Analysis_Handler * _ana) { p_analysis = _ana; }
    Analysis_Handler * GetAnalysisHandler() const   { return p_analysis; }
  };
}

#endif
