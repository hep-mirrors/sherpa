#ifndef Analysis_Phase_H
#define Analysis_Phase_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Sample_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class Analysis_Phase : public Event_Phase_Handler {
    int               m_status;
    Sample_Analysis * p_analysis;
    std::string       m_outputpath;
  public:
    Analysis_Phase();
    Analysis_Phase(Sample_Analysis *,const  std::string & =std::string(""));
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string & =std::string("./Results"));

    void SetSampleAnalysis(Sample_Analysis * _ana) { p_analysis = _ana; }
    Sample_Analysis * GetSampleAnalysis() const    { return p_analysis; }
  };
}

#endif
