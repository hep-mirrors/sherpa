#ifndef Analysis_Phase_H
#define Analysis_Phase_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class Analysis_Phase : public Event_Phase_Handler {
    int                          m_status;
    std::string                  m_btype;
    ATOOLS::Primitive_Analysis * p_analysis;

    std::string                  m_path, m_file;
  public:
    Analysis_Phase();
    Analysis_Phase(const std::string & =std::string(""),const std::string & =std::string("./"),const std::string & =std::string("Analysis.dat"));
    Analysis_Phase(ATOOLS::Primitive_Analysis *,const  std::string & =std::string(""));
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string & =std::string("./Results"));
    void SetPrimitiveAnalysis(ATOOLS::Primitive_Analysis * _ana) { p_analysis = _ana; }
    ATOOLS::Primitive_Analysis * GetPrimitiveAnalysis() const    { return p_analysis; }
  };
}

#endif
