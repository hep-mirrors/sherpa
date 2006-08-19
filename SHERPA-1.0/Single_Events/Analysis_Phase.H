#ifndef Analysis_Phase_H
#define Analysis_Phase_H

#include "Event_Phase_Handler.H"

namespace ANALYSIS { class Analysis_Handler; }

namespace SHERPA {

  class Analysis_Phase: public Event_Phase_Handler {
  private:

    ANALYSIS::Analysis_Handler *p_analysis;

  public:
    
    // constructors
    Analysis_Phase(ANALYSIS::Analysis_Handler *const analysis);

    // member functions
    bool Treat(ATOOLS::Blob_List *bloblist,double &weight);
    void CleanUp();

    void Finish(const std::string &path);

    // inline functions
    inline ANALYSIS::Analysis_Handler *const GetAnalysisHandler() const   
    { return p_analysis; }

  };// end of class Analysis_Phase

}// end of namespace SHERPA

#endif
