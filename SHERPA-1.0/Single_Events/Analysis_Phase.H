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
    Analysis_Phase(ANALYSIS::Analysis_Handler *const);

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp();

    void Finish(const std::string &);

    // inline functions
    inline ANALYSIS::Analysis_Handler *const GetAnalysisHandler() const   
    { return p_analysis; }

  };// end of class Analysis_Phase

}// end of namespace SHERPA

#endif
