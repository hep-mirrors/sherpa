#ifndef SHERPA_Single_Events_Analysis_Phase_H
#define SHERPA_Single_Events_Analysis_Phase_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {

  class Analysis_Interface;

  class Analysis_Phase: public Event_Phase_Handler {
  private:
    
    Analysis_Interface *p_analysis;

    size_t m_wit, m_init;

  public:
    
    // constructors
    Analysis_Phase(Analysis_Interface *const ana);

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp();

    void Finish(const std::string &);

    // inline functions
    inline void *const GetAnalysisHandler() const   
    { return p_analysis; }

  };// end of class Analysis_Phase

}// end of namespace SHERPA

#endif
