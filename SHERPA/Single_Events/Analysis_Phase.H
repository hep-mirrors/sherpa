#ifndef SHERPA_Single_Events_Analysis_Phase_H
#define SHERPA_Single_Events_Analysis_Phase_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/Tools/Analysis_Interface.H"

namespace SHERPA {

  class Analysis_Phase: public Event_Phase_Handler {
  private:
    
    Analysis_Vector *p_analyses;

    size_t m_wit;

    std::map<Analysis_Interface*,bool> m_inits;

  public:
    
    // constructors
    Analysis_Phase(Analysis_Vector *const ana);

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp(const size_t & mode=0);

    void Finish(const std::string &);

  };// end of class Analysis_Phase

}// end of namespace SHERPA

#endif
