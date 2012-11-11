#ifndef SHERPA_Single_Events_Output_Phase_H
#define SHERPA_Single_Events_Output_Phase_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/Tools/Output_Base.H"

namespace SHERPA {

  class Event_Handler;

  class Output_Phase: public Event_Phase_Handler {
  private:
    
    Output_Vector *p_outputs;

    size_t m_wit;

    Event_Handler *p_eventhandler;

  public:
    
    // constructors
    Output_Phase(Output_Vector *const ana,Event_Handler *const h);

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);

    void CleanUp(const size_t & mode=0);

    void Finish(const std::string &);

  };// end of class Output_Phase

}// end of namespace SHERPA

#endif
