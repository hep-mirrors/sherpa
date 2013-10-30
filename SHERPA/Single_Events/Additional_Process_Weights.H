#ifndef SHERPA_Single_Events_Additional_Process_Weights_H
#define SHERPA_Single_Events_Additional_Process_Weights_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {

  class Additional_Process_Weights: public Event_Phase_Handler {
  private:
    bool m_treated;
    ATOOLS::Blob * p_signal;
  public:
    Additional_Process_Weights();
    ~Additional_Process_Weights();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *,double &);
    void CleanUp(const size_t & mode=0);
    void Finish(const std::string &);
  };

}

#endif

