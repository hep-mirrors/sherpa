#ifndef SHERPA_Single_Events_Hard_Decays_H
#define SHERPA_Single_Events_Hard_Decays_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {
  class Hard_Decay_Handler;

  class Hard_Decays : public Event_Phase_Handler {
  private :
    Hard_Decay_Handler * p_dechandler;
  public :
    Hard_Decays(Hard_Decay_Handler *);
    ~Hard_Decays();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);
  };
}
#endif

