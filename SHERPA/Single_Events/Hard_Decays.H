#ifndef SHERPA_Single_Events_Hard_Decays_H
#define SHERPA_Single_Events_Hard_Decays_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"

namespace SHERPA {
  class Decay_Handler_Base;

  class Hard_Decays : public Event_Phase_Handler {
  private :
    Decay_Handler_Base* p_dechandler;
  public :
    Hard_Decays(Decay_Handler_Base *);
    ~Hard_Decays();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}
#endif

