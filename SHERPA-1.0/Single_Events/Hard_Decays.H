#ifndef Hard_Decays_h
#define Hard_Decays_h

#include "Event_Phase_Handler.H"
#include "Hard_Decay_Handler.H"

namespace SHERPA {
  class Hard_Decays : public Event_Phase_Handler {
  private :
    Hard_Decay_Handler * p_dechandler;
  public :
    Hard_Decays(Hard_Decay_Handler *);
    ~Hard_Decays();
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
  };
}
#endif

