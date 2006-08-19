#ifndef Hadron_Decays_H
#define Hadron_Decays_H

#include "Event_Phase_Handler.H"
#include "Hadron_Decay_Handler.H"

namespace SHERPA {
  class Hadron_Decays : public Event_Phase_Handler {
  private :
    HDHandlersMap * p_dechandlers;
  public :
    Hadron_Decays(HDHandlersMap *);
    ~Hadron_Decays();
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);
  };
}
#endif

