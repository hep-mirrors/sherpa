#ifndef Hadronization_h
#define Hadronization_h

#include "Event_Phase_Handler.H"
#include "Beam_Remnant_Handler.H"
#include "Fragmentation_Handler.H"

namespace SHERPA {
  class Hadronization : public Event_Phase_Handler {
  private :
    Beam_Remnant_Handler  * p_beamremnanthandler;
    Fragmentation_Handler * p_fragmentationhandler;
  public :
    Hadronization(Beam_Remnant_Handler *,Fragmentation_Handler *);
    ~Hadronization();
    bool Treat(APHYTOOLS::Blob_List *);
    void CleanUp();
  };
}
#endif

