#ifndef MC_Interface_h
#define MC_Interface_h

#include "Event_Phase_Handler.H"
#include "Pythia_Interface.H"

namespace SHERPA {
  class MC_Interface : public Event_Phase_Handler {
  private :
    Pythia_Interface * p_pythia;
  public :
    MC_Interface(Pythia_Interface *);
    ~MC_Interface();

    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);
  };
}
#endif

