#ifndef MC_Interface_h
#define MC_Interface_h

#include "Event_Phase_Handler.H"
#include "Lund_Interface.H"
#ifdef USING__MCatNLO
#include "Herwig_Interface.H"
#include "MCatNLO_Interface.H"
#endif

namespace SHERPA {
  class MC_Interface : public Event_Phase_Handler {
  private :
    Lund_Interface    * p_pythia;
#ifdef USING__MCatNLO
    Herwig_Interface  * p_herwig;
    MCatNLO_Interface * p_mcatnlo;
#endif
    int m_mode;
  public :
    MC_Interface(Lund_Interface *);
#ifdef USING__MCatNLO
    MC_Interface(Herwig_Interface *);
    MC_Interface(MCatNLO_Interface *);
#endif
    ~MC_Interface();

    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
    void Finish(const std::string &);
  };
}
#endif

