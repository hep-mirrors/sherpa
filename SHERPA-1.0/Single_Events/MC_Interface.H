#ifndef MC_Interface_h
#define MC_Interface_h

#include "CXXFLAGS.H"

#include "Event_Phase_Handler.H"
#include "Lund_Interface.H"

namespace SHERPA {
  class MC_Interface : public Event_Phase_Handler {
  private :
    Lund_Interface    * p_lund;
    int m_mode;
  public :
    MC_Interface(Lund_Interface *);
    ~MC_Interface();

    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List *, double &);
    void                       CleanUp();
    void                       Finish(const std::string &);
  };
}
#endif

