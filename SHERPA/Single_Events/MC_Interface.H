#ifndef SHERPA_Single_Events_MC_Interface_H
#define SHERPA_Single_Events_MC_Interface_H

#include "ATOOLS/Org/CXXFLAGS.H"

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/LundTools/Lund_Interface.H"

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

