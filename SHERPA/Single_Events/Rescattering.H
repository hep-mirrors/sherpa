#ifndef SHERPA_Single_Events_Rescattering_H
#define SHERPA_Single_Events_Rescattering_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/PerturbativePhysics/Rescattering_Handler.H"

namespace SHERPA {
  class Rescattering_Handler;
  
  class Rescattering : public Event_Phase_Handler {
  private :
    bool m_on;
    Rescattering_Handler * p_reschandler;
  public :
    Rescattering(Rescattering_Handler * reschandler=NULL);
    ~Rescattering();
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}

#endif
