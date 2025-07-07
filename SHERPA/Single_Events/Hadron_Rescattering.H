#ifndef SHERPA_Single_Events_Hadron_Rescattering_H
#define SHERPA_Single_Events_Hadron_Rescattering_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "ATOOLS/Org/Return_Value.H"

namespace ATOOLS {
  class Blob;
  class Particle;
  class Blob_List;
}

namespace SHERPA {
  class Hadron_Rescattering_Handler;

  class Hadron_Rescattering : public Event_Phase_Handler {
  private :
    Hadron_Rescattering_Handler * p_reschandler;
    
  public :
    Hadron_Rescattering(Hadron_Rescattering_Handler *);
    ~Hadron_Rescattering();
    
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}

#endif
