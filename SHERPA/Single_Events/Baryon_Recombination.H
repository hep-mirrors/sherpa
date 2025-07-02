#ifndef SHERPA_Single_Events_Baryon_Recombination_H
#define SHERPA_Single_Events_Baryon_Recombination_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "ATOOLS/Org/Return_Value.H"

namespace ATOOLS {
  class Blob;
  class Particle;
  class Blob_List;
}

namespace SHERPA {
  class Baryon_Recombination_Handler;
  
  class Baryon_Recombination : public Event_Phase_Handler {
  private :
    Baryon_Recombination_Handler * p_BRhandler;
  public :
    Baryon_Recombination(Baryon_Recombination_Handler * BRhandler);
    ~Baryon_Recombination();
    
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List*);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  };
}
#endif

