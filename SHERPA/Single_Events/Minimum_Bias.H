#ifndef SHERPA_Single_Events_Minimum_Bias_H
#define SHERPA_Single_Events_Minimum_Bias_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"

namespace PDF    { class Remnant_Base; }

namespace SHERPA {

  class Minimum_Bias: public Event_Phase_Handler {
  private:
    Soft_Collision_Handler * p_schandler;
    bool m_on;
  public:
    Minimum_Bias(Soft_Collision_Handler_Map * schandlers);
    ~Minimum_Bias();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List* blobs);
    void                       CleanUp(const size_t & mode=0);
    void                       Finish(const std::string &);
  }; // end of class Minimum_Bias

} // end of namespace SHERPA

#endif

