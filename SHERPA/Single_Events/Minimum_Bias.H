#ifndef SHERPA_Single_Events_Minimum_Bias_H
#define SHERPA_Single_Events_Minimum_Bias_H

#include "SHERPA/Single_Events/Event_Phase_Handler.H"
#include "SHERPA/SoftPhysics/Soft_Collision_Handler.H"

namespace PDF    { class Remnant_Base; }

namespace SHERPA {

  class Minimum_Bias: public Event_Phase_Handler {
  private:
    double m_weight;

    Soft_Collision_Handler * p_schandler;
    PDF::Remnant_Base      * p_remnants[2];
  public:
    Minimum_Bias(Soft_Collision_Handler * schandler);
    ~Minimum_Bias();

    // member functions
    ATOOLS::Return_Value::code Treat(ATOOLS::Blob_List * blobs, 
				     double & weight);
    void                       CleanUp();
    void                       Finish(const std::string &);
  }; // end of class Minimum_Bias

} // end of namespace SHERPA

#endif

