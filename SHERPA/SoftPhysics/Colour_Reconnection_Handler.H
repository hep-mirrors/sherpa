#ifndef SHERPA_SoftPhysics_Colour_Reconnection_Handler_H
#define SHERPA_SoftPhysics_Colour_Reconnection_Handler_H

#include "RECONNECTIONS/Main/Reconnection_Handler.H"
#include "ATOOLS/Org/Return_Value.H"

namespace SHERPA {
  class Colour_Reconnection_Handler {
  private:
    bool                                  m_on;
    RECONNECTIONS::Reconnection_Handler * p_reconnections;
    ATOOLS::Blob * p_blob;
  public :
    Colour_Reconnection_Handler();
    ~Colour_Reconnection_Handler();
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const);
    void CleanUp(const size_t & mode=0);
    void Output();
  };
}

#endif
