#ifndef RECONNECTIONS_Main_Reconnection_Handler_H
#define RECONNECTIONS_Main_Reconnection_Handler_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"

namespace RECONNECTIONS {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Reconnection_Handler {
  private:
    bool m_on;
    long unsigned int   m_nfails;
    Reconnection_Base * p_reconnector;
    void AddReconnectionBlob(ATOOLS::Blob_List * blobs);
 public:
    Reconnection_Handler(const bool & on);
    ~Reconnection_Handler();
    
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);

    void Initialize();
    void Reset();
  };
}

#endif
