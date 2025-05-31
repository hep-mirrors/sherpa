#ifndef NUCLEON_CLUSTERER_Main_Clusterer_Handler_H
#define NUCLEON_CLUSTERER_Main_Clusterer_Handler_H

#include "NUCLEON_CLUSTERER/Main/Clusterer_Base.H"

namespace NUCLEON_CLUSTERER {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Clusterer_Handler {
  private:
    bool m_on;
    long unsigned int   m_nfails;
    Clusterer_Base * p_reconnector;
    void AddClustererBlob(ATOOLS::Blob_List * blobs);
 public:
    Clusterer_Handler(const bool & on);
    ~Clusterer_Handler();
    
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);

    void Initialize();
    void Reset();
  };
}

#endif
