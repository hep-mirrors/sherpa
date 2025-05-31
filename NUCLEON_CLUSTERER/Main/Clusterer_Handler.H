#ifndef NUCLEON_CLUSTERER_Main_Clusterer_Handler_H
#define NUCLEON_CLUSTERER_Main_Clusterer_Handler_H

#include "NUCLEON_CLUSTERER/Main/Clusterer_Base.H"

namespace NUCLEON_CLUSTERER {
  // replace if nulceons interact, i.e. via hadron-hadron interactions
  typedef std::pair<unsigned int,unsigned int> nucleonpair;
  
  class Clusterer_Handler {
  private:
    bool m_on; // clustering on/off
    long unsigned int   m_nfails; // number of failed clustering attempts
    Clusterer_Base * p_clusterer; // pointer to the clusterer object
    void AddClustererBlob(ATOOLS::Blob_List * blobs); // add cluster data to blobs
  public:
    // Constructor and Destructor
    Clusterer_Handler(const bool & on);
    ~Clusterer_Handler();
    // perform clustering on input blobs and particle list (maybe make into nucelon list)
    ATOOLS::Return_Value::code operator()(ATOOLS::Blob_List *const blobs,
					  ATOOLS::Particle_List *const=NULL);
    // lifecycle management
    void Initialize();
    void Reset();
  };
}

#endif
