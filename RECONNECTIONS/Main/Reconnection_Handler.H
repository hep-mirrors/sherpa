#ifndef RECONNECTIONS_Main_Reconnection_Handler_H
#define RECONNECTIONS_Main_Reconnection_Handler_H

#include "RECONNECTIONS/Main/Reconnection_Base.H"
#include <fstream> // output

namespace RECONNECTIONS {
  typedef std::pair<unsigned int,unsigned int> colpair;
  
  class Reconnection_Handler {
  private:
    bool m_on;
    long unsigned int   m_nfails;
    Reconnection_Base * p_reconnector;
    size_t m_n_variations;
    size_t m_n_mpi_variations;
    double m_weight_cutoff;
    std::vector<size_t> m_cutoff_count;
    
    std::vector<double> m_sum_weights;
    std::vector<double> m_sum_weights_squared;
    size_t m_total_events;
    
    std::ofstream m_weight_file; // output
    
    void AddReconnectionBlob(ATOOLS::Blob_List * blobs);
    void PrintVariationStatistics();
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
