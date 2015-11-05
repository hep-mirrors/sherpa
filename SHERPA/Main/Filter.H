#ifndef SHERPA_Filter_H
#define SHERPA_Filter_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Phys/Blob_List.H"

namespace SHERPA {
  class Filter {
  private:
    size_t m_N_charged_min, m_N_charged_max;
    double m_etamax, m_ptmin;

    size_t NumberOfChargedFSParticles(ATOOLS::Blob_List * blobs) const;
  public:
    Filter(); 
    ~Filter();

    bool operator()(ATOOLS::Blob_List * blobs) const;
  };
}

#endif
