#ifndef Cluster_Decay_Analysis_H
#define Cluster_Decay_Analysis_H

#include "Blob.H"
#include "Cluster.H"
#include "Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Decay_Analysis {
  private:
    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;
  public:
    Cluster_Decay_Analysis();
    ~Cluster_Decay_Analysis();

    void AnalyseThis(ATOOLS::Blob *);
  };
}
#endif
