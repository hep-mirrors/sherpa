#ifndef AHADIC_Decays_Cluster_Part_H
#define AHADIC_Decays_Cluster_Part_H


#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Tools/Cluster_Splitter.H"
#include "ATOOLS/Math/Histogram.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Part {
  protected:
    bool               m_ana;
    Cluster_Splitter * p_csplitter;
    FlavDecayMap       m_options;
    long int           m_fails, m_att;

    std::map<std::string,ATOOLS::Histogram *> m_histograms;
  public:
    Cluster_Part(bool=false);
    ~Cluster_Part();
    bool TestDecay(Cluster * const);
  };
}
#endif
