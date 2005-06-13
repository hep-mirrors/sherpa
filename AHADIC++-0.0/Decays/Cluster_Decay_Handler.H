#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Cluster_Decayer_Base.H"
#include "Cluster_Transformer.H"
#include "Blob.H"
#include "Cluster.H"
#include "Histogram.H"
#include <map>


namespace AHADIC {
  struct cdm {
    enum code {
      QoverM_Isotropic        = 11,
      QoverM_Retain           = 12,
      FourFermions_Isotropic  = 21,
      FourFermions_Retain     = 22
    };
  };

  class Cluster_Decay_Handler {
  private:
    cdm::code              m_cdm;
    Cluster_Decayer_Base * p_decayer;
    Cluster_Transformer  * p_transformer;

    Part_List            * p_partlist;
    ATOOLS::Blob         * p_blob;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool DecayIt(Cluster *);
    void InitializeAnalysis();
    void AnalyseThis(ATOOLS::Blob *);
  public:
    Cluster_Decay_Handler(Cluster_Transformer *,bool=false);
    ~Cluster_Decay_Handler();
    void DecayClusters(Cluster_List *,ATOOLS::Blob *);
    Part_List    * GetParticles() { return p_partlist; }
  };
}
#endif
