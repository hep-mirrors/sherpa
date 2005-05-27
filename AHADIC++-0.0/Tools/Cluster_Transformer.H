#ifndef Cluster_Transformer_H
#define Cluster_Transformer_H

#include "Transitions.H"
#include "Cluster.H"
#include "Blob.H"

namespace AHADIC {
  struct ctrans {
    enum code {
      nothing      = 0,
      forceddecay  = 1,
      pi0emission  = 2
    };
  };

  class Cluster_Transformer {
  private:
    ctrans::code             m_mode;
    double                   m_offset;
    All_Single_Transitions * p_transitions;

    void DecayCluster(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,ATOOLS::Blob *);
    void DecayCluster(Cluster *,ATOOLS::Flavour &,ATOOLS::Flavour &,Part_List *);
  public:
    Cluster_Transformer(double=0.75);
    ~Cluster_Transformer();
    bool TreatSingleCluster(Cluster_List *,ATOOLS::Blob *);
    bool TreatSingleCluster(Cluster *,Part_List *);
    void TreatClusterList(Cluster_List *,ATOOLS::Blob *);
  };
}
#endif
