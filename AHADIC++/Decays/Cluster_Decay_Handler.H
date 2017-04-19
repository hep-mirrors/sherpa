#ifndef AHADIC_Decays_Cluster_Decay_Handler_H
#define AHADIC_Decays_Cluster_Decay_Handler_H

#include "AHADIC++/Decays/Cluster_Part.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "AHADIC++/Tools/Cluster.H"
#include "AHADIC++/Decays/Cluster_Decay_Analysis.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    Soft_Cluster_Handler   * p_softclusters;
    Cluster_Part           * p_clus;
    Cluster_List           * p_clulist;

    Cluster_Decay_Analysis * p_analysis;
 
    ATOOLS::Blob * DecayIt(Cluster *);
  public:
    Cluster_Decay_Handler(Cluster_List *,bool=false);
    ~Cluster_Decay_Handler();

    int DecayClusters(ATOOLS::Blob *);
  };
}
#endif
