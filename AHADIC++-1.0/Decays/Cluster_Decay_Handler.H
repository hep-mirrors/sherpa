#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Cluster_Part.H"
#include "Soft_Cluster_Handler.H"
#include "Blob_List.H"
#include "Cluster.H"
#include "Cluster_Decay_Analysis.H"
#include "Return_Value.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    bool                     m_cib;

    Soft_Cluster_Handler   * p_softclusters;
    Cluster_Part           * p_clus;

    Cluster_Decay_Analysis * p_analysis;
 
    ATOOLS::Blob *  DecayIt(Cluster *);
    void            AnalyseThis(ATOOLS::Blob *);
  public:
    Cluster_Decay_Handler(Soft_Cluster_Handler *,bool=false,bool=false);
    ~Cluster_Decay_Handler();

    ATOOLS::Return_Value::code DecayClusters(Cluster_List *,ATOOLS::Blob_List * = NULL);
  };
}
#endif
