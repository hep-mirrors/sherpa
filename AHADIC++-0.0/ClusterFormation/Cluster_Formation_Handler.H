#ifndef Cluster_Formation_Handler_H
#define Cluster_Formation_Handler_H

#include "Clusters_2_Hadrons.H"
#include "Gluon_Decayer.H"
#include "Cluster_Former.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Blob.H"


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    double             * p_newmasses, * p_oldens2, * p_ens, m_xmt;
    int                  m_number;
    Part_List          * p_partlist;
    Clusters_2_Hadrons * p_transformer;
    Gluon_Decayer      * p_gludecayer;
    Cluster_Former     * p_cformer;
    
    void  Shift(ATOOLS::Blob_List * bl,ATOOLS::Blob * blob);
  public:
    Cluster_Formation_Handler(Clusters_2_Hadrons *);
    ~Cluster_Formation_Handler();

    ATOOLS::Blob * FormClusters(ATOOLS::Blob_List * bl=NULL);
    Cluster_List * GetClusters() { return p_cformer->GetClusters(); }
  };
}
#endif
