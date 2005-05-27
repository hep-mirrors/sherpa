#ifndef Cluster_Formation_Handler_H
#define Cluster_Formation_Handler_H

#include "Gluon_Decayer.H"
#include "Cluster_Former.H"
#include "Colour_Reconnections.H"
#include "Cluster_Transformer.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                         m_single_cr,m_double_cr;
    double                       m_offset;

    Gluon_Decayer              * p_gludecayer;
    Cluster_Former             * p_cformer;
    Colour_Reconnections       * p_recons;
    Cluster_Transformer        * p_ctransformer;

    ATOOLS::Blob               * p_blob;
    Cluster_List               * p_clulist;
    std::vector<Part_List *>     m_partlists;
    std::vector<Cluster_List *>  m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void Reset();
    void ExtractSinglets(ATOOLS::Blob_List * bl);
    void FormOriginalClusters();
    void ApplyColourReconnections();
    void ClustersToHadrons();
    void MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(bool=false);
    ~Cluster_Formation_Handler();

    ATOOLS::Blob * FormClusters(ATOOLS::Blob_List * bl=NULL);

    Cluster_Transformer * GetClusterTransformer() { return p_ctransformer; }
    Cluster_List        * GetClusters()           { return p_clulist;  }
    ATOOLS::Blob        * GetFragmentationBlob()  { return p_blob; }
  };
}
#endif
