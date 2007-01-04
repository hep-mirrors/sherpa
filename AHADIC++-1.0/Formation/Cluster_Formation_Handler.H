#ifndef Cluster_Formation_Handler_H
#define Cluster_Formation_Handler_H

#include "Gluon_Decayer.H"
#include "Cluster_Former.H"
#include "Colour_Reconnections.H"
#include "Soft_Cluster_Handler.H"
#include "Cluster.H"
#include "Blob_List.H"
#include "Blob.H"
#include "Histogram.H"
#include "Return_Value.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                                m_single_cr,m_double_cr;

    Gluon_Decayer                     * p_gludecayer;
    Cluster_Former                    * p_cformer;
    Colour_Reconnections              * p_recons;
    Soft_Cluster_Handler              * p_softclusters;

    ATOOLS::Blob                      * p_blob;
    Cluster_List                      * p_clulist;
    VectorOfPPLs                        m_partlists;
    std::vector<Cluster_List *>         m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    void Reset();
    ATOOLS::Return_Value::code ExtractSinglets();
    ATOOLS::Return_Value::code FormOriginalClusters();
    ATOOLS::Return_Value::code ApplyColourReconnections();
    ATOOLS::Return_Value::code ClustersToHadrons(ATOOLS::Blob_List *);
    ATOOLS::Return_Value::code MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(bool=false);
    ~Cluster_Formation_Handler();

    ATOOLS::Return_Value::code FormClusters(ATOOLS::Blob *,
					    ATOOLS::Blob_List *);

    //! general access routines
    //@{
    Soft_Cluster_Handler * GetSoftClusterHandler() { return p_softclusters; }
    Cluster_List         * GetClusters()           { return p_clulist;  }
    ATOOLS::Blob         * GetFragmentationBlob()  { return p_blob; }
    //@}
  };
}
#endif
