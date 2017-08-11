#ifndef AHADIC_Formation_Cluster_Formation_Handler_H
#define AHADIC_Formation_Cluster_Formation_Handler_H

#include "AHADIC++/Formation/Gluon_Decayer.H"
#include "AHADIC++/Formation/Cluster_Former.H"
#include "AHADIC++/Formation/Colour_Reconnections.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Cluster.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Histogram.H"
#include <map>


namespace AHADIC {
  class Cluster_Formation_Handler {
  private:
    bool                                m_single_cr,m_double_cr;

    Gluon_Decayer                     * p_gludecayer;
    Cluster_Former                    * p_cformer;
    Colour_Reconnections              * p_recons;
    Soft_Cluster_Handler              * p_softclusters;

    Cluster_List                      * p_clulist;
    ListOfPPLs                          m_partlists;
    std::vector<Cluster_List*>          m_clulists;

    bool                                      m_analyse;
    std::map<std::string,ATOOLS::Histogram *> m_histograms;

    bool ShiftOnMassShells();
    bool ShiftList(Proto_Particle_List *);
    Proto_Particle_List * SelectFromList(ListOfPPLs * lppl,
					 Proto_Particle_List * ppl=NULL);
    bool ExtractSinglets(ATOOLS::Blob *);
    void SetInfoTagForPrimaryParticle(Proto_Particle * proto) const;
    bool FormOriginalClusters();
    bool ApplyColourReconnections();
    bool ClustersToHadrons(ATOOLS::Blob *);
    bool MergeClusterListsIntoOne();
  public:
    Cluster_Formation_Handler(Cluster_List*, bool=false);
    ~Cluster_Formation_Handler();

    int FormClusters(ATOOLS::Blob *);
    void Reset();
  };
}
#endif
