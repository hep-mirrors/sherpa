#ifndef NUCLEON_CLUSTERER_Main_Cluster_By_Singlet_H
#define NUCLEON_CLUSTERER_Main_Cluster_By_Singlet_H

#include "NUCLEON_CLUSTERER/Main/Clusterer_Base.H"
#include "NUCLEON_CLUSTERER/Main/Clusterer_Weights.H"

namespace NUCLEON_CLUSTERER {
  class Cluster_By_Singlet : public Clusterer_Base {
  private:    
    Clusterer_Weights m_weights;
    int    m_Pmode;
    double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_restring;
    std::list<ATOOLS::Part_List *> m_singlets;    
    
    void SetParameters();
    void MakeSinglets();
    void FlattenSinglets();
    ATOOLS::Particle * FindStart();
    ATOOLS::Particle * FindNext(const size_t & col);
    void ReshuffleSinglets();
    bool ReshuffleSinglet(ATOOLS::Part_List * singlet);
    void ClusterSinglets();
    void SpliceSinglets(ATOOLS::Part_List * sing1,ATOOLS::Part_List * sing2,
			ATOOLS::Part_Iterator & pit1,ATOOLS::Part_Iterator & pit2);
    void AftermathOfSlicing(ATOOLS::Particle * part11,ATOOLS::Particle * part12,
			    ATOOLS::Particle * part21,ATOOLS::Particle * part22);
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  public:
    Cluster_By_Singlet();
    ~Cluster_By_Singlet();

    int    operator()(ATOOLS::Blob_List *const blobs);
    void   Reset();
    double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  };
}
#endif
