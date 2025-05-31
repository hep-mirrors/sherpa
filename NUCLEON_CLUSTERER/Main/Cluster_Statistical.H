#ifndef NUCLEON_CLUSTERER_Main_Cluster_Statistical_H
#define NUCLEON_CLUSTERER_Main_Cluster_Statistical_H

#include "NUCLEON_CLUSTERER/Main/Clusterer_Base.H"

namespace NUCLEON_CLUSTERER {
  class Cluster_Statistical : public Clusterer_Base {
  private:
    std::vector<unsigned int> m_collist;
    int    m_Pmode;
    double m_Q02, m_etaQ, m_R02, m_etaR, m_reshuffle, m_kappa, m_norm;

    void   SetParameters();
    void   FixPMode(const std::string & pm);
    bool   SelectColourPair(const size_t & N,unsigned int & col1, unsigned int & col2);
    bool   AttemptSwap(const unsigned int col[2]);
    void   UpdateColours();
    
    double MomDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double PosDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double ColDistance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
    double TotalLength();
  public:
    Cluster_Statistical();
    ~Cluster_Statistical();

    int    operator()(ATOOLS::Blob_List *const blobs);
    void   Reset();
    double Distance(ATOOLS::Particle * trip,ATOOLS::Particle * anti);
  };
}
#endif
