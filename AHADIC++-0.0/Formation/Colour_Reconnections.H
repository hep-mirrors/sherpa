#ifndef Colour_Reconnections_H
#define Colour_Reconnections_H

#include "Cluster.H"

namespace AHADIC {
  class Colour_Reconnections {
  private:
    int    m_kinmode, m_colmode;
    double m_t0,m_w12,m_w34,m_w14,m_w23;

    bool   TestClusters(Cluster * cl1,Cluster * cl2,int gen=1);
    double KinematicWeight(const ATOOLS::Vec4D & mom1,const ATOOLS::Vec4D & mom2,
			   const ATOOLS::Vec4D & mom3,const ATOOLS::Vec4D & mom4);
    double ColourWeight(int gen=1);
  public:
    Colour_Reconnections();
    ~Colour_Reconnections();
    void   Singlet_CR(Cluster_List * clin);
    void   Two_Singlet_CR(Cluster_List * cl1,Cluster_List * cl2);
  };
}
#endif
