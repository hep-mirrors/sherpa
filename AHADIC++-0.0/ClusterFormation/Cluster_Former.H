#ifndef Cluster_Former_H
#define Cluster_Former_H

#include "Cluster.H"
#include "Particle.H"

namespace AHADIC {
  typedef std::pair<ATOOLS::Particle *,ATOOLS::Particle *> precluster;

  class Cluster_Former {
  private:
    int                       m_kinmode, m_colmode;
    double                    m_t0;
    Part_List               * p_pl;
    Cluster_List            * p_cl;
    std::vector<precluster *> m_preclusters;

    void   Reset(); 
    void   ConstructPreClusters();
    void   ReshuffleClusters();
    void   FormFinalClusters();
    double KinematicWeight(const ATOOLS::Vec4D & mom1,const ATOOLS::Vec4D & mom2,
			   const ATOOLS::Vec4D & mom3,const ATOOLS::Vec4D & mom4);
    double ColourWeight();
  public:
    Cluster_Former();
    ~Cluster_Former();
    void   FormClusters(Part_List * pl);
    Cluster_List * GetClusters() { return p_cl; }
  };
}
#endif
