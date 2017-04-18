#ifndef AHADIC_Decays_Cluster_Splitter_H
#define AHADIC_Decays_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    double m_z1min, m_z1max, m_z2min, m_z2max, m_R12, m_R21;

    bool MakeLongitudinalMomenta();
    void CalculateLimits();

    void FillParticlesInLists();
    bool CheckIfAllowed();
    bool CheckKinematics() { return true; }
    Cluster * MakeCluster(size_t i);
  public:
    Cluster_Splitter(std::list<Cluster *> * cluster_list,
		     Soft_Cluster_Handler * softclusters) :
      Splitter_Base(cluster_list,softclusters) {}  
  };
}

#endif
