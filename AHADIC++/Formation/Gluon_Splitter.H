#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    double m_z1min, m_z1max, m_z2min, m_z2max, m_x, m_y;

    bool      MakeLongitudinalMomenta();
    void      CalculateLimits();
    bool      CalculateXY();

    void      FillParticlesInLists();
    bool      CheckKinematics();
    Cluster * MakeCluster();
    void      UpdateSpectator();
  public:
    Gluon_Splitter(std::list<Cluster *> * cluster_list,
		   Soft_Cluster_Handler * softclusters) :
    Splitter_Base(cluster_list,softclusters) {}  
  };
}

#endif
