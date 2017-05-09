#ifndef AHADIC_Formation_Gluon_Splitter_H
#define AHADIC_Formation_Gluon_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Gluon_Splitter : public Splitter_Base {
  private:
    double m_zmin, m_zmax, m_z, m_x, m_y;

    bool      MakeLongitudinalMomenta();
    void      CalculateLimits();
    bool      CalculateXY();

    bool      FillParticlesInLists();
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
