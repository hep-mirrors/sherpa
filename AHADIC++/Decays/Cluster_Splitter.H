#ifndef AHADIC_Decays_Cluster_Splitter_H
#define AHADIC_Decays_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    double m_alpha[3], m_beta[3], m_gamma[3];
    size_t m_flcnt;

    double m_z1min, m_z1max, m_z2min, m_z2max, m_R12, m_R21;
    double m_m12min, m_m22min;

    bool MakeLongitudinalMomenta();
    double DeltaM2(const double & MC2,const double & mt2,
		   const ATOOLS::Flavour & flav);
    void CalculateLimits();

    bool FillParticlesInLists();
    bool CheckIfAllowed();
    bool CheckKinematics() { return true; }
    Cluster * MakeCluster(size_t i);
  public:
    Cluster_Splitter(std::list<Cluster *> * cluster_list,
		     Soft_Cluster_Handler * softclusters);
    void Init();
  };
}

#endif
