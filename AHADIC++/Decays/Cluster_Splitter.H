#ifndef AHADIC_Decays_Cluster_Splitter_H
#define AHADIC_Decays_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    int    m_mode;
    double m_alpha[4], m_beta[4], m_gamma[4];
    double m_a[2], m_b[2], m_c[2], m_kt02;

    double m_R2[2], m_mean[2], m_sigma[2], m_m2min[2], m_masses;
    bool   m_output;

    bool   MakeLongitudinalMomenta();
    bool   MakeLongitudinalMomentaZ();
    bool   MakeLongitudinalMomentaZSimple();
    bool   MakeLongitudinalMomentaMass();
    bool   MakeLongitudinalMomentaMassSimple();
    double DeltaM(const size_t & cl);
    void   FixCoefficients();
    void   CalculateLimits();
    bool   FillParticlesInLists();
    bool   RecalculateZs();
    bool   CheckKinematics() { return true; }
    Cluster * MakeCluster(size_t i);
  public:
    Cluster_Splitter(std::list<Cluster *> * cluster_list,
		     Soft_Cluster_Handler * softclusters);
    void   Init(const bool & isgluon=false);
    double WeightFunction(const double & z,
			  const double & zmin,const double & zmax,
			  const unsigned int & cnt);
    void   SetOutput(const bool & out) { m_output = out; }
  };
}

#endif
