#ifndef AHADIC_Decays_Cluster_Splitter_H
#define AHADIC_Decays_Cluster_Splitter_H

#include "AHADIC++/Tools/Splitter_Base.H"

namespace AHADIC {
  class Cluster_Splitter : public Splitter_Base {
  private:
    int    m_mode;
    double m_alpha[4], m_beta[4], m_gamma[4], m_a, m_b, m_c, m_kt02;
    size_t m_flcnt;

    double m_R12, m_R21;
    double m_m12min, m_m22min, m_masses2;
    double m_mean12, m_mean21, m_sigma12, m_sigma21, m_sigma;
    bool   m_output;

    bool   MakeLongitudinalMomenta();
    double DeltaM2(const size_t & cl);
    void   FixCoefficients(const ATOOLS::Flavour & flav1,
			   const ATOOLS::Flavour & flav2);
    void   CalculateLimits();
    bool   FillParticlesInLists();
    bool   CheckIfAllowed();
    bool   CheckKinematics() { return true; }
    Cluster * MakeCluster(size_t i);
  public:
    Cluster_Splitter(std::list<Cluster *> * cluster_list,
		     Soft_Cluster_Handler * softclusters);
    void   Init(const bool & isgluon=false);
    double WeightFunction(const double & z,
			  const double & zmin,const double & zax);
    void   SetOutput(const bool & out) { m_output = out; }
  };
}

#endif
