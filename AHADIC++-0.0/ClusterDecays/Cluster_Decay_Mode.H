#ifndef Cluster_Decay_Mode_H
#define Cluster_Decay_Mode_H

#include "Clusters_2_Hadrons.H"
#include "Hadron_List.H"
#include "Cluster.H"
#include "Blob.H"
#include "Particle.H"
#include "Poincare.H"
#include "Message.H"
#include <deque>

namespace AHADIC {
  class Cluster_Decay_Mode {
  protected:
    Clusters_2_Hadrons             * p_transformer;
    std::deque<ATOOLS::Particle *> * p_hadrons;
    FlavCCMap                        m_dicelist;
    double                           m_totalwt;
  public:
    Cluster_Decay_Mode();
    Cluster_Decay_Mode(Clusters_2_Hadrons *);
    virtual ~Cluster_Decay_Mode();

    virtual void Decay(Cluster * cl)=0;
    virtual void Reset();
    virtual void FillHadronsInBlob(ATOOLS::Blob * blob);
  };

  inline void Cluster_Decay_Mode::Reset() {
    if (p_hadrons->empty()) return;
    p_hadrons->clear();
  }
  inline void Cluster_Decay_Mode::FillHadronsInBlob(ATOOLS::Blob * blob) {
    if (p_hadrons->empty()) return;
    while (!p_hadrons->empty()) {
      blob->AddToOutParticles(p_hadrons->back());
      p_hadrons->pop_back();
    } 
    p_hadrons->clear();
  }


  class Simple_Cluster_Fission : public Cluster_Decay_Mode {
  private:
    double            m_Q, * p_masses;
    Cluster        ** p_clusters;
    ATOOLS::Vec4D     m_inmoms[2], m_outmoms[4];
    ATOOLS::Flavour   m_flav;
    WFcomponent       m_options;
    WFcompiter        m_optiter;

    void              SelectFlavour(double maxmass);
    void              SelectHadronPair(Cluster * cluster);
    void              ConstructTestClusters(Cluster * cluster);
    void              IsotropicDecay(Cluster * cluster);
    void              LeadingParticle(ATOOLS::Particle *,double &,double &);
    void              RearrangeMomenta(Cluster * cluster);
  public:
    Simple_Cluster_Fission(Clusters_2_Hadrons *);
    ~Simple_Cluster_Fission();
    
    void Decay(Cluster * cl);
  };
}
#endif
