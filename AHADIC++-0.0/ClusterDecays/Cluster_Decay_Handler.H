#ifndef Cluster_Decay_Handler_H
#define Cluster_Decay_Handler_H

#include "Clusters_2_Hadrons.H"
#include "Cluster_Decay_Mode.H"
#include "Cluster.H"
#include "Blob.H"
#include "Particle.H"


namespace AHADIC {
  class Cluster_Decay_Handler {
  private:
    int                  m_decaymode;
    Cluster_Decay_Mode * p_decayer;
  public:
    Cluster_Decay_Handler(Clusters_2_Hadrons *);
    ~Cluster_Decay_Handler();

    void DecayThem(Cluster_List * cl=NULL,ATOOLS::Blob * blob=NULL);
  };
}
#endif
