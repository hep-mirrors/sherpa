#ifndef AHADIC_Decays_Gluon_Decayer_H
#define AHADIC_Decays_Gluon_Decayer_H

#include "AHADIC++/Decays/Cluster_Splitter.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Proto_Particle.H"
#include <list>

namespace AHADIC {
  class Cluster_Decayer {
  private:
    std::list<Cluster *> * p_cluster_list; 
    Soft_Cluster_Handler * p_softclusters;
    Cluster_Splitter m_splitter;

    bool Treat(Cluster * cluster);
  public:
    Cluster_Decayer(std::list<Cluster *> * cluster_list,
		    Soft_Cluster_Handler * softclusters);
    ~Cluster_Decayer();
    void Init();
    
    bool operator()();
  };
}

#endif
