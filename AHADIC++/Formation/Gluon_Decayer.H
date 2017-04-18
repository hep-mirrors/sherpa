#ifndef AHADIC_Formations_Gluon_Decayer_H
#define AHADIC_Formations_Gluon_Decayer_H

#include "AHADIC++/Formation/Gluon_Splitter.H"
#include "AHADIC++/Formation/Trivial_Splitter.H"
#include "AHADIC++/Tools/Singlet_Tools.H"
#include "AHADIC++/Tools/Soft_Cluster_Handler.H"
#include "AHADIC++/Tools/Proto_Particle.H"
#include <list>

namespace AHADIC {
  class Gluon_Decayer : public Singlet_Tools {
  private:
    std::list<Cluster *> * p_cluster_list; 
    Soft_Cluster_Handler * p_softclusters;
    Gluon_Splitter   m_splitter;
    Trivial_Splitter m_breaker;
    
    Proto_Particle * p_part1, * p_part2;

    bool             SplitGluonRing();
    Proto_Particle * FirstGluon();

    bool Trivial(Proto_Particle * part1,Proto_Particle * part2);
    int  Step(Proto_Particle * part1,Proto_Particle * part2,
	      Proto_Particle * part3=NULL);
    bool LastStep();
  public:
    Gluon_Decayer(std::list<Cluster *> * cluster_list,
		  Soft_Cluster_Handler * softclusters);
    ~Gluon_Decayer();
    void Init();
    
    bool operator()(Singlet * singlet);
  };
}

#endif
