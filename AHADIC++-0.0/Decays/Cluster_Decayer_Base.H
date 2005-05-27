#ifndef Cluster_Decayer_Base_H
#define Cluster_Decayer_Base_H

#include "Cluster_Part.H"
#include "Hadron_Part.H"
#include "Transitions.H"
#include "Cluster.H"


namespace AHADIC {
  class Cluster_Decayer_Base {
  private:
    All_Single_Transitions * p_stransitions;
    Cluster_Part           * p_cdecs;
    Hadron_Part            * p_chads;
 
    int             m_test;
    double          m_offset;
    ATOOLS::Flavour m_had1, m_had2;
    
    void TestOffSprings(Cluster *);
    void TreatHadDecay(Cluster *,Part_List *);
  public:
    Cluster_Decayer_Base(Cluster_Part *,Hadron_Part *);
    ~Cluster_Decayer_Base();
    bool Treat(Cluster *,Part_List *);
  };
}
#endif
