#ifndef Clusters_2_Hadrons_H
#define Clusters_2_Hadrons_H

#include "Hadron_List.H"
#include "Cluster.H"
#include "Blob.H"
#include "Particle.H"


namespace AHADIC {
  class Clusters_2_Hadrons {
  private:
    double            m_Q, * p_newmasses, * p_oldens2, * p_ens, m_xmt;
    ATOOLS::Flavour * p_newflavs;
    ATOOLS::Vec4D   * p_newmoms, m_cms;
    void Init(int);
    void ShiftBack(size_t);
    void Shift(size_t);
    void AddToBlob(ATOOLS::Blob * blob,const int);
  public:
    Clusters_2_Hadrons();
    ~Clusters_2_Hadrons();
    void Transition(Cluster_List * cl=NULL,ATOOLS::Blob * blob=NULL);
    int  Transition(Cluster * in);
  };
}
#endif
