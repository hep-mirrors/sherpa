#ifndef Cluster_Former_H
#define Cluster_Former_H

#include "Cluster.H"

namespace AHADIC {
  class Cluster_Former {
  private:
    void   Reset(); 
  public:
    Cluster_Former();
    ~Cluster_Former();
    void ConstructClusters(Part_List * plin,Cluster_List * clout);
  };
}
#endif
