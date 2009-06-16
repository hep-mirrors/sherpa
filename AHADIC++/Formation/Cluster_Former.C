#include "AHADIC++/Formation/Cluster_Former.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Cluster_Former::Cluster_Former() { }

Cluster_Former::~Cluster_Former() { }

void Cluster_Former::ConstructClusters(SP(Proto_Particle_List) plin, Cluster_List * clout)
{
  SP(Cluster) cluster(NULL);
  int       lead(0);
  PPL_Iterator pit1,pit2;
  while (!plin->empty()) {
    pit1 = plin->begin();
    pit2 = pit1;pit2++;
    cluster = new Cluster((*pit1),(*pit2));
#ifdef memchecker
    std::cout<<"@@@ New cluster "<<cluster<<" from "<<METHOD<<"."<<std::endl;
#endif
    clout->push_back(cluster);
    pit1 = plin->erase(pit1);
    pit1 = plin->erase(pit1);
  }
}
