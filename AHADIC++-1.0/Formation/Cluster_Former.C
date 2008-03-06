#include "Cluster_Former.H"
#include "Message.H"

using namespace AHADIC;
using namespace ATOOLS;


Cluster_Former::Cluster_Former() { }

Cluster_Former::~Cluster_Former() { }

void Cluster_Former::ConstructClusters(Proto_Particle_List * plin, Cluster_List * clout)
{
  Cluster * cluster(NULL);
  int       lead(0);
  PPL_Iterator pit1,pit2;
  while (!plin->empty()) {
    pit1=plin->begin();
    pit2 = pit1;pit2++;
    cluster = new Cluster((*pit1),(*pit2));
    clout->push_back(cluster);

    Particle * self = new Particle(-1,Flavour(kf_cluster),cluster->Momentum());
    self->SetNumber();
    self->SetStatus(part_status::active);
    self->SetInfo('C');
    self->SetFinalMass(cluster->Mass());
    control::s_AHAparticles++;
    cluster->SetSelf(self);

    plin->pop_front(); plin->pop_front();
  }
}
