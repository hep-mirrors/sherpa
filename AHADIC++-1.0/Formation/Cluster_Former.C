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
    cluster = new Cluster(pit1->m_flav,pit1->m_mom,
			  pit2->m_flav,pit2->m_mom);
    lead = 0;
    if (pit1->m_info=='L') lead+=1; 
    if (pit2->m_info=='L') lead+=2;
    if (lead>0) cluster->SetLeads(ltp::code(lead));
    clout->push_back(cluster);

    Particle * self = new Particle(-1,Flavour(kf::cluster),cluster->Momentum());
    self->SetNumber();
    self->SetStatus(part_status::active);
    self->SetInfo('C');
    cluster->SetSelf(self);
    plin->pop_front();
    plin->pop_front();
  }
}
