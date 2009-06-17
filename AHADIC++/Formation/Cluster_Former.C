#include "AHADIC++/Formation/Cluster_Former.H"
#include "ATOOLS/Org/Message.H"

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
    pit1 = plin->begin();
    pit2 = pit1;pit2++;
    cluster = new Cluster((*pit1),(*pit2));
#ifdef memchecker
    std::cout<<"*** New cluster ("<<cluster->Number()<<"/"<<cluster<<") with "
	     <<"pps ("<<(*pit1)<<"/"<<(*pit2)<<") from "<<METHOD<<"."<<std::endl;
#endif
    clout->push_back(cluster);
    (*pit1)->m_kt2max = (*pit2)->m_kt2max = ATOOLS::Max((*pit1)->m_kt2max,(*pit2)->m_kt2max);
    if (IsZero((*pit1)->m_kt2max)) 
      (*pit1)->m_kt2max = (*pit2)->m_kt2max = 
	cluster->GetTrip()->m_mom.PPerp2(cluster->GetAnti()->m_mom);
    if (IsZero((*pit1)->m_kt2max)) 
      (*pit1)->m_kt2max = (*pit2)->m_kt2max = 
	cluster->Mass2()-sqr(cluster->GetTrip()->m_flav.HadMass()+
			     cluster->GetAnti()->m_flav.HadMass());
    pit1 = plin->erase(pit1);
    pit1 = plin->erase(pit1);

  }
  EstablishRelations(clout);
}

void Cluster_Former::EstablishRelations(Cluster_List * clist) {
  Cluster * cluster(NULL);
  Proto_Particle * hook;
  Cluster_Iterator clu,check;
  for (clu=clist->begin();clu!=clist->end();clu++) {
    cluster = (*clu);
    hook = cluster->GetTrip();
    for (check=clist->begin();check!=clist->end();check++) {
      if ((*check)==cluster) continue;
      if (hook->p_partner==(*check)->GetAnti()) {
	cluster->SetNBTrip((*check));
	(*check)->SetNBAnti(cluster);
      }
    }
    hook = cluster->GetAnti();
    for (check=clist->begin();check!=clist->end();check++) {
      if ((*check)==cluster) continue;
      if (hook->p_partner==(*check)->GetTrip()) {
	cluster->SetNBAnti((*check));
	(*check)->SetNBTrip(cluster);
      }
    }
  }
}
