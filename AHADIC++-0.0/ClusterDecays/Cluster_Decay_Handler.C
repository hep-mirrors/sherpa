#include "Cluster_Decay_Handler.H"

using namespace AHADIC;
using namespace ATOOLS;

Cluster_Decay_Handler::Cluster_Decay_Handler(Clusters_2_Hadrons * _transformer) :
  m_decaymode(0), p_decayer(NULL)
{
  switch (m_decaymode) {
  default: p_decayer = new Simple_Cluster_Fission(_transformer); break;
  }
}


Cluster_Decay_Handler::~Cluster_Decay_Handler() 
{
  if (p_decayer!=NULL) { delete p_decayer; p_decayer=NULL; }
}

void Cluster_Decay_Handler::DecayThem(Cluster_List * cl,ATOOLS::Blob * blob)
{
  p_decayer->Reset();
  for (Cluster_Iterator cit=cl->begin();cit!=cl->end();cit++) p_decayer->Decay((*cit));
  p_decayer->FillHadronsInBlob(blob);

  int sizecl = cl->size();
  for (int i=0;i<sizecl;i++) { delete cl->back(); cl->pop_back(); }
}
