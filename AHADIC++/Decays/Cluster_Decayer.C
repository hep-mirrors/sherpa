#include "AHADIC++/Decays/Cluster_Decayer.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decayer::Cluster_Decayer(list<Cluster *> * cluster_list,
			     Soft_Cluster_Handler * softclusters) :
  p_cluster_list(cluster_list), p_softclusters(softclusters),
  m_splitter(Cluster_Splitter(cluster_list,softclusters))
{}

Cluster_Decayer::~Cluster_Decayer() {
}

void Cluster_Decayer::Init() {
  m_splitter.Init();
}

bool Cluster_Decayer::operator()() {
  while (!p_cluster_list->empty()) {
    if (!Treat(p_cluster_list->front())) return false;
    p_cluster_list->pop_front();
  }
  return true;
}

bool Cluster_Decayer::Treat(Cluster * cluster) {
  if (!(p_softclusters->MustPromptDecay(cluster)) &&
      m_splitter((*cluster)[0],(*cluster)[1])) {
    delete cluster;
    return true;
  }
  if (p_softclusters->Treat(cluster,true)) {
    cluster->Clear();
    delete cluster;
    return true;
  }
  msg_Error()<<METHOD<<" throws error for:\n"<<(*cluster)<<"\n";
  return false;
}
