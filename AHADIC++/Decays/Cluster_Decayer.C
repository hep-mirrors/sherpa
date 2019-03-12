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

Cluster_Decayer::~Cluster_Decayer() {}

void Cluster_Decayer::Init() { m_splitter.Init(); }

void Cluster_Decayer::Reset() {}

bool Cluster_Decayer::operator()() {
  //if (!p_cluster_list->empty()) msg_Out()<<METHOD<<"\n"<<(*p_cluster_list)<<"\n\n";
  while (!p_cluster_list->empty()) {
    if (!Treat(p_cluster_list->front())) {
      return false;
    }
    p_cluster_list->pop_front();
  }
  return true;
}

bool Cluster_Decayer::Treat(Cluster * cluster) {
  bool mustdecay = p_softclusters->MustPromptDecay(cluster);
  if (!mustdecay && m_splitter((*cluster)[0],(*cluster)[1])) {
    delete cluster;
    return true;
  } 
  switch (p_softclusters->Treat(cluster,true)) {
  case -1:
    // cluster cannot decay into anything - return false (triggers new event)
    msg_Error()<<METHOD<<"("<<mustdecay<<") throws error for: "<<cluster<<"\n"
	       <<(*cluster)<<"\n";
    cluster->Clear();
    delete cluster;
    return false;
  case 1:
    // cluster decayed into hadrons - delete it and carry on.
    cluster->Clear();
    delete cluster;
    return true;
  case 0:
  default:
    //cluster should have decayed into clusters - throw error
    break;
  }
  msg_Error()<<METHOD<<" throws error for:\n"<<(*cluster)<<"\n";
  return false;
}
