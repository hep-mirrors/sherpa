#include "AHADIC++/Decays/Cluster_Decay_Handler.H"
#include "AHADIC++/Decays/Cluster_Part.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_List * clulist,bool ana) :
  p_softclusters(hadpars.GetSoftClusterHandler()),
  p_clus(new Cluster_Part(hadpars.GetSplitter(),ana)),
  p_clulist(clulist),
  p_analysis(ana?new Cluster_Decay_Analysis():NULL)
{ }



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (p_clus)     { delete p_clus;     p_clus=NULL;     }
  if (p_analysis) { delete p_analysis; p_analysis=NULL; }
}

int Cluster_Decay_Handler::DecayClusters(Blob * blob)
{
  SP(Cluster) cluster;
  Cluster_List clist;
  //std::cout<<METHOD<<" for "<<p_clulist->size()<<std::endl;
  while (!p_clulist->empty()) {
    cluster = p_clulist->front();
    if (cluster->Active()) {
      if (!p_clus->TestDecay(cluster)) return -1;
      clist.push_back(cluster->GetLeft());
      clist.push_back(cluster->GetRight());
      if (!p_softclusters->TreatClusterList(&clist,blob)) {
	msg_Tracking()<<"Error in "<<METHOD<<" : "<<std::endl
		      <<"   Did not find a kinematically allowed "
		      <<"solution for the cluster list."<<std::endl
		      <<"   Will trigger retrying the event."<<std::endl;
	return -1;
      }
      while (!clist.empty()) {
	p_clulist->push_back(clist.back());
	clist.pop_back();
      }
      cluster->SetActive(false);
    }
    p_clulist->pop_front();
  }
  if (p_analysis) p_analysis->AnalyseThis(blob);  

  return 1;
}

ATOOLS::Blob * Cluster_Decay_Handler::ClusterDecayBlob(Cluster * cluster,Cluster_List * p_clulist) {
  Blob * decblob(cluster->ConstructDecayBlob());
#ifdef AHAmomcheck
  if (dabs(decblob->CheckMomentumConservation().Abs2())>1.e-12) {
    msg_Out()<<METHOD<<" : Momentum violation at cluster decay blob : "
	     <<decblob->CheckMomentumConservation()<<std::endl
	     <<(*decblob)<<std::endl;
  }
  else {
    msg_Debugging()<<METHOD<<" : Momentum conservation at cluster decay blob : "
		   <<decblob->CheckMomentumConservation().Abs2()<<std::endl;
  }
#endif
  if (cluster->GetLeft()!=NULL && cluster->GetLeft()->GetFlav()==Flavour(kf_cluster)) {
    p_clulist->push_back(cluster->GetLeft());
  }
  if (cluster->GetRight()!=NULL && cluster->GetRight()->GetFlav()==Flavour(kf_cluster)) {
    p_clulist->push_back(cluster->GetRight());
  }
  if (cluster) {
#ifdef memchecker
    std::cout<<"@@@ Delete cluster "<<cluster<<" in "<<METHOD<<"."<<std::endl;
#endif
    cluster->SetActive(false);
  }
  return decblob;
}


