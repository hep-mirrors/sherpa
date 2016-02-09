#include "AHADIC++/Decays/Cluster_Decay_Handler.H"
#include "AHADIC++/Decays/Cluster_Part.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_List * clulist,bool ana) :
  p_softclusters(hadpars->GetSoftClusterHandler()),
  p_clus(new Cluster_Part(ana)),
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
  Cluster * cluster;
  Cluster_Iterator cit(p_clulist->begin());
  //msg_Out()<<":::::: "<<METHOD<<" with "<<p_clulist->size()<<" clusters.\n";
  while (!p_clulist->empty()) {
    cluster = p_clulist->front();
    if (!cluster->Active()) return -1;
    if (p_clus->TestDecay(cluster)) {
      //msg_Out()<<":::::: "<<METHOD<<": decay ok for\n"<<(*cluster);
      Cluster_List * clist(cluster->GetClusters());
      if (!p_softclusters->TreatClusterList(clist,blob)) {
	msg_Error()<<"Error in "<<METHOD<<" : \n"
		   <<"   Did not find a kinematically allowed "
		   <<"solution for the cluster list.\n"
		   <<"   Will trigger retrying the event.\n";
	return -1;
      }
      while (!clist->empty()) {
	p_clulist->push_back(clist->front());
	clist->pop_front();
      }
      //msg_Out()<<":::::: "<<p_clulist->size()<<" clusters now.\n";
    }
    else {
      //msg_Out()<<":::::: "<<METHOD<<":\n"
      //       <<"   Enter soft cluster treatment for undecayed cluster.\n"
      //       <<(*cluster);
      Cluster_List clist;
      clist.push_back(cluster);
      if (!p_softclusters->TreatClusterList(&clist,blob)) {
	msg_Error()<<"Error in "<<METHOD<<".\n";
	return -1;
      }
    }
    delete (p_clulist->front()->GetTrip());
    delete (p_clulist->front()->GetAnti());
    delete (p_clulist->front());
    p_clulist->pop_front();
  }
  if (p_analysis) p_analysis->AnalyseThis(blob);  

  return 1;
}
