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
  Cluster * cluster;
  Cluster_List clist;
  msg_Tracking()<<"+++++++++++++++++++++++++++++++++++++++++++++++++++++"<<std::endl
		<<"+++ "<<METHOD<<" for "<<p_clulist->size()<<" clusters."<<std::endl;
  Cluster_Iterator cit(p_clulist->begin());
  while (!p_clulist->empty()) {
    cluster = p_clulist->front();
    if (!cluster->Active()) {
      msg_Error()<<"Error in "<<METHOD<<" inactive cluster:"<<std::endl;
      cluster->Print();
      msg_Error()<<"   Will continue and hope for the best."<<std::endl;
      return -1;
    }
    msg_Tracking()<<"+++ Test cluster ["
		  <<cluster->GetTrip()->m_flav<<"("<<cluster->GetTrip()->m_info<<"), "
		  <<cluster->GetAnti()->m_flav<<"("<<cluster->GetAnti()->m_info<<")]."<<std::endl;
    if (p_clus->TestDecay(cluster)) {
      clist.empty();
      clist.push_back(cluster->GetLeft());
      clist.push_back(cluster->GetRight());
      msg_Tracking()<<"++++ From "<<cluster->Number()<<": "
		    <<cluster->GetLeft()->Number()<<" ("<<cluster->GetLeft()->Mass()<<") "
		    <<cluster->GetRight()->Number()<<" ("<<cluster->GetRight()->Mass()<<"), "
		    <<"popped "<<cluster->GetLeft()->GetAnti()->m_flav<<"."<<std::endl;
      if (!p_softclusters->TreatClusterDecay(&clist,blob)) {
	msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
		   <<"   Did not find a kinematically allowed "
		   <<"solution for the cluster list."<<std::endl
		   <<"   Will trigger retrying the event."<<std::endl;
	return -1;
      }
      while (!clist.empty()) {
	p_clulist->push_back(clist.front());
	clist.pop_front();
      }
    }
    else {
      if (!p_softclusters->EnforcedDecay(cluster,blob,true)) {
	msg_Tracking()<<"  no decay."<<std::endl
		      <<"++++ TestDecay did not work out, will return -1."<<std::endl;
	return -1;
      }
    }
    delete (p_clulist->front()->GetTrip());
    delete (p_clulist->front()->GetAnti());
    delete (p_clulist->front());
    p_clulist->pop_front();
    msg_Tracking()<<"++++ After popping front cluster list has now size = "
		  <<p_clulist->size()<<"."<<std::endl;
  }
  if (p_analysis) p_analysis->AnalyseThis(blob);  

  return 1;
}

ATOOLS::Blob * Cluster_Decay_Handler::ClusterDecayBlob(Cluster * cluster,Cluster_List * p_clulist) {
  Blob * decblob(cluster->ConstructDecayBlob());
  if (cluster->GetLeft()!=NULL && cluster->GetLeft()->GetFlav()==Flavour(kf_cluster)) {
    p_clulist->push_back(cluster->GetLeft());
  }
  if (cluster->GetRight()!=NULL && cluster->GetRight()->GetFlav()==Flavour(kf_cluster)) {
    p_clulist->push_back(cluster->GetRight());
  }
  if (cluster) cluster->SetActive(false);
  return decblob;
}


