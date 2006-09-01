#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadronisation_Parameters.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Soft_Cluster_Handler * softclusters,
					     bool cib,bool ana) :
  m_cib(cib), 
  p_softclusters(softclusters),
  p_analysis(NULL)
{ 
  p_clus = new Cluster_Part();
  if (ana) p_analysis = new Cluster_Decay_Analysis();
}



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (p_clus)     { delete p_clus;     p_clus=NULL;     }
  if (p_analysis) { delete p_analysis; p_analysis=NULL; }
}

Return_Value::code Cluster_Decay_Handler::DecayClusters(Cluster_List * clusters,
							Blob_List * blobs)
{
  Cluster    * cluster;
  Blob       * blob;
  Cluster_List clist;
  Cluster_Iterator cit=clusters->begin();
  while (!clusters->empty()) {
    cluster = (*cit);
    //cout<<METHOD<<" ############################### "<<clusters->size()<<endl<<(**cit)<<endl;
    blob    = DecayIt(cluster);
    blobs->push_back(blob);
    clist.push_back(cluster->GetLeft());
    clist.push_back(cluster->GetRight());
    p_softclusters->TreatClusterList(&clist,blob);
    Cluster_Iterator dcit=clist.begin();
    //cout<<"   ---- before loop "<<clist.size()<<endl;
    while (!clist.empty()) {
      if ((*dcit)) {
	//cout<<"  add to blob . "<<(*dcit)->GetLeft()<<" && "<<(*dcit)->GetRight()<<endl;
	blob->AddToOutParticles((*dcit)->GetSelf());
	if (!(*dcit)->GetLeft() && !(*dcit)->GetRight()) {
	  clusters->push_back((*dcit));
	}
	else {
	  //cout<<"                     before fill blob ."<<endl;
	  blobs->push_back((*dcit)->CHHDecayBlob());
	  //cout<<"                      after fill blob ."<<endl;
	  if ((*dcit)->GetLeft()) {
	    if ((*dcit)->GetLeft()->GetSelf()->Flav()==Flavour(kf::none) ||
		(*dcit)->GetLeft()->GetSelf()->Flav()==Flavour(kf::cluster)) {
	      clusters->push_back((*dcit)->GetLeft());
	    }
	  }
	  if ((*dcit)->GetRight()) {
	    if ((*dcit)->GetRight()->GetSelf()->Flav()==Flavour(kf::none) ||
		(*dcit)->GetRight()->GetSelf()->Flav()==Flavour(kf::cluster)) {
	      clusters->push_back((*dcit)->GetRight());
	    }
	  }
	  if (*dcit) delete (*dcit);
	}
      }
      dcit = clist.erase(dcit);
    }
    if (cluster) delete cluster;
    cit = clusters->erase(cit);
  }
  if (blob!=NULL && p_analysis) p_analysis->AnalyseThis(blob);  

  return Return_Value::Success;
}


Blob * Cluster_Decay_Handler::DecayIt(Cluster * cluster)
{
  Blob * blob = new Blob();
  control::s_AHAblobs++;
  blob->SetType(btp::Cluster_Decay);
  blob->SetTypeSpec("AHADIC-1.0");
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetId();
  blob->AddToInParticles(cluster->GetSelf());
  cluster->GetSelf()->SetStatus(part_status::decayed);
  cluster->GetSelf()->ProductionBlob()->UnsetStatus(blob_status::needs_hadrondecays);

  p_clus->TestDecay(cluster);

  return blob;
}

