#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadronisation_Parameters.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(bool ana) :
  p_softclusters(hadpars.GetSoftClusterHandler()),
  p_clus(new Cluster_Part(hadpars.GetSplitter())),
  p_analysis(ana?new Cluster_Decay_Analysis():NULL)
{ }



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (p_clus)     { delete p_clus;     p_clus=NULL;     }
  if (p_analysis) { delete p_analysis; p_analysis=NULL; }
}

int Cluster_Decay_Handler::DecayClusters(Cluster_List * clusters,Blob_List * blobs)
{
  Cluster    * cluster;
  Blob       * blob(NULL);
  Cluster_List clist;
  Cluster_Iterator cit=clusters->begin();
  while (!clusters->empty()) {
    cluster = (*cit);
    blob    = DecayIt(cluster);
    blobs->push_back(blob);
    clist.push_back(cluster->GetLeft());
    clist.push_back(cluster->GetRight());
    p_softclusters->TreatClusterList(&clist,blob);
    Cluster_Iterator dcit=clist.begin();
    while (!clist.empty()) {
      if ((*dcit)) {
	blob->AddToOutParticles((*dcit)->GetSelf());
	if (!(*dcit)->GetLeft() && !(*dcit)->GetRight()) {
	  clusters->push_back((*dcit));
	}
	else {
	  Blob * decblob((*dcit)->ConstructDecayBlob());
#ifdef AHAmomcheck
	  if (dabs(decblob->CheckMomentumConservation().Abs2())>1.e-12) {
	    msg_Out()<<METHOD<<" : Momentum violation at cluster decay blob : "
		     <<decblob->CheckMomentumConservation()<<std::endl
		     <<(*decblob)<<std::endl;
	  }
	  else {
	    msg_Out()<<METHOD<<" : Momentum conservation at cluster decay blob : "
		     <<decblob->CheckMomentumConservation().Abs2()<<std::endl;
	  }
#endif
	  blobs->push_back(decblob);
	  if ((*dcit)->GetLeft()) {
	    if ((*dcit)->GetLeft()->GetSelf()->Flav()==Flavour(kf_none) ||
		(*dcit)->GetLeft()->GetSelf()->Flav()==Flavour(kf_cluster)) {
	      clusters->push_back((*dcit)->GetLeft());
	    }
	  }
	  if ((*dcit)->GetRight()) {
	    if ((*dcit)->GetRight()->GetSelf()->Flav()==Flavour(kf_none) ||
		(*dcit)->GetRight()->GetSelf()->Flav()==Flavour(kf_cluster)) {
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

  return 1;
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

