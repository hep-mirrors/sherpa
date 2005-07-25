#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadron_Part.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_Transformer * transformer,bool ana) :
  m_cdm(cdm::RunningQoverM_Retain), 
  p_decayer(NULL), p_analysis(NULL), 
  p_transformer(transformer),
  p_partlist(new Part_List)
{ 
  Cluster_Part * cp = NULL;
  Hadron_Part  * hp = NULL;
  switch (int(m_cdm/10)) {
  case 4:
    cp = new Four_Fermion();
    break;
  case 2:
    cp = new Running_Q_over_M();
    break;
  case 1:
  default:
    cp = new Simple_Q_over_M();
  }
  switch (int(m_cdm%10)) {
  case 2:
    hp = new Retain();
    break;
  case 1:
  default:
    hp = new Isotropic();
  }
  p_decayer = new Cluster_Decayer_Base(cp,hp);
  if (ana) p_analysis = new Cluster_Decay_Analysis();
}



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (p_decayer)  { delete p_decayer;  p_decayer=NULL;  }
  if (p_analysis) { delete p_analysis; p_analysis=NULL; }
  if (p_partlist) { delete p_partlist; p_partlist=NULL; }
}

void Cluster_Decay_Handler::DecayClusters(Cluster_List * clusters,Blob * blob)
{
  p_partlist->clear();
  msg.Tracking()<<"Decay the clusters ------------------------------------------------"<<endl;
  //cout<<"#################################################################################"<<endl<<(*clusters)<<endl;
  Cluster_Iterator cit;
  Vec4D clumom = Vec4D(0.,0.,0.,0.), partmom = Vec4D(0.,0.,0.,0.);
  for (cit=clusters->begin();cit!=clusters->end();) {
    clumom += (*cit)->Momentum();
    if (DecayIt((*cit))) cit++;
    else cit=clusters->erase(cit);
  }
  msg.Tracking()<<"Add "<<p_partlist->size()
		<<" particles to the blob --------------------------------------"<<endl;
  for (Part_Iterator pit=p_partlist->begin();pit!=p_partlist->end();) {
    blob->AddToOutParticles((*pit));
    partmom += (*pit)->Momentum();
    pit = p_partlist->erase(pit);
  }

  if (dabs(blob->CheckMomentumConservation().Abs2())>1.e-9) {
    msg.Tracking()<<"Check this : "
		  <<blob->CheckMomentumConservation()<<", "<<blob->CheckMomentumConservation().Abs2()<<endl
		  <<"   Compare with "<<clumom<<" -> "<<partmom<<" = "<<clumom-partmom<<endl
		  <<(*blob)
		  <<"----------------------------------------------------------"<<endl;
  }
  if (p_analysis) p_analysis->AnalyseThis(blob);
}




bool Cluster_Decay_Handler::DecayIt(Cluster * cluster)
{
  if (p_decayer->Treat(cluster,p_partlist)) {
    //cout<<"Decay "<<endl<<(*cluster)<<endl;
    if (cluster->GetLeft())  DecayIt(cluster->GetLeft());
    if (cluster->GetRight()) DecayIt(cluster->GetRight());
    return true;
  }
  //cout<<"Treat "<<endl<<(*cluster)<<endl;
  p_transformer->TreatSingleCluster(cluster,p_partlist);
  if (cluster->GetPrev()!=NULL) cout<<"------------- Found prev."<<endl;
  return false;
}
