#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadron_Part.H"
#include "Hadronisation_Parameters.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_Transformer * transformer,
					     bool cib,bool ana) :
  m_cib(cib), 
  p_transformer(transformer),p_transitions(hadpars.GetSingleTransitions()),
  p_analysis(NULL), p_clusters(NULL), p_blob(NULL),  
  m_offset1(hadpars.Get(string("Offset_C->H"))),
  m_offset2(hadpars.Get(string("Offset_C->HH")))
{ 
  p_clus = new Cluster_Part();
  p_hads = new Hadron_Part();
  if (ana) p_analysis = new Cluster_Decay_Analysis();
}



Cluster_Decay_Handler::~Cluster_Decay_Handler()
{ 
  if (p_clus)     { delete p_clus;     p_clus=NULL;     }
  if (p_hads)     { delete p_hads;     p_hads=NULL;     }
  if (p_analysis) { delete p_analysis; p_analysis=NULL; }
}

Return_Value::code Cluster_Decay_Handler::DecayClusters(Cluster_List * clusters,
							Blob_List * blobs)
{
  p_clusters   = clusters;
  Blob * blob(NULL);
  Cluster_Iterator cit=clusters->begin();
  while (!clusters->empty()) {
    //cout<<"################################ "<<(*cit)->Number()<<" #####################"<<endl;
    if (DecayIt((*cit),blob)) {
      blobs->push_back(blob);
      if ((*cit)) { delete (*cit); (*cit)=NULL; }
      cit=clusters->erase(cit); 
    }
    else return Return_Value::Error;
  }
  
  if (blob!=NULL && p_analysis) p_analysis->AnalyseThis(blob);  

  return Return_Value::Success;
}


bool Cluster_Decay_Handler::DecayIt(Cluster * cluster,Blob *& blob)
{
  InitDecayBlob(cluster,blob);
  int test(3);
  Flavour had1=Flavour(kf::none),had2=Flavour(kf::none);
  if (p_hads->MustTransit(cluster,had1,had2,m_offset2)) {
    p_hads->FixHHDecay(cluster,had1,had2);
    test = 3;
  }
  else if (p_clus->TestDecay(cluster)) {
    test = TestOffSprings(cluster);
    if (test!=0) p_clus->UpdateDecay(cluster,test);
  }
  else {
    return false;
  }
  FillDecayBlob(cluster,blob,test);
  return true;
}

int Cluster_Decay_Handler::TestOffSprings(Cluster * cluster)
{
  Flavour had;
  int test = int(p_transitions->MustTransit(cluster->GetLeft(),had,m_offset1,false));
  if (test==1) {
    cluster->GetLeft()->GetSelf()->SetFlav(had);
    cluster->GetLeft()->GetSelf()->SetInfo('P');
    cluster->GetLeft()->GetSelf()->SetFinalMass(had.PSMass());
  }
  test    += 2*int(p_transitions->MustTransit(cluster->GetRight(),had,m_offset1,false));
  if (test>=2) {
    cluster->GetRight()->GetSelf()->SetFlav(had);
    cluster->GetRight()->GetSelf()->SetInfo('P');
    cluster->GetRight()->GetSelf()->SetFinalMass(had.PSMass());
  }
  return test;
}

void Cluster_Decay_Handler::InitDecayBlob(Cluster * cluster,Blob *& blob)
{
  blob = new Blob();
  control::s_AHAblobs++;
  blob->SetType(btp::Cluster_Decay);
  blob->SetTypeSpec("AHADIC-1.0");
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetId();
  blob->AddToInParticles(cluster->GetSelf());
  cluster->GetSelf()->SetStatus(part_status::decayed);
  cluster->GetSelf()->ProductionBlob()->UnsetStatus(blob_status::needs_hadrondecays);
}

void Cluster_Decay_Handler::FillDecayBlob(Cluster * cluster,Blob * blob,const int mode)
{
  blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
  blob->AddToOutParticles(cluster->GetRight()->GetSelf());
  switch (mode) {
  case 0:
    p_clusters->push_back(cluster->GetLeft());
    p_clusters->push_back(cluster->GetRight());
    break;
  case 1:
    if (cluster->GetLeft()) { 
      cluster->DeleteLeft(); 
    }
    p_clusters->push_back(cluster->GetRight());
    break;
  case 2:
    if (cluster->GetRight()) { 
      cluster->DeleteRight(); 
    }
    p_clusters->push_back(cluster->GetLeft());
    break;
  case 3:
    if (cluster->GetLeft())  cluster->DeleteLeft();
    if (cluster->GetRight()) cluster->DeleteRight();
    break;
  }
  //cout<<METHOD<<endl<<(*blob)<<endl;
}
