#include "Cluster_Decay_Handler.H"
#include "Cluster_Part.H"
#include "Hadron_Part.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decay_Handler::Cluster_Decay_Handler(Cluster_Transformer * transformer,
					     bool cib,bool ana) :
  m_cib(cib), m_cdm(cdm::SchwingerUni_Retain),
  p_stransitions(hadpars.GetSingleTransitions()),
  p_dtransitions(hadpars.GetDoubleTransitions()),
  p_transformer(transformer),p_clus(NULL),p_hads(NULL),
  p_analysis(NULL), p_clusters(NULL), p_blob(NULL),  
  m_offset(hadpars.Get(string("Offset")))
{ 
  switch (int(m_cdm/10)) {
  case 8:
    p_clus = new Schwinger_Uniform();
    break;
  case 4:
    p_clus = new Four_Fermion();
    break;
  case 2:
    p_clus = new Running_Q_over_M();
    break;
  case 1:
  default:
    p_clus = new Simple_Q_over_M();
  }
  switch (int(m_cdm%10)) {
  case 4:
    p_hads = new Keep_PPerpY();
    break;
  case 2:
    p_hads = new Retain();
    break;
  case 1:
  default:
    p_hads = new Isotropic();
  }
  if (ana) {
    p_analysis = new Cluster_Decay_Analysis();
  }
  p_clus->SetDecayAnalysisOn();
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
  cout<<METHOD<<endl
      <<"#######################################################"<<endl
      <<(*clusters)<<endl<<"#######################################################"<<endl;
  p_clusters   = clusters;
  Vec4D clumom = Vec4D(0.,0.,0.,0.), partmom = Vec4D(0.,0.,0.,0.);
  Blob * blob(NULL);
  Cluster_Iterator cit=clusters->begin();
  while (!clusters->empty()) {
    clumom += (*cit)->Momentum();
    cout<<"#######################################################"<<endl
	<<METHOD<<"    Before DecayIt, list length = "<<clusters->size()<<endl;
    cout<<"   Try to decay "<<(*cit)<<endl<<(**cit)<<endl;
    switch (int(DecayIt((*cit),blob))) {
      case int(Return_Value::Success) : 
        blobs->push_back(blob);
        if ((*cit)) { delete (*cit); (*cit)=NULL; }
        cit=clusters->erase(cit); 
	cout<<METHOD<<" : DecayIt yields success, list length = "<<clusters->size()<<endl;
        break;
      case int(Return_Value::Error) :
        return Return_Value::Retry_Method; 
      default:
        msg.Error()<<"Error in "<<METHOD<<": "<<endl
		   <<"   Unknown return value."<<endl;
        abort();
    }
  }
  
  if (blob!=NULL) {
    if (p_analysis) p_analysis->AnalyseThis(blob);
  
    if (dabs(blob->CheckMomentumConservation().Abs2())>1.e-6) {
      msg.Tracking()<<"Check this : "
		    <<blob->CheckMomentumConservation()<<", "
		    <<blob->CheckMomentumConservation().Abs2()<<endl
		    <<"   Compare with "<<clumom<<" -> "<<partmom<<" = "<<clumom-partmom<<endl
		    <<(*blob)
		    <<"----------------------------------------------------------"<<endl;
      rvalue.IncWarning(METHOD);
      return Return_Value::Warning;
    }
  }
  msg.Tracking()<<METHOD<<": Success"<<endl;
  return Return_Value::Success;
}

Return_Value::code Cluster_Decay_Handler::DecayIt(Cluster * cluster,Blob *& blob)
{
  Flavour had1=Flavour(kf::none),had2=Flavour(kf::none);

  switch (int(p_clus->TestDecay(cluster))) {
  case Return_Value::Success :
    cout<<METHOD<<" TestDecay yields Success,"<<endl
	<<"         momenta : "
	<<cluster->GetLeft()->Momentum()<<" ("<<cluster->GetLeft()->Momentum().Abs2()<<"), "
	<<cluster->GetRight()->Momentum()<<" ("<<cluster->GetRight()->Momentum().Abs2()<<")."<<endl;
    InitDecayBlob(cluster,blob);
    int test(TestOffSprings(cluster,had1,had2));
    if (TreatHadDecay(cluster,blob,test,had1,had2)) {
      FillDecayBlob(cluster,blob,test);
      return Return_Value::Success;
    }
    break;
  case Return_Value::Nothing :
    cout<<METHOD<<" TestDecay yields Nothing."<<endl;
    InitDecayBlob(cluster,blob);
    if (p_transformer->TreatSingleCluster(cluster,blob)==Return_Value::Success) {
      return Return_Value::Success;
    }
   break;
  case Return_Value::Error :
    break;
  default:
    msg.Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Unknown return value."<<endl;
    abort();
    break;
  }
  return Return_Value::Error;
}

bool Cluster_Decay_Handler::TreatHadDecay(Cluster * cluster,Blob *& blob,int & mode,
					  Flavour & had1,Flavour & had2)
{
  if (mode==0) return true;
  cout<<METHOD<<" mode = "<<mode<<endl;
  switch (int(p_hads->RedoDecay(cluster,blob,mode,had1,had2))) {
  case int(Return_Value::Success) :
    cout<<"Again in "<<METHOD<<", mode = "<<mode<<endl;
    return true;
  case int(Return_Value::Error) : break;
    return false;
  default:
    msg.Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Unknown return value."<<endl;
    abort();
  }
  return false;
}

int Cluster_Decay_Handler::TestOffSprings(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  cout<<METHOD<<", Masses = "<<cluster->Mass(0)<<"  --> "
      <<cluster->GetLeft()->Mass(0)<<" + "<<cluster->GetRight()->Mass(0)<<endl;
  if (p_dtransitions->MustTransit(cluster,had1,had2,m_offset)) return 3;
  int test =   int(p_stransitions->MustTransit(cluster->GetLeft(),had1,m_offset));
  test    += 2*int(p_stransitions->MustTransit(cluster->GetRight(),had2,m_offset));
  if (test==3) {
    p_dtransitions->MustTransit(cluster,had1,had2,m_offset);
  }
  cout<<"         yields "<<test<<", select new hadrons : "<<had1<<" & "<<had2<<endl;
  return test;
}

void Cluster_Decay_Handler::InitDecayBlob(Cluster * cluster,Blob *& blob)
{
  blob = new Blob();
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
  cout<<METHOD<<" mode = "<<mode<<endl<<(*cluster)<<endl;
  switch (mode) {
  case 0:
    p_clusters->push_back(cluster->GetLeft());
    p_clusters->push_back(cluster->GetRight());
    blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
    blob->AddToOutParticles(cluster->GetRight()->GetSelf());
    break;
  case 1:
    if (cluster->GetLeft()) {
      if (cluster->GetLeft()->GetSelf())  delete cluster->GetLeft()->GetSelf();  
      cluster->DeleteLeft();
    } 
    blob->AddToOutParticles(cluster->GetRight()->GetSelf());
    p_clusters->push_back(cluster->GetRight());
    break;
  case 2:
    if (cluster->GetRight()) {
      if (cluster->GetRight()->GetSelf()) delete cluster->GetRight()->GetSelf(); 
      cluster->DeleteRight();
    }
    blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
    p_clusters->push_back(cluster->GetLeft());
    break;
  case 3:
    cluster->DeleteLeft(); 
    cluster->DeleteRight();
    break;
  }
  cout<<METHOD<<" : "<<mode<<", added to outparticles : "<<endl<<(*blob)<<endl;
}
