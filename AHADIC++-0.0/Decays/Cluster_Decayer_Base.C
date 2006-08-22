#include "Cluster_Decayer_Base.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"



using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decayer_Base::Cluster_Decayer_Base(Cluster_Part * decs,Hadron_Part * hads,
					   Cluster_Transformer * transformer) :
  p_stransitions(hadpars.GetSingleTransitions()),
  p_dtransitions(hadpars.GetDoubleTransitions()),
  p_cdecs(decs), p_chads(hads), p_transformer(transformer),
  m_offset(hadpars.Get(string("Offset")))
{ }

Cluster_Decayer_Base::~Cluster_Decayer_Base()
{
  if (p_cdecs) { delete p_cdecs; p_cdecs = NULL; }
  if (p_chads) { delete p_chads; p_chads = NULL; }
}
 
Return_Value::code Cluster_Decayer_Base::Treat(Cluster * cluster,Blob *& blob)
{
  Flavour had1=Flavour(kf::none),had2=Flavour(kf::none);
  switch (int(p_cdecs->TestDecay(cluster))) {
  case Return_Value::Success :
    cout<<METHOD<<" TestDecay yields Success, momenta : "
	<<cluster->GetLeft()->Momentum()<<" "<<cluster->GetRight()->Momentum()<<endl;
    InitDecayBlob(cluster,blob);
    int test(TestOffSprings(cluster,had1,had2));
    if (test==0) {
      blob->AddToOutParticles(cluster->GetLeft()->GetSelf());
      blob->AddToOutParticles(cluster->GetRight()->GetSelf());
      return Return_Value::Success;
    }
    else {
      return TreatHadDecay(cluster,blob,test,had1,had2); 
    }
  case Return_Value::Nothing :
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

int Cluster_Decayer_Base::TestOffSprings(Cluster * cluster,Flavour & had1,Flavour & had2)
{
  if (p_dtransitions->MustTransit(cluster,had1,had2,m_offset)) return 3;
  int test = int(p_stransitions->MustTransit(cluster->GetLeft(),had1,m_offset));
  test    += 2*int(p_stransitions->MustTransit(cluster->GetRight(),had2,m_offset));
  if (test==3) {
    p_dtransitions->MustTransit(cluster,had1,had2,m_offset);
  }
  cout<<METHOD<<" test = "<<test<<", select new hadrons : "<<had1<<" & "<<had2<<endl;
  return test;
}

void Cluster_Decayer_Base::InitDecayBlob(Cluster * cluster,Blob *& blob)
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

Return_Value::code Cluster_Decayer_Base::TreatHadDecay(Cluster * cluster,Blob *& blob,int mode,
						       Flavour & had1,Flavour & had2)
{
  switch (int(p_chads->RedoDecay(cluster,blob,mode,had1,had2))) {
  case int(Return_Value::Success) :
    if (mode&1) { delete cluster->GetLeft()->GetSelf();  cluster->DeleteLeft();  }
    if (mode&2) { delete cluster->GetRight()->GetSelf(); cluster->DeleteRight(); }
    return Return_Value::Success;
  case int(Return_Value::Error) :
    return Return_Value::Error;
  default:
    msg.Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Unknown return value."<<endl;
    abort();
  }
  return Return_Value::Undefined;
}
