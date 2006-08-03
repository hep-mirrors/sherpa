#include "Cluster_Decayer_Base.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decayer_Base::Cluster_Decayer_Base(Cluster_Part * decs,Hadron_Part * hads) :
  p_stransitions(hadpars.GetSingleTransitions()),
  p_cdecs(decs), p_chads(hads), 
  m_test(0), m_offset(hadpars.Get(string("Offset")))
{ }

Cluster_Decayer_Base::~Cluster_Decayer_Base()
{
  if (p_cdecs) { delete p_cdecs; p_cdecs = NULL; }
  if (p_chads) { delete p_chads; p_chads = NULL; }
}
 
Return_Value::code Cluster_Decayer_Base::Treat(Cluster * cluster,Part_List * pl)
{
  switch (int(p_cdecs->TestDecay(cluster,pl))) {
  case Return_Value::Success :
    //cout<<METHOD<<" 1: Success."<<endl;
    TestOffSprings(cluster);
    if (m_test>0) return TreatHadDecay(cluster,pl); 
    return Return_Value::Nothing;
  case Return_Value::Nothing :
    //cout<<METHOD<<" 1: Nothing."<<endl;
    return p_chads->ForcedDecay(cluster,pl);
  case Return_Value::Error :
    //cout<<METHOD<<" 1: Error."<<endl;
    break;
  default:
    msg.Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Unknown return value."<<endl;
    abort();
    break;
  }
  //cout<<METHOD<<" 2: Error."<<endl;
  return Return_Value::Error;
}

void Cluster_Decayer_Base::TestOffSprings(Cluster * cluster)
{
  m_test  =   int(p_stransitions->MustTransit(cluster->GetLeft(),m_had1,m_offset));
  m_test += 2*int(p_stransitions->MustTransit(cluster->GetRight(),m_had2,m_offset));
}

Return_Value::code Cluster_Decayer_Base::TreatHadDecay(Cluster * cluster,Part_List * pl)
{
  switch (int(p_chads->RedoDecay(cluster,pl,m_test,m_had1,m_had2))) {
  case int(Return_Value::Success) :
    if (m_test&1) cluster->DeleteLeft();
    if (m_test&2) cluster->DeleteRight();
    //cout<<METHOD<<" 1: Success : "<<m_test<<endl;
    return Return_Value::Success;
  case int(Return_Value::Error) :
    //cout<<METHOD<<" 1: Error."<<endl;
    return Return_Value::Error;
  default:
    msg.Error()<<"Error in "<<METHOD<<": "<<endl
	       <<"   Unknown return value."<<endl;
    abort();
  }
  //cout<<METHOD<<" 1: Success."<<endl;
  return Return_Value::Undefined;
}
