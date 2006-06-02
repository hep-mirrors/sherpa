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
 
bool Cluster_Decayer_Base::Treat(Cluster * cluster,Part_List * pl)
{
  //cout<<"   Produce a test decay -------------------------------------------"<<endl
  //   <<"   "<<cluster->Mass()<<","<<cluster->GetFlav(1)<<" "<<cluster->GetFlav(2)<<endl;
  if (p_cdecs->TestDecay(cluster,pl)) {
    //cout<<"   Test the offsprings "<<cluster<<endl
    //	<<(*cluster)<<" ----------------------------------"<<endl<<endl;
    TestOffSprings(cluster);
    if (m_test>0) {
      //cout<<"   Treat the hadronic decay "
      //<<m_test<<" ----------------------------------"<<endl;
      TreatHadDecay(cluster,pl); 
    }
    return true;
  }
  if (p_chads->ForcedDecay(cluster,pl)) {
    //cout<<"Forced hadronic decay ----------------------------------"<<endl;
    return true;
  }
  return false;
}

void Cluster_Decayer_Base::TestOffSprings(Cluster * cluster)
{
  //  cout<<"Test decays of cluster : "<<cluster->Mass()
  //    <<" ("<<cluster->GetFlav(1)<<", "<<cluster->GetFlav(2)<<")"<<endl;
  m_test  =   int(p_stransitions->MustTransit(cluster->GetLeft(),m_had1,m_offset));
  m_test += 2*int(p_stransitions->MustTransit(cluster->GetRight(),m_had2,m_offset));


  //   if (m_test>0) {
  //     cout<<"+++ Must decay("<<m_test<<") : Masses : "<<cluster->Mass()<<" -> "
  // 	<<cluster->GetLeft()->Mass()<<" + "<<cluster->GetRight()->Mass()<<"  -> "
  // 	<<m_had1.Mass()<<"("<<m_had1<<") + "<<m_had2.Mass()<<"("<<m_had2<<") = "
  // 	<<m_had1.Mass()+m_had2.Mass()<<endl
  // 	<<"                                            "
  // 	<<"{"<<cluster->GetLeft()->GetFlav(1)<<","<<cluster->GetLeft()->GetFlav(2)<<"}    "
  // 	<<"{"<<cluster->GetRight()->GetFlav(1)<<","<<cluster->GetRight()->GetFlav(2)<<"}"<<endl;
  //   }
}

void Cluster_Decayer_Base::TreatHadDecay(Cluster * cluster,Part_List * pl)
{
  //   cout<<"Check Treat1(m_test = "<<m_test<<"): "<<m_had1<<"/"<<m_had2<<" : "
  //       <<cluster->GetLeft()<<"/"<<cluster->GetRight()<<endl;
  p_chads->RedoDecay(cluster,pl,m_test,m_had1,m_had2);
  if (m_test&1) cluster->DeleteLeft();
  if (m_test&2) cluster->DeleteRight();
  //   cout<<"Check Treat2(m_test = "<<m_test<<") "
  //       <<cluster->GetLeft()<<"/"<<cluster->GetRight()<<endl;
}
