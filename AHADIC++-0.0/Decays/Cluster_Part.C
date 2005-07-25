#include "Cluster_Part.H"
#include "Hadronisation_Parameters.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Simple Q over M
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Simple_Q_over_M::Simple_Q_over_M() :
  m_popper(Pair_Popper()), m_Q(hadpars.Get(string("Q_breakup")))
{ }

Simple_Q_over_M::~Simple_Q_over_M() { }



bool Simple_Q_over_M::TestDecay(Cluster * cluster,Part_List * pl)
{
  cluster->BoostInCMS();
  double M    = sqrt(cluster->Momentum(0).Abs2());
  double m1   = cluster->Mass(1),         m2    = cluster->Mass(2);
  Vec4D orig1 = cluster->Momentum(1),     orig2 = cluster->Momentum(2);

  Vec4D mom1  = (1.-m_Q/M)*orig1,         mom2  = m_Q/M*orig2;
  Vec4D mom3  = m_Q/M*orig1,              mom4  = (1.-m_Q/M)*orig2;

  double E1   = sqrt((mom1+mom2).Abs2()), E2    = sqrt((mom3+mom4).Abs2()); 

  Flavour flav;
  bool massflag = true;
  do {
    if (m_popper.Pop(flav)) {
      if (m1+hadpars.GetConstituents()->Mass(flav)<E1 &&
	  m2+hadpars.GetConstituents()->Mass(flav)<E2) {
	massflag = false;
      }
    }
    else { 
      cluster->BoostBack();
      return false;
    }
  } while (massflag);
  
  Cluster * left  = new Cluster(cluster->GetFlav(1),mom1,flav.Bar(),mom2);
  Cluster * right = new Cluster(flav,mom3,cluster->GetFlav(2),mom4);
  if (cluster->GetLeads()==1 || cluster->GetLeads()==3) left->SetLeads(ltp::leadingtrip);
  if (cluster->GetLeads()==2 || cluster->GetLeads()==3) right->SetLeads(ltp::leadinganti);
  
  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  cluster->BoostBack();
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//  Running Q over M
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Running_Q_over_M::Running_Q_over_M() :
  m_popper(Pair_Popper()), m_Q(hadpars.Get(string("Q_breakup")))
{ }

Running_Q_over_M::~Running_Q_over_M() { }


double Running_Q_over_M::SelectQ(const double M)
{
  double Q, norm = M*M;
  do { Q = m_Q+ran.Get()*sqrt(m_Q*M); } while (exp(-(Q*Q)/norm)>ran.Get());
  return Q;
  return (m_Q*M)/(m_Q+M);
}

bool Running_Q_over_M::TestDecay(Cluster * cluster,Part_List * pl)
{
  cluster->BoostInCMS();
  double M    = sqrt(cluster->Momentum(0).Abs2());
  double m1   = cluster->Mass(1),         m2    = cluster->Mass(2);
  Vec4D orig1 = cluster->Momentum(1),     orig2 = cluster->Momentum(2);
  double Q    = SelectQ(M);

  Vec4D mom1  = (1.-Q/M)*orig1,         mom2  = Q/M*orig2;
  Vec4D mom3  = Q/M*orig1,              mom4  = (1.-Q/M)*orig2;

  double E1   = sqrt((mom1+mom2).Abs2()), E2    = sqrt((mom3+mom4).Abs2()); 

  Flavour flav;
  bool massflag = true;
  do {
    if (m_popper.Pop(flav)) {
      if (m1+hadpars.GetConstituents()->Mass(flav)<E1 &&
	  m2+hadpars.GetConstituents()->Mass(flav)<E2) {
	massflag = false;
      }
    }
    else { 
      cluster->BoostBack();
      return false;
    }
  } while (massflag);
  
  Cluster * left  = new Cluster(cluster->GetFlav(1),mom1,flav.Bar(),mom2);
  Cluster * right = new Cluster(flav,mom3,cluster->GetFlav(2),mom4);
  if (cluster->GetLeads()==1 || cluster->GetLeads()==3) left->SetLeads(ltp::leadingtrip);
  if (cluster->GetLeads()==2 || cluster->GetLeads()==3) right->SetLeads(ltp::leadinganti);
  
  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  cluster->BoostBack();
  return true;
}
