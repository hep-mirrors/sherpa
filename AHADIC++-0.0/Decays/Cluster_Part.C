#include "Cluster_Part.H"
#include "Message.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;


Q_over_M::Q_over_M() :
  m_popper(Pair_Popper()), m_Q(1.)
{ }

Q_over_M::~Q_over_M() { }



bool Q_over_M::TestDecay(Cluster * cluster,Part_List * pl)
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

  cluster->SetLeft(left);
  cluster->SetRight(right);
  left->SetPrev(cluster);
  right->SetPrev(cluster);

  cluster->BoostBack();
  return true;
}
