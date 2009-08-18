#include "AHADIC++/Decays/Cluster_Part.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(Dipole_Splitter * splitter,bool ana) :
  m_ana(ana),
  p_splitter(splitter) 
{ 
  if (m_ana) {
    m_histograms[string("PT_Cluster")]      = new Histogram(0,0.,1.5,150);
  }
}

Cluster_Part::~Cluster_Part()
{
  if (m_ana) {
    Histogram * histo;
    string name;
    for (map<string,Histogram *>::iterator hit=m_histograms.begin();
	 hit!=m_histograms.end();hit++) {
      histo = hit->second;
      name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}

bool Cluster_Part::TestDecay(Cluster * const cluster)
{
  msg_Tracking()<<":::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl
		<<"::: "<<METHOD<<" : Try "<<cluster->Number()<<" ("	   
		<<cluster->GetTrip()->m_flav<<" "<<cluster->GetAnti()->m_flav<<", "
		<<"mass = "<<cluster->Mass()<<") --> "<<std::endl;
  if (!p_splitter->SplitCluster(cluster)) {
    msg_Tracking()<<"::: Warning in "<<METHOD<<":"<<std::endl
		  <<":::   Could not split cluster ("<<cluster->Number()<<"): "
		  <<cluster->GetTrip()->m_flav<<"/"<<cluster->GetAnti()->m_flav<<", "
		  <<"mass = "<<cluster->Mass()<<","<<std::endl
		  <<":::   try to enforce cluster splitting into hadrons."<<std::endl;
    return false;
  }
  if (m_ana) {
    Vec4D lmom(cluster->GetLeft()->Momentum());
    double pt = sqrt(sqr(lmom[1]) + sqr(lmom[2]));
    Histogram* histo((m_histograms.find(std::string("PT_Cluster")))->second);
    histo->Insert(pt);
  }
#ifdef AHAmomcheck
  cluster->CheckConsistency(msg_Error(),METHOD);
#endif
  msg_Tracking()<<"::: "<<METHOD<<" : decay of "<<cluster->Momentum()<<" succeded."<<std::endl
		<<(*cluster->GetLeft())<<(*cluster->GetRight())
		<<":::::::::::::::::::::::::::::::::::::::::::::::"<<std::endl;	   
  return true;
}
