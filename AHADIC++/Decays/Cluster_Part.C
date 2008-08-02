#include "Cluster_Part.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(Dipole_Splitter * splitter,bool ana) :
  m_ana(ana), p_splitter(splitter)
{ 
  if (m_ana) {
    m_histograms[string("Flavour_Cluster")] = new Histogram(0,0.,15.,15);
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

bool Cluster_Part::TestDecay(SP(Cluster) const cluster)
{
  if (!cluster->Active()) return true;
#ifdef AHAmomcheck
  Vec4D checkbef = cluster->Momentum();
#endif
  cluster->BoostInCMSAndRotateOnZ();
  if (!p_splitter->SplitCluster(cluster)) {
    msg_Tracking()<<"ERROR in "<<METHOD<<":"<<std::endl
		  <<"   Could not split cluster, may lead to new event."<<std::endl;
    return false;
  }
  if (m_ana) {
    Vec4D lmom(cluster->GetLeft()->Momentum());
    double pt = sqrt(sqr(lmom[1]) + sqr(lmom[2]));
    Histogram* histo((m_histograms.find(std::string("PT_Cluster")))->second);
    histo->Insert(pt);
  }

  cluster->RotateAndBoostBack();
#ifdef AHAmomcheck
  Vec4D checkaft = cluster->GetLeft()->Momentum()+cluster->GetRight()->Momentum();
  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"    Four-momentum not conserved: "
	       <<checkbef<<" vs. "<<checkaft<<" : "<<(checkbef-checkaft).Abs2()<<"."<<std::endl;
  }
#endif
  return true;
}
