#include "AHADIC++/Decays/Cluster_Part.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(Dipole_Splitter * splitter,bool ana) :
  m_ana(ana), p_splitter(splitter), m_fails(0), m_att(0)
{ 
  if (m_ana) {
    m_histograms[string("PT_Cluster")] = new Histogram(0,0.,1.5,150);
  }
}

Cluster_Part::~Cluster_Part()
{
  if (m_ana) {
    msg_Out()<<METHOD<<" (had "<<m_fails<<" fails in "<<m_att<<" attempts).\n";
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
  m_att++;
  if (!p_splitter->SplitCluster(cluster)) {
    m_fails++;
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
  return true;
}
