#include "Cluster_Part.H"
#include "Message.H"
#include "Random.H"


using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Part::Cluster_Part(Dipole_Splitter * splitter) :
  m_ana(true), p_splitter(splitter)
{ 
  if (m_ana) {
    m_histograms[string("YStar_by_YStarMax")] = new Histogram(0,-1.,1.,100);
    m_histograms[string("YStar")]             = new Histogram(0,-5.,5.,100);
    m_histograms[string("YBar_by_YBarMax")]   = new Histogram(0,-1.,1.,100);
    m_histograms[string("YBar")]              = new Histogram(0,-5.,5.,100);
    m_histograms[string("Flavour")]           = new Histogram(0,0.,15.,15);
    m_histograms[string("PT")]                = new Histogram(0,0.,1.5,150);
    m_histograms[string("SQQ")]               = new Histogram(0,0.,5.,100);
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
      name  = string("Analysis/")+hit->first+string(".dat");
      histo->Output(name);
      delete histo;
    }
    m_histograms.clear();
  }
}

bool Cluster_Part::TestDecay(Cluster * const cluster)
{
  Vec4D checkbef = cluster->Momentum();
  cluster->BoostInCMSAndRotateOnZ();
  if (!ClusterDecay(cluster)) return false;
  cluster->RotateAndBoostBack();
  Vec4D checkaft = cluster->GetLeft()->Momentum()+cluster->GetRight()->Momentum();

  if (dabs((checkbef-checkaft).Abs2())>1.e-12) {
    msg_Error()<<"Error in "<<METHOD<<" : "<<std::endl
	       <<"    Four-momentum not conserved: "
	       <<checkbef<<" vs. "<<checkaft<<" : "<<(checkbef-checkaft).Abs2()<<"."<<std::endl;
  }
  Particle * part = new Particle(-1,Flavour(kf_cluster),cluster->GetLeft()->Momentum()); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass(cluster->GetLeft()->Mass());
  control::s_AHAparticles++;
  cluster->GetLeft()->SetSelf(part);
  
  part = new Particle(-1,Flavour(kf_cluster),cluster->GetRight()->Momentum()); 
  part->SetNumber();
  part->SetStatus(part_status::active);
  part->SetInfo('C');
  part->SetFinalMass(cluster->GetRight()->Mass());
  control::s_AHAparticles++;
  cluster->GetRight()->SetSelf(part);

  return true;
}


bool Cluster_Part::ClusterDecay(Cluster * const cluster)
{
  return p_splitter->SplitCluster(cluster);
}



