#include "AHADIC++/Decays/Cluster_Decayer.H"
#include "AHADIC++/Tools/Hadronisation_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace AHADIC;
using namespace ATOOLS;
using namespace std;

Cluster_Decayer::Cluster_Decayer(list<Cluster *> * cluster_list,
				 Soft_Cluster_Handler * softclusters,
				 Flavour_Selector     * flavourselector) :
  p_cluster_list(cluster_list), p_softclusters(softclusters),
  m_splitter(Cluster_Splitter(cluster_list,softclusters,flavourselector))
{}

Cluster_Decayer::~Cluster_Decayer() {}

void Cluster_Decayer::Init() {
  m_splitter.Init();
}

void Cluster_Decayer::Reset() {}

bool Cluster_Decayer::operator()(bool breakit) {
  while (!p_cluster_list->empty()) {
    if (!Treat(p_cluster_list->front())) {
      return false;
    }
    p_cluster_list->pop_front();
    if (breakit) return true;
  }
  return true;
}

bool Cluster_Decayer::Treat(Cluster * cluster) {
  bool mustdecay = p_softclusters->MustPromptDecay(cluster);
  if (!mustdecay && m_splitter((*cluster)[0],(*cluster)[1])) {
    delete cluster;
    return true;
  } 
  switch (p_softclusters->Treat(cluster,true)) {
  case -1:
    // cluster cannot decay into anything - return false (triggers new event)
    msg_Error()<<METHOD<<"("<<mustdecay<<") throws error for: "<<cluster<<"\n"
	       <<(*cluster)<<"\n";
    cluster->Clear();
    delete cluster;
    return false;
  case 1:
    // cluster decayed into hadrons - delete it and carry on.
    cluster->Clear();
    delete cluster;
    return true;
  case 0:
  default:
    //cluster should have decayed into clusters - throw error
    break;
  }
  return false;
}

void Cluster_Decayer::Test(const double & Q, const Flavour & flav,
			   bool clustermasses) {
  map<string,Histogram *> histos;
  if (clustermasses) {
    histos[string("M_clusters")] = new Histogram(0,0.,Q,int(10.*Q));
    histos[string("x_clusters")] = new Histogram(0,0.,1.,200);
  }
  else {
    histos[string("N_hadrons")]  = new Histogram(0,0.,int(2.*Q),int(2.*Q));
    histos[string("pt_hadrons")] = new Histogram(0,0.,2.,200);
    histos[string("y_hadrons")]  = new Histogram(0,-log(Q/0.1),log(Q/0.1),100);
    histos[string("x_hadrons")]  = new Histogram(0,0.,1.,200);
    histos[string("xb_hadrons")] = new Histogram(0,0.,1.,200);
  }
  double momz    = 0.5*sqrt(Q*Q-sqr(hadpars->GetConstituents()->Mass(flav)));
  for (long int i=0;i<1000000;i++) {
    if (i%10000==0) msg_Out()<<"* "<<int(i/10000)<<" M clusters.\n";
    Proto_Particle * q    = new Proto_Particle(flav,
					       Vec4D(Q/2.,0.,0.,momz));
    Proto_Particle * qbar = new Proto_Particle(flav.Bar(),
					       Vec4D(Q/2.,0.,0.,-momz));
    p_cluster_list->push_back(new Cluster(q,qbar));
    if ((*this)(clustermasses)) {
      if (clustermasses) {
	while (!p_cluster_list->empty()) {
	  Vec4D mom = p_cluster_list->front()->Momentum();
	  histos[string("M_clusters")]->Insert(sqrt(mom.Abs2()));
	  histos[string("x_clusters")]->Insert(mom.PSpat()/momz);
	  delete p_cluster_list->front();
	  p_cluster_list->pop_front();
	}
      }
      else {
	list<Proto_Particle *> * hadrons = p_softclusters->GetHadrons();
	histos[string("N_hadrons")]->Insert(hadrons->size()+0.5);
	while (!hadrons->empty()) {
	  Vec4D mom = hadrons->front()->Momentum();
	  histos[string("pt_hadrons")]->Insert(mom.PPerp());
	  histos[string("y_hadrons")]->Insert(mom.Y());
	  if (hadrons->front()->Flavour().IsB_Hadron() ||
	      hadrons->front()->Flavour().IsC_Hadron())
	    histos[string("xb_hadrons")]->Insert(mom.PSpat()/momz);
	  else
	    histos[string("x_hadrons")]->Insert(mom.PSpat()/momz);
	  delete hadrons->front();
	  hadrons->pop_front();
	}
      }
    }
  }
  Histogram * histo;
  string name;
  for (map<string,Histogram *>::iterator hit=histos.begin();
       hit!=histos.end();hit++) {
    histo = hit->second;
    name  = string("Fragmentation_Analysis/")+hit->first+string(".dat");
    histo->Output(name);
    delete histo;
  }
  histos.clear();
}
