#include "Cluster_Former.H"
#include "Message.H"
#include "Random.H"

using namespace AHADIC;
using namespace ATOOLS;


Cluster_Former::Cluster_Former() : 
  m_kinmode(1), m_colmode(1), m_t0(1.), p_pl(NULL), p_cl(new Cluster_List) 
{ }

Cluster_Former::~Cluster_Former() 
{
  p_cl->clear(); delete p_cl;
}

void Cluster_Former::FormClusters(Part_List * pl)
{
  Reset();
  p_pl = pl;
  ConstructPreClusters();
  ReshuffleClusters();
  FormFinalClusters();
}


void Cluster_Former::Reset() 
{
  for (int i=0;i<m_preclusters.size();i++) {
    if (m_preclusters[i]!=NULL) { delete m_preclusters[i]; m_preclusters[i]=NULL; }
  }
  m_preclusters.clear();
  p_cl->clear();
}

void Cluster_Former::ConstructPreClusters()
{
  int col1;
  precluster * prec;
  Particle   * part;

  for (Part_Iterator pit1=p_pl->begin();pit1!=p_pl->end();++pit1) {
    col1 = (*pit1)->GetFlow(1);
    if (col1!=0) {
      if ((*pit1)->GetFlow(2)!=0) {
	msg.Error()<<"ERROR in Cluster_Former::FormClusters : "<<std::endl
		   <<"   Colour octet left in particle list."<<std::endl
		   <<"   "<<(*pit1)->Number()<<" : "<<(*pit1)->Flav()
		   <<" ("<<(*pit1)->GetFlow(1)<<","<<(*pit1)->GetFlow(2)
		   <<"), abort the run."<<std::endl;
	abort();
      }
      for (Part_Iterator pit2=p_pl->begin();pit2!=p_pl->end();++pit2) {
	if ((*pit2)->GetFlow(2)==col1) {
	  if ((*pit2)->GetFlow(1)!=0) {
	    msg.Error()<<"ERROR in Cluster_Former::FormClusters : "<<std::endl
		       <<"   Colour octet left in particle list."<<std::endl
		       <<"   "<<(*pit2)->Number()<<" : "<<(*pit2)->Flav()
		       <<" ("<<(*pit2)->GetFlow(1)<<","<<(*pit2)->GetFlow(2)
		       <<"), abort the run."<<std::endl;
	    abort();
	  }
	  prec         = new precluster;
	  prec->first  = (*pit1);
	  prec->second = (*pit2);
	  m_preclusters.push_back(prec);
	  (*pit1)->SetStatus(2);
	  (*pit2)->SetStatus(2);
	}
      }
    }
  }
}

void Cluster_Former::ReshuffleClusters()
{
  if (m_preclusters.size()<2) return;
  Particle * help;
  int        col;
  double kinweight, colweight=ColourWeight();
  for (int i=0;i<m_preclusters.size()-1;i++) {
    for (int j=i+1;j<m_preclusters.size();j++) {
      kinweight = KinematicWeight(m_preclusters[i]->first->Momentum(),
				  m_preclusters[i]->second->Momentum(),
				  m_preclusters[j]->first->Momentum(),
				  m_preclusters[j]->second->Momentum());
      if (kinweight*colweight>ran.Get()) {
	help                     = m_preclusters[i]->second;
	m_preclusters[i]->second = m_preclusters[j]->second;
	m_preclusters[j]->second = help;
	m_preclusters[i]->second->SetFlow(2,m_preclusters[i]->first->GetFlow(1));
	m_preclusters[j]->second->SetFlow(2,m_preclusters[j]->first->GetFlow(1));
      }
    }
  }
}

void Cluster_Former::FormFinalClusters()
{
  for (int i=0;i<m_preclusters.size();i++) {
    p_cl->push_back(new Cluster(m_preclusters[i]->first,m_preclusters[i]->second));
  }
}


double Cluster_Former::KinematicWeight(const Vec4D & mom1,const Vec4D & mom2,
				       const Vec4D & mom3,const Vec4D & mom4)
{
  double w12, w34, w14, w23;
  switch (m_kinmode) {
  case 2:
  case 1:
    w12 = sqrt((mom1+mom2).Abs2());
    w34 = sqrt((mom3+mom4).Abs2());
    w14 = sqrt((mom1+mom4).Abs2());
    w23 = sqrt((mom2+mom3).Abs2());
    break;
  default:
    w12 = w34 = 0;
    w14 = w23 = 1.e64;
    break;
  }
  double w1234 = m_t0/(m_t0+4.*sqr(w12+w34));
  double w1423 = m_t0/(m_t0+4.*sqr(w14+w23));

  return w1423/(w1234+w1423);
}


double Cluster_Former::ColourWeight() { 
  switch (m_colmode) {
  case 1:  return 1./9.; 
  default: return 0.;
  }
}
