#include "Jet_Y.H"

#include "Scaling.H"

using namespace AMISIC;

Jet_Y::Jet_Y(int _m_type,double _m_min,double _m_max,int _m_nbins,int _m_njets):
  Primitive_Observable_Base(_m_type,_m_min,_m_max,_m_nbins,NULL),
  p_histogram(new std::vector<ATOOLS::Histogram*>(_m_njets)),
  p_jets(new std::vector<int>(_m_njets)),
  p_ys(new std::vector<double>()),
  m_njets(_m_njets),
  m_mode(0)
{
  histo = new ATOOLS::Histogram(type,xmin,xmax,nbins);
  for (unsigned int i=0;i<m_njets;++i) {
    (*p_histogram)[i] = new ATOOLS::Histogram(type,xmin,xmax,nbins);
    (*p_jets)[i]=i;
  }
  name=std::string("Y_jet");
}

Jet_Y::~Jet_Y()
{
  while (p_histogram->size()>0) {
    delete p_histogram->front();
    p_histogram->erase(p_histogram->begin());
  }
  delete p_histogram;
  delete p_jets;
  delete p_ys;
}

void Jet_Y::SortJetPT(std::vector<ATOOLS::Vec4D> &jetmomenta)
{
  ATOOLS::Vec4D help;
  for (unsigned int i=jetmomenta.size()-1;i>0;--i) {
    for (unsigned int j=i;j<jetmomenta.size();++j) {
      if ((jetmomenta[j-1][1]*jetmomenta[j-1][1]+jetmomenta[j-1][2]*jetmomenta[j-1][2])<
	  (jetmomenta[j][1]*jetmomenta[j][1]+jetmomenta[j][2]*jetmomenta[j][2])) {
	help=jetmomenta[j-1];
	jetmomenta[j-1]=jetmomenta[j];
	jetmomenta[j]=help;
      }
      else {
	break;
      }
    }
  }
}

void Jet_Y::Evaluate(const ATOOLS::Particle_List &particles,double weight) 
{
  std::vector<ATOOLS::Vec4D> jetmomenta;
  if (m_mode>0) {
    p_jetfinder->ConstructJets(&particles,(*p_jets),(*p_ys));
    jetmomenta=p_jetfinder->JetMomenta();
  }
  else {
    jetmomenta.resize(particles.size());
    for (unsigned int i=0;i<particles.size();++i) jetmomenta[i]=particles[i]->Momentum();
  }
  SortJetPT(jetmomenta);
  for (unsigned int i=0;i<ATOOLS::Min(p_histogram->size(),jetmomenta.size());i+=2) {
    ATOOLS::Vec4D &p=jetmomenta[i];
    histo->Insert(log((p[1]+p[3])/(p[1]-p[3])),weight);
    (*p_histogram)[i/2]->Insert(log((p[1]+p[3])/(p[1]-p[3])),weight);
  }
}

void Jet_Y::Evaluate(const ATOOLS::Blob_List &blobs,double weight) 
{
  ATOOLS::Particle_List particles;
  for (ATOOLS::Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) { 
    for (int i=0;i<(*bit)->NOutP();++i) {
      ATOOLS::Particle *cur=(*bit)->OutParticle(i);
      if (cur->Status()==1) particles.push_back(cur);
    }
  }
  Evaluate(particles,weight);
}

void Jet_Y::EndEvaluation() 
{
  ATOOLS::msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  histo->Finalize();
  histo->Output();
  for (unsigned int i=0; i<p_histogram->size();++i) {
    (*p_histogram)[i]->Finalize();
    (*p_histogram)[i]->Output();
  }
  ATOOLS::msg.Debugging()<<std::endl;
}

void Jet_Y::Output(std::string pathname) 
{
  mkdir(pathname.c_str(),448); 
  std::string filename=pathname+std::string("/")+name;
  histo->Output((filename+std::string("_.ath")).c_str());
  for (unsigned int i=0; i<p_histogram->size();++i) {
    (*p_histogram)[i]->Output((filename+ATOOLS::ToString(i)+std::string(".ath")).c_str());
  }
}

