#include "Jet_Multiplicity.H"

using namespace AMISIC;

Jet_Multiplicity::Jet_Multiplicity(int _m_type,double _m_min,double _m_max,
				   int _m_nbins,double _m_ycut,double _m_ptcut):
  Primitive_Observable_Base(_m_type,_m_min,_m_max,_m_nbins,NULL),
  p_jets(new std::vector<int>((unsigned int)_m_max)),
  p_ys(new std::vector<double>()),
  m_njets((unsigned int)_m_max),
  m_mode(0),
  m_ycut(_m_ycut),
  m_ptcut(_m_ptcut)
{
  histo = new ATOOLS::Histogram(type,xmin,xmax,nbins);
  for (unsigned int i=0;i<m_njets;++i) (*p_jets)[i]=i;
  name=std::string("Jet_Multiplicity");
}

Jet_Multiplicity::~Jet_Multiplicity()
{
  delete p_jets;
  delete p_ys;
}

void Jet_Multiplicity::Evaluate(const ATOOLS::Particle_List &particles,double weight) 
{
  if (m_mode>0) {
    p_jetfinder->ConstructJets(&particles,(*p_jets),(*p_ys));
    for (unsigned int i=0;i<p_ys->size();++i) {
      if ((*p_ys)[i]>=m_ycut) {
	histo->Insert((*p_jets)[i],weight);
	return;
      }
    }
  }
  else {
    unsigned int jets=0;
    for (unsigned int i=0;i<particles.size();++i) {
      ATOOLS::Vec4D p=particles[i]->Momentum();
      if (sqrt(p[1]*p[1]+p[2]*p[2])>m_ptcut) ++jets;
    }
//     std::cout<<" found "<<jets<<std::endl;
    if (jets==5) {
      for (unsigned int i=0;i<particles.size();std::cout<<particles[i++]<<std::endl); 
    }
    histo->Insert(jets,weight);
  }
}

void Jet_Multiplicity::Evaluate(const ATOOLS::Blob_List &blobs,double weight) 
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

void Jet_Multiplicity::EndEvaluation() 
{
  ATOOLS::msg.Debugging()<<"  "<<name<<" : "<<std::endl;
  histo->Finalize();
  histo->Output();
  ATOOLS::msg.Debugging()<<std::endl;
}

void Jet_Multiplicity::Output(std::string pathname) 
{
  mkdir(pathname.c_str(),448); 
  std::string filename=pathname+std::string("/")+name+std::string(".ath");
  histo->Output(filename.c_str());
}
