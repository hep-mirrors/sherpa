#include "Primordial_KT.H"

#include "Scaling.H"

using namespace AMISIC;

Primordial_KT::Primordial_KT(int _m_type,double _m_min,double _m_max,int _m_nbins):
  Primitive_Observable_Base(_m_type,_m_min,_m_max,_m_nbins,NULL)
{
  histo = new ATOOLS::Histogram(type,xmin,xmax,nbins);
  name=std::string("Primordial_KT");
}

Primordial_KT::~Primordial_KT() {}

void Primordial_KT::Evaluate(const ATOOLS::Particle_List &particles,double weight) 
{
  for (unsigned int i=0;i<particles.size();++i) {
    if (particles[i]->ProductionBlob()!=NULL) {
      if (particles[i]->ProductionBlob()->Type().find("Beam Remnant")!=std::string::npos) {
	ATOOLS::Vec4D kt=particles[i]->Momentum();
	histo->Insert(sqrt(kt[1]*kt[1]+kt[2]*kt[2]),weight);
      }
    }
  }
}

void Primordial_KT::Evaluate(const ATOOLS::Blob_List &blobs,double weight) 
{
  ATOOLS::Particle_List particles;
  for (ATOOLS::Blob_Const_Iterator bit=blobs.begin();bit!=blobs.end();++bit) { 
    if ((*bit)->Type().find("Beam Remnant")!=std::string::npos) {
      for (int i=0;i<(*bit)->NOutP();++i) particles.push_back((*bit)->OutParticle(i));
    }
  }
  Evaluate(particles,weight);
}

void Primordial_KT::EndEvaluation() 
{
  histo->Finalize();
  histo->Output();
}

void Primordial_KT::Output(std::string pathname) 
{
  mkdir(pathname.c_str(),448); 
  std::string filename=pathname+std::string("/")+name;
  histo->Output((filename+std::string(".ath")).c_str());
}

