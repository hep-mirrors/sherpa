#include "Primitive_Detector.H"
#include "Particle_List.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Detector_Element::Primitive_Detector_Element() : 
  m_nx(-1), m_ny(-1), m_name(std::string("Unspecified")), p_cells(NULL) { }

Primitive_Detector_Element::Primitive_Detector_Element(std::string name) : 
  m_nx(-1), m_ny(-1), m_name(name), p_cells(NULL) { }

Primitive_Detector_Element::Primitive_Detector_Element(int nx,int ny,std::string name) :
  m_nx(nx), m_ny(ny), m_name(name)
{
  p_cells = new double*[m_nx];
  for (int i=0; i<m_nx;++i) p_cells[i] = new double[m_ny];
}

Primitive_Detector_Element::~Primitive_Detector_Element()
{
  if (p_cells) {
    for (int i=0;i<m_nx;++i) delete [] p_cells[i];
    p_cells=NULL;
  }
}

Primitive_Detector::Primitive_Detector() { m_name = std::string("Full Detector"); }

Primitive_Detector::~Primitive_Detector()
{
  for (String_DetectorElement_Iter sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();++sdeiter) {
    if (sdeiter->second!=NULL) delete [] sdeiter->second;
  }
  m_elements.clear();
}

void Primitive_Detector::Add(Primitive_Detector_Element * pde)
{
  m_elements[pde->Name()] = pde;
}

void Primitive_Detector::Print() {
  if (!ATOOLS::msg.LevelIsInfo()) return;
  ATOOLS::msg.Out()<<"==================================================="<<std::endl
		   <<m_name<<" with "<<m_elements.size()<<" components : "<<std::endl;
  int i=1;
  std::string name;
  for (String_DetectorElement_Iter sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();sdeiter++) {
    ATOOLS::msg.Out()<<"Element "<<i<<": "<<sdeiter->second->Name()<<std::endl;
  }
  ATOOLS::msg.Out()<<"==================================================="<<std::endl;
}

void Primitive_Detector::Fill(const ATOOLS::Blob_List * bl)
{
  Particle_List * pl = new Particle_List; 
  for (Blob_Const_Iterator blit=bl->begin();blit!=bl->end();++blit) {
    for (int i=0;i<(*blit)->NOutP();++i) {
      Particle * p = (*blit)->OutParticle(i);
      if (p->DecayBlob()==NULL) pl->push_back(new Particle(*p));
    }
  }
  Fill(pl);
  for (Particle_List::iterator pit=pl->begin();pit!=pl->end();++pit) delete (*pit);
  delete pl;
} 

void Primitive_Detector::Fill(const ATOOLS::Particle_List * pl) 
{
  for (String_DetectorElement_Iter sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();++sdeiter) {
    if (sdeiter->second!=NULL) sdeiter->second->Fill(pl);
  }
}

Primitive_Detector_Element * Primitive_Detector::GetElement(std::string name)
{
  String_DetectorElement_Iter sdeiter=m_elements.find(name);
  if (sdeiter==m_elements.end()) {
    msg.Error()<<"Potential Error in Primitive_Detector::GetElement("<<name<<") :"<<std::endl
	       <<"   Element not found, return NULL and hope for the best."<<std::endl;
    return NULL;
  }
  return sdeiter->second;
}
