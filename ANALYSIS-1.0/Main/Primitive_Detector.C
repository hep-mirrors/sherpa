#include "Primitive_Detector.H"

#include "Primitive_Calorimeter.H"
#include "Primitive_Analysis.H"


#include "MyStrStream.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

#include <iomanip>

DECLARE_GETTER(Primitive_Detector_Getter,"Detector",
	       Primitive_Observable_Base,String_Matrix);

void Primitive_Detector_Getter::PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"{\n"
     <<std::setw(width+7)<<" "<<"InList  list\n"
     <<std::setw(width+7)<<" "<<"OutList list\n"
     <<std::setw(width+7)<<" "<<"HadCal  etamin etamax etacells phicells\n"
     <<std::setw(width+7)<<" "<<"CalCone etmin etamin etamax deltar [bjets]\n"
     <<std::setw(width+4)<<" "<<"}";
}

Primitive_Observable_Base *
Primitive_Detector_Getter::operator()(const String_Matrix &parameters) const
{
  std::string inlist="FinalState", outlist="Detected";
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="InList" && cur.size()>1) inlist=cur[1];
    else if (cur[0]=="OutList" && cur.size()>1) outlist=cur[1];
  }
  Primitive_Detector *detector = new Primitive_Detector(inlist,outlist);
  for (size_t i=0;i<parameters.size();++i) {
    const std::vector<std::string> &cur=parameters[i];
    if (cur[0]=="HadCal" && cur.size()>4) {
      detector->Add(new Primitive_Calorimeter(ATOOLS::ToType<double>(cur[1]),
					      ATOOLS::ToType<double>(cur[2]),
					      ATOOLS::ToType<int>(cur[3]),
					      ATOOLS::ToType<int>(cur[4])));
    }
    else if (cur[0]=="CalCone" && cur.size()>4) {
      msg.Out()<<"WARNING CalCone   no longer supported by Primitive Detector ! "<<std::endl;
    }
  }
  return detector;
}

#include "Final_Selector.H"

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

Primitive_Detector::Primitive_Detector(const std::string &inlist,
				       const std::string &outlist): 
  m_inlistname(inlist), m_outlistname(outlist)
{ 
  m_name = std::string("Full Detector"); 
  m_splitt_flag = false;
}

Primitive_Detector::~Primitive_Detector()
{
  for (String_DetectorElement_Iter sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();++sdeiter) {
    if (sdeiter->second!=NULL) {
      m_elements.erase(sdeiter--);
    }
  }
  m_elements.clear();
}

Primitive_Observable_Base* Primitive_Detector::Copy() const
{
  std::cout<<"WARNING: Potential error in Primitive_Detector: "
	   <<"No appropriate Copy() method.\n"
	   <<"   Continue and hope for the best"<<std::endl;
  
  Primitive_Detector *detector = 
    new Primitive_Detector(m_inlistname,m_outlistname);

  for (String_DetectorElement_Map::const_iterator sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();++sdeiter) {
    if (sdeiter->second!=NULL) {
      detector->Add(sdeiter->second->Copy());
    }
  }

  return detector;
}

void Primitive_Detector::Evaluate(const ATOOLS::Blob_List &bloblist,
				  double value,int ncount)
{
  Particle_List *inparticles=p_ana->GetParticleList(m_inlistname);
  if (inparticles==NULL) {
    ATOOLS::msg.Error()<<"Primitive_Detector::Evaluate(..): "
		       <<"Particle list '"<<m_inlistname<<"' not found."<<std::endl;
    return;
  }
  Fill(inparticles);
  Particle_List *outparticles = new Particle_List();
  p_ana->AddParticleList(m_outlistname,outparticles);
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
  for (Blob_List::const_iterator blit=bl->begin();blit!=bl->end();++blit) {
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
