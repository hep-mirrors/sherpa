#include "Primitive_Detector.H"

#include "Primitive_Calorimeter.H"
#include "Primitive_Analysis.H"


#include "MyStrStream.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;


Primitive_Detector::Primitive_Detector(const det_mode::code mode,const std::string &inlist): 
  m_mode(mode), 
  m_inlistname(inlist), m_internalparticles(new Particle_List)
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

  for (String_ObjectDefinition_Iter soditer=m_objectdefs.begin();
       soditer!=m_objectdefs.end();++soditer) {
    if (soditer->second!=NULL) {
      m_objectdefs.erase(soditer--);
    }
  }
  m_objectdefs.clear();

  for (std::deque<Particle*>::iterator pit=m_internalparticles->begin();
       pit!=m_internalparticles->end();pit++) {
    if (*pit!=NULL) { delete (*pit); (*pit)=NULL; }
  }
  m_internalparticles->clear();
  delete m_internalparticles;
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

void Primitive_Detector::SetName(std::string name) { m_name = name; }

std::string Primitive_Detector::Name() const { return m_name; }

void Primitive_Detector::AddDetectorElement(Primitive_Detector_Element * pde)
{
  m_elements[pde->Name()] = pde;
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

void Primitive_Detector::AddObjectDefinition(Object_Definition_Base * obj)
{
  m_objectdefs[obj->Name()] = obj;
}

Object_Definition_Base * Primitive_Detector::GetObjectDefinition(std::string name)
{
  String_ObjectDefinition_Iter soditer=m_objectdefs.find(name);
  if (soditer==m_objectdefs.end()) {
    msg.Error()<<"Potential Error in Primitive_Detector::GetObjectDefinition("<<name<<") :"<<std::endl
	       <<"   Element not found, return NULL and hope for the best."<<std::endl;
    return NULL;
  }
  return soditer->second;
}

void Primitive_Detector::Reset() {
  for (String_DetectorElement_Iter sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();sdeiter++) {
    sdeiter->second->Reset();
  }
  for (std::deque<Particle*>::iterator pit=m_internalparticles->begin();
       pit!=m_internalparticles->end();pit++) {
    if (*pit!=NULL) { delete (*pit); (*pit)=NULL; }
  }
  m_internalparticles->clear();
}

Primitive_Observable_Base* Primitive_Detector::Copy() const
{
  std::cout<<"WARNING: Potential error in Primitive_Detector: "
	   <<"No appropriate Copy() method.\n"
	   <<"   Continue and hope for the best"<<std::endl;
  
  Primitive_Detector *detector(new Primitive_Detector(m_mode,m_inlistname));

  for (String_DetectorElement_Map::const_iterator sdeiter=m_elements.begin();
       sdeiter!=m_elements.end();++sdeiter) {
    if (sdeiter->second!=NULL) {
      detector->AddDetectorElement(dynamic_cast<Primitive_Detector_Element *>(sdeiter->second->GetCopy()));
    }
  }
  for (String_ObjectDefinition_Map::const_iterator soditer=m_objectdefs.begin();
       soditer!=m_objectdefs.end();++soditer) {
    if (soditer->second!=NULL) {
      detector->AddObjectDefinition(dynamic_cast<Object_Definition_Base *>(soditer->second->GetCopy()));
    }
  }

  return detector;
}

void Primitive_Detector::Evaluate(const ATOOLS::Blob_List &bloblist,
				  double value,int ncount)
{
  Particle_List * anaparticles=p_ana->GetParticleList(m_inlistname);
  if (anaparticles==NULL) {
    ATOOLS::msg.Error()<<"Primitive_Detector::Evaluate(..): "
		       <<"Particle list '"<<m_inlistname<<"' not found."<<std::endl;
    return;
  }
  for (std::deque<Particle*>::iterator pit=anaparticles->begin();
       pit!=anaparticles->end();pit++) {
    m_internalparticles->push_back(new Particle(**pit));
  }
  switch (m_mode) {
  case det_mode::object_reco:
    FillDetector();
    ProduceReconstructedObjectLists();
    break;
  case det_mode::simple_det:
    FillDetector();
    ProduceSimpleDetectorObjectLists();
    break;
  case det_mode::MC_truth:
  default:
    ProduceMCTruthObjectLists();
    break;
  }
}

void Primitive_Detector::FillDetector() {
  if (GetElement("Tracker"))      GetElement("Tracker")->Fill(m_internalparticles);
  if (GetElement("ECal"))         GetElement("ECal")->Fill(m_internalparticles);
  if (GetElement("HCal"))         GetElement("HCal")->Fill(m_internalparticles);
  if (GetElement("MuonChambers")) GetElement("MuonChambers")->Fill(m_internalparticles);
}


void Primitive_Detector::ProduceMCTruthObjectLists() {
  for (String_ObjectDefinition_Iter soditer=m_objectdefs.begin();
       soditer!=m_objectdefs.end();soditer++) {
    soditer->second->FillMCTruthList(m_internalparticles,this);
  }
}

void Primitive_Detector::ProduceSimpleDetectorObjectLists() {
  for (String_ObjectDefinition_Iter soditer=m_objectdefs.begin();
       soditer!=m_objectdefs.end();soditer++) {
    soditer->second->FillSimpleDetectorList(m_internalparticles,this);
  }
}

void Primitive_Detector::ProduceReconstructedObjectLists() {
}


