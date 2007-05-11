#include "Object_Definition_Base.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"
#include "MyStrStream.H"
#include "Exception.H"


using namespace ANALYSIS;
using namespace ATOOLS;

Object_Definition_Base::
Object_Definition_Base(Primitive_Analysis * ana,const std::string name,const std::string mode) :
  m_elements(false), m_name(name), p_myparticles(NULL),
  p_tracker(NULL), p_ECal(NULL), p_HCal(NULL), p_chambers(NULL)
{
  p_ana = ana;
  Analysis_Object * det = p_ana->GetObject("Detector");
  if (det==NULL) {
    det = new Detector();
    p_ana->AddObject(det);
  }
  dynamic_cast<Detector *>(det)->AddObjectDefinition(this);

  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
}

Object_Definition_Base::~Object_Definition_Base() {
  if (p_order)       { delete p_order; p_order = NULL; }
  if (p_myparticles) { p_myparticles->Clear(); delete p_myparticles; p_myparticles = NULL; }
}

void Object_Definition_Base::SetOrdering(const std::string mode) {
  if (p_order) delete p_order;
  p_order = Order_Getter::GetObject(mode,"");
  if (p_order==NULL) 
    THROW(fatal_error,"Invalid ordering mode '"+mode+"'");
}

void Object_Definition_Base::GetElements() {
  if (m_elements) return;
  if (p_tracker==NULL || p_ECal==NULL || p_HCal==NULL) {
    Detector * det = (dynamic_cast<Detector *>(p_ana->GetObject("Detector")));
    p_tracker  = dynamic_cast<Tracker *>(det->GetElement("Tracker"));
    p_ECal     = dynamic_cast<ElMag_Calorimeter *>(det->GetElement("ECal"));
    p_HCal     = dynamic_cast<Had_Calorimeter *>(det->GetElement("HCal"));
    p_chambers = dynamic_cast<Muon_Chambers *>(det->GetElement("Muon_Chambers"));
  }
  if (p_tracker==NULL || p_ECal==NULL || p_HCal==NULL) {
    msg.Error()<<"ERROR in "<<METHOD<<":"<<std::endl
	       <<"   Not all detector components found.  Will abort."<<std::endl;
    abort();
  }
  m_elements = true;
}

void Object_Definition_Base::Reset() {
  //p_myparticles->Clear();
}

void Object_Definition_Base::DropUsedCells() {
  std::vector<Cell *> * cells;
  std::list<Cell *> * ecalcells(NULL), * hcalcells(NULL);
  std::list<Cell *>::iterator cit;
  if (p_ECal) ecalcells = p_ECal->GetHitCells();
  if (p_HCal) hcalcells = p_HCal->GetHitCells();

  Cell * cell;
  bool found;
  for (ObjectListIterator olit=m_objects.begin();olit!=m_objects.end();olit++) {
    cells = (*olit)->GetCells();
    if (!cells || cells->size()==0) continue;
    else {
      for (size_t c=0;c<cells->size();c++) {
	found = false;
	cell  = (*cells)[c];
	for (cit=ecalcells->begin();cit!=ecalcells->end();cit++) {
	  if (cell==(*cit)) { ecalcells->erase(cit); found = true; break; }
	}
	if (found) continue;
	for (cit=hcalcells->begin();cit!=hcalcells->end();cit++) {
	  if (cell==(*cit)) { hcalcells->erase(cit); found = true; break; }
	}
	if (found) continue;
      }
    }
  }
}
