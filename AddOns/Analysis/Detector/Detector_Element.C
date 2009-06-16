#include "AddOns/Analysis/Detector/Detector_Element.H"
#include "AddOns/Analysis/Main/Primitive_Analysis.H"
#include "AddOns/Analysis/Detector/Detector.H"
#include "ATOOLS/Org/Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Detector_Element::Detector_Element(Primitive_Analysis * ana,const std::string name) : 
  m_name(name),m_istested(false) 
{
  p_ana = ana;
  Analysis_Object * det = p_ana->GetObject("Detector");
  if (det==NULL) {
    det = new Detector(p_ana);
    p_ana->AddObject(det);
  }
  dynamic_cast<Detector *>(det)->AddDetectorElement(this);
}

Detector_Element::Detector_Element(Detector_Segment * seg) 
{
  m_segments.insert(seg);
}

Detector_Element::~Detector_Element() {
  while (m_segments.size()>0) {
    delete (*m_segments.begin());
    m_segments.erase(m_segments.begin());
  }
  m_segments.clear();
}

void Detector_Element::AddDetectorSegment(Detector_Segment * seg) { 
  m_segments.insert(seg);
  std::set<Detector_Segment *,DS_Order>::iterator ds,ds1;
  for (ds=m_segments.begin();ds!=m_segments.end();ds++) {
    if ((*ds)==seg) break;
  }
  if (ds==m_segments.end() || m_segments.size()<2) return;
  //std::cout<<METHOD<<" : "<<seg<<" vs. "<<(*ds)<<" in "<<m_segments.size()<<std::endl;
  ds1=ds; 
  if (ds1++!=m_segments.end()) {
    (*ds)->GetLast()->SetPlus((*ds1)->GetFirst());
    (*ds1)->GetFirst()->SetMinus((*ds)->GetLast());
  }
  if (ds!=m_segments.begin()) {
    ds1=ds;
    ds1--;
    (*ds1)->GetLast()->SetPlus((*ds)->GetFirst());
    (*ds)->GetFirst()->SetMinus((*ds1)->GetLast());
  }
}

void Detector_Element::Evaluate(const ATOOLS::Blob_List &,double,int) {}

std::string Detector_Element::Name() { return m_name; }

void Detector_Element::PrintHits() const {
  msg_Out()<<"Hits in "<<m_name<<" : "<<std::endl;
  for (std::list<Cell *>::const_iterator cit=m_hitcells.begin();
       cit!=m_hitcells.end();cit++) (*cit)->Print();
}

Cell * Detector_Element::GetCell(const double eta,const double phi) {
  double etamin,etamax;
  for (std::set<Detector_Segment *,DS_Order>::iterator ds=m_segments.begin();
       ds!=m_segments.end();ds++) {
    (*ds)->Dimensions(etamin,etamax);
    if ((*ds)->InEtaRange(eta)) return (*ds)->LocateCell(eta,phi);
  }
  return NULL;
}

void Detector_Element::GetTracks(std::list<Track *> & tracks,const double eta,const double phi,
				 const double R2,kf_code kfcode) {
  tracks.clear();
  for (std::list<Track *>::iterator trit=m_tracks.begin();trit!=m_tracks.end();trit++) {
    if ( (eta==0&&phi==0&&R2==0) ||
	 (sqr((*trit)->eta-eta)+sqr((*trit)->phi-phi)<=R2) ) tracks.push_back(*trit);
  }
}
