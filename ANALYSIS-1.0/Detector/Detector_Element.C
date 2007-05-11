#include "Detector_Element.H"
#include "Primitive_Analysis.H"
#include "Detector.H"
#include "Message.H"

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
  msg.Out()<<"Hits in "<<m_name<<" : "<<std::endl;
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

std::list<Track *> * Detector_Element::GetTracks(const double eta,const double phi,
						 const double R2,kf::code kfcode) {
  std::list<Track *> * tracks(new std::list<Track *>);
  //std::cout<<"   ... "<<METHOD<<" : from "<<this<<", tracks =  "
  //<<tracks<<", local = "<<(&m_tracks)<<std::endl;
  if (eta==0&&phi==0&&R2==0) (*tracks) = m_tracks;
  else {
    //std::cout<<METHOD<<": "<<m_tracks.size()<<std::endl;
    
    for (std::list<Track *>::iterator trit=m_tracks.begin();trit!=m_tracks.end();trit++) {
      //if (kfcode!=kf::none && (*trit)->flav.Kfcode()==kfcode) {
      // std::cout<<"  "<<(*trit)->flav<<"  "
      //	       <<(*trit)->eta<<" ("<<eta<<"), "<<(*trit)->phi<<" ("<<phi<<")"
      //	       <<"  --> "<<(sqr((*trit)->eta-eta)+sqr((*trit)->phi-phi))<<" vs. "<<R2<<std::endl;
      //}
      if (sqr((*trit)->eta-eta)+sqr((*trit)->phi-phi)>R2) continue;
      tracks->push_back(*trit);
    }
  }
  return tracks;
}
