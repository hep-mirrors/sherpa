#include "AddOns/Analysis/Detector/Detector_Segment.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Detector_Segment::Detector_Segment(const double etamin,const double etamax,
				   const long int neta,const long int nphi) :
  m_etamin(etamin),m_etamax(etamax),m_neta(neta),m_nphi(nphi),
  p_first(NULL),p_last(NULL)
{
  double deltaeta((m_etamax-m_etamin)/double(m_neta));
  p_first = new Etastrip(m_etamin,m_etamin+deltaeta,nphi,this);
  Etastrip * etastrip(NULL),* pref(p_first);

  double eta;
  for (long int i=1;i<m_neta;i++) {
    eta      = etamin+i*deltaeta;
    etastrip = new Etastrip(eta,eta+deltaeta,nphi,this);
    etastrip->SetMinus(pref);
    pref->SetPlus(etastrip);
    pref = etastrip;
  }
  p_last = etastrip;
}

Detector_Segment::~Detector_Segment() {
  if (p_first==p_last) { delete p_first; return; }
  Etastrip * etastrip(p_first);
  do {
    etastrip = etastrip->GetPlus();
    delete etastrip->GetMinus();
  } while (etastrip!=p_last);
  delete p_last;
}

Analysis_Object * Detector_Segment::GetCopy() const {
  return new Detector_Segment(m_etamin,m_etamax,m_neta,m_nphi);
}

void Detector_Segment::Dimensions(double & etamin,double & etamax) const {
  etamin = m_etamin;
  etamax = m_etamax;
}

void Detector_Segment::Reset() {
  Etastrip * etastrip(p_first);
  do {
    etastrip->Reset();
    if (etastrip==p_last) break;
    etastrip = etastrip->GetPlus();
  } while (true);
}

Cell * Detector_Segment::LocateCell(const double eta,const double phi) { 
  return p_first->LocateCell(eta,phi); 
}

Cell * Detector_Segment::AddDeposit(const double eta,const double phi,const double dep) {
  return p_first->AddDeposit(eta,phi,dep);
}

Cell * Detector_Segment::AddParticle(const double eta,const double phi,
				     Particle * part,double dep) {
  return p_first->AddParticle(eta,phi,part,dep);
}
