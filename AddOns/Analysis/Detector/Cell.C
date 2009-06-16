#include "AddOns/Analysis/Detector/Cell.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Cell::Cell(const double* dims,Etastrip * etastrip) :
  m_etamin(dims[0]),m_etamax(dims[1]),m_phimin(dims[2]),m_phimax(dims[3]),
  m_eta((m_etamin+m_etamax)/2.),m_phi((m_phimin+ m_phimax)/2.),
  m_costheta((exp(m_eta)-exp(-m_eta))/(exp(m_eta)+exp(-m_eta))),
  m_sintheta(sqrt(1.-m_costheta*m_costheta)),m_summedE(0.),
  m_used(false),
  p_etastrip(etastrip),p_up(NULL),p_down(NULL),
  m_direction(Vec4D(1.,cos(m_phi)*m_sintheta,sin(m_phi)*m_sintheta,m_costheta))
{
}

Cell::~Cell() {
  //exh->GenerateStackTrace(std::cout);
  //msg_Info();
  //PRINT_INFO(this);
}

void Cell::Reset() {
  if (!m_energydeposits.empty()) {
    for (std::map<Particle *,double>::iterator part=m_energydeposits.begin();
	 part!=m_energydeposits.end();part++) {
      if (part->first) { delete part->first; }
    }
    m_energydeposits.clear();
  }
  m_summedE = 0.;
  m_used    = false;
}

void Cell::Dimensions(double * dims) const {
  dims[0] = m_etamin; dims[1] = m_etamax;
  dims[2] = m_phimin; dims[3] = m_phimax;
}

void Cell::Centroid(double & eta,double & phi) const {
  eta = m_eta;
  phi = m_phi;
}

void Cell::AddParticle(Particle * part,double deposit) {
  if (deposit<0.) deposit = part->Momentum()[0];
  Particle * prt = new Particle(*part);
  prt->SetOriginalPart(part);
  m_summedE += m_energydeposits[prt] = deposit;
}

Cell * Cell::AddDeposit(const double phi,const double dep) {
  Cell * cell(LocateCell(phi));
  if (cell) cell->AddDeposit(dep);
  return cell;
}

Cell * Cell::AddParticle(const double phi,Particle * part,double dep) {
  Cell * cell(LocateCell(phi));
  if (cell) cell->AddParticle(part,dep);
  return cell;
}

Cell * Cell::LocateCell(double phi) {
  while (phi>M_PI)  phi-=2.*M_PI;
  while (phi<-M_PI) phi+=2.*M_PI;
  if (phi>=m_phimin) {
    if (phi<=m_phimax) return this;
    if (p_up) return p_up->LocateCell(phi);
         else return NULL;
  }
  if (p_down) return p_down->LocateCell(phi);
  return NULL;
}

double Cell::R2(const double eta,const double phi) {
  return sqr(eta-m_eta)+sqr(phi-m_phi);
}

Vec4D Cell::TrueMom() const {
  Vec4D truemom(Vec4D(0.,0.,0.,0.));
  for (std::map<Particle *,double>::const_iterator part=m_energydeposits.begin();
       part!=m_energydeposits.end();part++) {
    truemom += part->first->Momentum();
  }
  return truemom;
}

void Cell::Print() const {
  if (m_energydeposits.size()>0) {
    msg_Out()<<"Cell ("<<this<<" : "<<m_eta<<","<<m_phi<<"), E = "<<m_summedE
	     << " from ("<<m_energydeposits.size()<<"):"<<std::endl;
    for (std::map<Particle *,double>::const_iterator part=m_energydeposits.begin();
	 part!=m_energydeposits.end();part++) {
      msg_Out()<<"   "<<part->first<<std::endl;
      msg_Out()<<"   "<<part->first->Flav()<<" "<<part->first->Momentum()
	       <<" ["<<part->first->Momentum().Eta()<<","<<part->first->Momentum().Phi()<<"]"
	       <<" -> "<<part->second<<" = "<<part->first->Momentum()[0]/part->second<<"."<<std::endl;
    }
  }
  else
    msg_Out()<<"Cell ("<<this<<" : "<<m_eta<<","<<m_phi<<"), is emtpy. "
	     <<"E = "<<m_summedE<<"."<<std::endl;
}


Etastrip::Etastrip(const double etamin,const double etamax,const long int nphi,
		   Detector_Segment * segment) :
  m_etamin(etamin),m_etamax(etamax),m_nphi(nphi),
  p_segment(segment),p_plus(NULL),p_minus(NULL)
{
  double deltaphi(2.*M_PI/double(m_nphi));
  double * dims = new double[4];
  dims[0] = m_etamin;
  dims[1] = m_etamax;
  dims[2] = -M_PI;
  dims[3] = dims[2]+deltaphi;
  p_zero  = new Cell(dims,this);

  Cell * cell,* pref(p_zero);
  for (long int i=1;i<m_nphi;i++) {
    dims[2] = -M_PI+i*deltaphi;
    dims[3] = dims[2]+deltaphi;
    cell    = new Cell(dims,this);
    cell->SetDown(pref);
    pref->SetUp(cell);
    pref = cell;
  }

  pref->SetUp(p_zero);
  p_zero->SetDown(pref);

  delete [] dims;
}

Etastrip::~Etastrip() {
  Cell * cell(p_zero), * stop(p_zero->GetDown());
  do {
    cell = cell->GetUp();
    delete cell->GetDown();
  } while (cell!=stop);
  delete stop;
}

void Etastrip::Reset() {
  Cell * cell(p_zero);
  do {
    cell->Reset();
    cell = cell->GetUp();
  } while (cell!=p_zero);
}

void Etastrip::Dimensions(double & etamin,double & etamax) const {
  etamin = m_etamin;
  etamax = m_etamax;
}

Cell * Etastrip::AddDeposit(const double eta,const double phi,const double dep) {
  Cell * cell(LocateCell(eta,phi));
  if (cell) cell->AddDeposit(dep);
  return cell;
}

Cell * Etastrip::AddParticle(const double eta,const double phi,Particle * part,double dep) {
  Cell * cell(LocateCell(eta,phi));
  if (cell) cell->AddParticle(part,dep);
  return cell;
}

Cell * Etastrip::LocateCell(const double eta,const double phi) {
  Etastrip * etastrip(LocateEtastrip(eta));
  if (etastrip) return etastrip->GetZero()->LocateCell(phi);
  return NULL;
}

Etastrip * Etastrip::LocateEtastrip(const double eta) {
  if (eta>=m_etamin) {
    if (eta<=m_etamax) return this;
    if (p_plus) return p_plus->LocateEtastrip(eta);
         else return NULL;
  }
  if (p_minus) return p_minus->LocateEtastrip(eta);
  return NULL;
}


