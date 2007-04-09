#include "Primitive_Detector_Element.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Detector_Element(const int neta,const int nphi,
			   const double etamin=-5.,const double etamax=5.,
			   const double phimin=0.,const double phimax=2.*M_PI) :
  m_neta(neta), m_nphi(nphi), 
  m_etamin(etamin), m_etamax(etamax),   m_phimin(phimin), m_phimax(phimax),
  m_delta_eta((m_etamax-m_etamin)/double(m_neta)), 
  m_delta_phi((m_phimax-m_phimin)/double(m_nphi)),
  p_qualifier(NULL)
{
  p_cells = new double*[m_neta];
  for (int i=0; i<m_neta;++i) p_cells[i] = new double[m_nphi];
}

Primitive_Detector_Element::~Primitive_Detector_Element()
{
  if (p_cells) {
    for (int i=0;i<m_neta;++i) delete [] p_cells[i];
    p_cells=NULL;
  }
  if (p_qualifier) { delete p_qualifier; p_qualifier=NULL; }
}

void Primitive_Calorimeter::GetDimensions(int & neta,int & nphi,
						 double & mineta, double & maxeta,
						 double & minphi, double & maxphi) { 
  neta     = m_neta;     nphi     = m_nphi;
  mineta   = m_mineta;   maxeta   = m_maxeta;
  mineta   = m_minphi;   maxeta   = m_maxphi;
}

double Primitive_Detector_Element::Cell(const int i,const int j) const {
  if (i>-1&&j>-1&&i<m_neta&&j<m_nphi) return p_cells[i][j];
  else ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
			  <<"   GetCell("<<i<<","<<j<<") out of bounds, return 0."<<std::endl;
  return 0.;
}

void Primitive_Detector_Element::MatchCell(const double eta,const double phi,
						  int & i, int & j) const {
  i = int(m_neta*(eta-m_etamin)/(m_etamax-m_etamin));
  j = int(m_nphi*(phi-m_phimin)/(m_phimax-m_phimin));
  if (i<0 || i>m_neta) i=-1;
  if (j<0 || j>m_neta) j=-1;
}

void Primitive_Detector_Element::MatchCell(const ATOOLS::Vec4D p,
						  int & i,int & j) const {
  MatchCell(p.Eta(),p.Phi(),i,j);
}

void Detector_Element::PseudoRapidityNAzimuthalAngle(const Vec4D & p,
							    double & eta,double & phi) {
  eta = p.Eta();
  phi = p.Phi();
}

void Detector_Element::PseudoRapidityNAzimuthalAngle(const int i,const int j,
							    double & eta, double & phi) {
  if (i<0||i>m_neta||j<0||j>m_nphi) {
    ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
		       <<"   PseudoRapidityNAzimuthalAngle("<<i<<","<<j<<") "
		       <<"out of bounds, continue and leave eta/phi unchanged."<<std::endl;
    return;
  }
  eta = m_mineta+m_deltaeta/2.+i*m_deltaeta;
  phi = m_minphi+m_deltaphi/2.+j*m_deltaphi;    
}

ATOOLS::Vec4D Detector_Element::ReconstructMasslessFourMom(const int i,const int j) {
  if (i<0||i>m_neta||j<0||j>m_nphi) {
    ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
		       <<"   PseudoRapidityNAzimuthalAngle("<<i<<","<<j<<") "
		       <<"out of bounds, continue and leave eta/phi unchanged."<<std::endl;
    return;
  }
  double E(Cell(i,j)),eta,phi; 
  PseudoRapidityNAzimuthalAngle(i,j,eta,phi);
  double plong(E*(exp(eta)-exp(-eta))/(exp(eta)+exp(-eta))),pperp(sqrt(E*E-pz*pz));
  return Vec4D(E,pperp*cos(phi),pperp*sin(phi),plong);
}


void SetName(std::string name) { m_name = name; }

std::string Primitive_Detector_Element::Name() const { return m_name; }
