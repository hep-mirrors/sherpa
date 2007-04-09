#include "Primitive_Detector_Element.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;

Primitive_Detector_Element::
Primitive_Detector_Element(const int neta,const int nphi,
			   const double etamin,const double etamax,
			   const double phimin,const double phimax) :
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

void Primitive_Detector_Element::GetDimensions(int & neta,int & nphi,
					       double & mineta, double & maxeta,
					       double & minphi, double & maxphi) { 
  neta     = m_neta;     nphi     = m_nphi;
  mineta   = m_etamin;   maxeta   = m_etamax;
  mineta   = m_phimin;   maxeta   = m_phimax;
}

double Primitive_Detector_Element::Cell(const int i,const int j) const {
  if (i>-1&&j>-1&&i<m_neta&&j<m_nphi) return p_cells[i][j];
  else msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
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

void Primitive_Detector_Element::PseudoRapidityNAzimuthalAngle(const Vec4D & p,
							       double & eta,double & phi) {
  eta = p.Eta();
  phi = p.Phi();
}

void Primitive_Detector_Element::PseudoRapidityNAzimuthalAngle(const int i,const int j,
							       double & eta, double & phi) {
  if (i<0||i>m_neta||j<0||j>m_nphi) {
    ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
		       <<"   PseudoRapidityNAzimuthalAngle("<<i<<","<<j<<") "
		       <<"out of bounds, continue and leave eta/phi unchanged."<<std::endl;
    return;
  }
  eta = m_etamin+m_delta_eta/2.+i*m_delta_eta;
  phi = m_phimin+m_delta_phi/2.+j*m_delta_phi;    
}

ATOOLS::Vec4D Primitive_Detector_Element::ReconstructMasslessFourMom(const int i,const int j) {
  if (i<0||i>m_neta||j<0||j>m_nphi) {
    ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
		       <<"   PseudoRapidityNAzimuthalAngle("<<i<<","<<j<<") "
		       <<"out of bounds, continue and leave eta/phi unchanged."<<std::endl;
    return Vec4D(0.,0.,0.,0.);
  }
  double E(Cell(i,j)),eta,phi; 
  PseudoRapidityNAzimuthalAngle(i,j,eta,phi);
  double plong(E*(exp(eta)-exp(-eta))/(exp(eta)+exp(-eta))),pperp(sqrt(E*E-plong*plong));
  return Vec4D(E,pperp*cos(phi),pperp*sin(phi),plong);
}


void Primitive_Detector_Element::SetName(std::string name) { m_name = name; }

std::string Primitive_Detector_Element::Name() const { return m_name; }
