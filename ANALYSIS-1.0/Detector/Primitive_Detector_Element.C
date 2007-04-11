#include "Primitive_Detector_Element.H"
#include "Message.H"

using namespace ANALYSIS;
using namespace ATOOLS;



Primitive_Detector_Element::
Primitive_Detector_Element(const pde::code mode,const long int neta,const long int nphi,
			   const double etamin,const double etamax,
			   const double phimin,const double phimax) :
  m_mode(mode),m_neta(neta), m_nphi(nphi), 
  m_etamin(etamin), m_etamax(etamax),   m_phimin(phimin), m_phimax(phimax),
  m_delta_eta((m_etamax-m_etamin)/double(m_neta)), 
  m_delta_phi((m_phimax-m_phimin)/double(m_nphi)),
  p_qualifier(NULL)
{
  switch (m_mode) {
    case pde::tracks:
      break;
    case pde::cells:
    default:
      m_cells.resize(m_neta);
      for (m_etastrip=m_cells.begin();m_etastrip!=m_cells.end();m_etastrip++)
	(*m_etastrip).resize(m_nphi);
      break;
  }
}

Primitive_Detector_Element::~Primitive_Detector_Element()
{
  if (p_qualifier) { delete p_qualifier; p_qualifier=NULL; }
}

void Primitive_Detector_Element::Reset()
{
  switch (m_mode) {
    case pde::tracks:
      m_tracks.clear();
      break;
    case pde::cells:
    default:
      for (m_etastrip=m_cells.begin();m_etastrip!=m_cells.end();m_etastrip++) {
	for (std::vector<double>::iterator etaphi=(*m_etastrip).begin();
	     etaphi!=(*m_etastrip).end();etaphi++) 
	  (*etaphi) = 0.;
      }
      break;
  }
  AddNoise();
}

void Primitive_Detector_Element::GetNumbersOfCells(long int & neta,long int & nphi) {
  neta     = m_neta;     nphi     = m_nphi;
}

void Primitive_Detector_Element::GetDimensions(long int & neta,long int & nphi,
					       double & mineta, double & maxeta,
					       double & minphi, double & maxphi) { 
  neta     = m_neta;     nphi     = m_nphi;
  mineta   = m_etamin;   maxeta   = m_etamax;
  mineta   = m_phimin;   maxeta   = m_phimax;
}

double Primitive_Detector_Element::Cell(const long int i,const long int j) const {
  if (i>-1&&j>-1&&i<m_neta&&j<m_nphi) return m_cells[i][j];
  else msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
		  <<"   GetCell("<<i<<","<<j<<") out of bounds, return 0."<<std::endl;
  return 0.;
}

void Primitive_Detector_Element::MatchCell(const double eta,const double phi,
					   long int & i, long int & j) const {
  i = (long int)(m_neta*(eta-m_etamin)/(m_etamax-m_etamin));
  j = (long int)(m_nphi*(phi-m_phimin)/(m_phimax-m_phimin));
  if (i<0 || i>m_neta) i=-1;
  if (j<0 || j>m_neta) j=-1;
}

void Primitive_Detector_Element::MatchCell(const ATOOLS::Vec4D p,
					   long int & i,long int & j) const {
  MatchCell(p.Eta(),p.Phi(),i,j);
}

void Primitive_Detector_Element::PseudoRapidityNAzimuthalAngle(const Vec4D & p,
							       double & eta,double & phi) {
  eta = p.Eta();
  phi = p.Phi();
}

void Primitive_Detector_Element::PseudoRapidityNAzimuthalAngle(const long int i,const long int j,
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

ATOOLS::Vec4D Primitive_Detector_Element::
ReconstructMasslessFourMom(const long int i,const long int j) {
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
