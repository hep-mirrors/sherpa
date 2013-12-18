#include "SHRiMPS/Beam_Remnants/Continued_PDF.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Continued_PDF::Continued_PDF(PDF::PDF_Base * pdf,
			     const Flavour & bunch) :
  p_pdf(pdf), m_bunch(bunch), 
  m_xmin(p_pdf->XMin()), m_xmax(p_pdf->XMax()), m_Q02(p_pdf->Q2Min()),
  m_geta(2.), m_glambda(-0.25)
{
  m_pdfpartons.push_back(Flavour(kf_u));
  m_pdfpartons.push_back(Flavour(kf_d));
  m_pdfpartons.push_back(Flavour(kf_s));
  m_pdfpartons.push_back(Flavour(kf_c));
  m_pdfpartons.push_back(Flavour(kf_b));
  m_pdfpartons.push_back(Flavour(kf_gluon));
  m_pdfpartons.push_back(Flavour(kf_u).Bar());
  m_pdfpartons.push_back(Flavour(kf_d).Bar());
  m_pdfpartons.push_back(Flavour(kf_s).Bar());
  m_pdfpartons.push_back(Flavour(kf_c).Bar());
  m_pdfpartons.push_back(Flavour(kf_b).Bar());
  for (std::list<Flavour>::iterator flit=m_pdfpartons.begin();
       flit!=m_pdfpartons.end();flit++) {
    m_xpdfmax[(*flit)] = m_xmaxpdf[(*flit)] = 0.;
  }
  CalculateNorms();
}

Continued_PDF::~Continued_PDF() {}

void Continued_PDF::CalculateNorms() 
{
  Sea_Kernel sea(p_pdf,m_bunch,&m_pdfpartons,m_Q02);
  Gauss_Integrator sintegrator(&sea);
  m_Snorm = sintegrator.Integrate(m_xmin,m_xmax,0.0001,1);
  Valence_Kernel val(p_pdf,m_bunch,&m_pdfpartons,m_Q02);
  Gauss_Integrator vintegrator(&val);
  m_Vnorm = vintegrator.Integrate(m_xmin,m_xmax,0.0001,1);
  m_gnorm = 
    exp(Gammln(m_geta+1.))*exp(Gammln(m_glambda+1.))/
    exp(Gammln(m_geta+m_glambda+2.));
  //std::cout<<"Gamma("<<0.5<<") = "<<exp(Gammln(0.5))<<", "
  //	   <<"Gamma("<<2<<") = "<<exp(Gammln(2))<<", "
  //	   <<"Gamma("<<4<<") = "<<exp(Gammln(4))<<".\n";
  //exit(1);
}

double Continued_PDF::XPDF(const Flavour & flav,const bool & defmax) {
  if (flav.IsDiQuark()) {
    if (m_bunch==Flavour(kf_p_plus) && !flav.IsAnti()) {
      return XPDF(Flavour(kf_u));
    }
    if (m_bunch==Flavour(kf_p_plus).Bar() && flav.IsAnti()) {
      return XPDF(Flavour(kf_u).Bar());
    }
    return 0.;
  }
  if (m_x<m_xmin) return 0.;
  if (m_Q2>m_Q02) return p_pdf->GetXPDF(flav);
  double seapart(0.), valpart(0.);
  if (m_bunch==Flavour(kf_p_plus)) {
    if (flav==Flavour(kf_u) || flav==Flavour(kf_d)) {
      seapart = p_pdf->GetXPDF(flav.Bar())*(m_Q2/m_Q02); 
      valpart = p_pdf->GetXPDF(flav)-p_pdf->GetXPDF(flav.Bar());
    }
    else if (flav==Flavour(kf_gluon)) {
      seapart = p_pdf->GetXPDF(flav)*(m_Q2/m_Q02); 
      valpart = 1./m_gnorm * pow(1.-m_x,m_geta) * pow(m_x,m_glambda) *
	m_Snorm * (1.-m_Q2/m_Q02);
      /*
	((p_pdf->GetXPDF(Flavour(kf_u))-p_pdf->GetXPDF(Flavour(kf_u).Bar()) +
	p_pdf->GetXPDF(Flavour(kf_d))-p_pdf->GetXPDF(Flavour(kf_d).Bar()))/  
	m_Vnorm) *
      */ 
    }
    else
      seapart = p_pdf->GetXPDF(flav)*(m_Q2/m_Q02);
  }
  else if (m_bunch==Flavour(kf_p_plus).Bar()) {
    if (flav==Flavour(kf_u).Bar() || flav==Flavour(kf_d).Bar()) {
      seapart = p_pdf->GetXPDF(flav.Bar())*(m_Q2/m_Q02); 
      valpart = p_pdf->GetXPDF(flav)-p_pdf->GetXPDF(flav.Bar());
    }
    else if (flav==Flavour(kf_gluon)) {
      seapart = p_pdf->GetXPDF(flav)*(m_Q2/m_Q02); 
      valpart = 
	(p_pdf->GetXPDF(Flavour(kf_u).Bar())-p_pdf->GetXPDF(Flavour(kf_u)) +
	 p_pdf->GetXPDF(Flavour(kf_d).Bar())-p_pdf->GetXPDF(Flavour(kf_d))) *  
	m_Snorm/m_Vnorm * (1.-m_Q2/m_Q02);
    }
    else
      seapart = p_pdf->GetXPDF(flav)*(m_Q2/m_Q02);
  }
  double total = seapart+valpart;
  if (defmax && total>m_xpdfmax[flav]) {
    m_xmaxpdf[flav] = m_x;
    m_xpdfmax[flav] = total;
  }
  return total;
}


///////////////////////////////////////////////////////////////////////////////
//
// Kernels for integration - will yield the norm of sea/valence at Q0^2
//
///////////////////////////////////////////////////////////////////////////////

double Continued_PDF::Sea_Kernel::operator()(double x) {
  if (x<m_xmin || x>m_xmax) return 0.;
  p_pdf->Calculate(x,m_Q02);
  double xpdf(0.);
  for (std::list<Flavour>::iterator flit=p_pdfpartons->begin();
       flit!=p_pdfpartons->end();flit++) {
    if (m_bunch==Flavour(kf_p_plus)) {
      if ((*flit)==Flavour(kf_u) || (*flit)==Flavour(kf_d)) 
	continue;
      else if ((*flit)==Flavour(kf_u).Bar() || (*flit)==Flavour(kf_d).Bar())
	xpdf += 2.*p_pdf->GetXPDF((*flit));
      else 
	xpdf += p_pdf->GetXPDF((*flit));
    }
    else if (m_bunch==Flavour(kf_p_plus).Bar()) {
      if ((*flit)==Flavour(kf_u).Bar() || (*flit)==Flavour(kf_d).Bar()) 
	continue;
      else if ((*flit)==Flavour(kf_u) || (*flit)==Flavour(kf_d))
	xpdf += 2.*p_pdf->GetXPDF((*flit));
      else 
	xpdf += p_pdf->GetXPDF((*flit));
    }
  }
  return xpdf;
}

double Continued_PDF::Valence_Kernel::operator()(double x) {
  if (x<m_xmin || x>m_xmax) return 0.;
  p_pdf->Calculate(x,m_Q02);
  double xpdf(0.);
  for (std::list<Flavour>::iterator flit=p_pdfpartons->begin();
       flit!=p_pdfpartons->end();flit++) {
    if (m_bunch==Flavour(kf_p_plus)) {
      if ((*flit)==Flavour(kf_u) || (*flit)==Flavour(kf_d))
	xpdf += p_pdf->GetXPDF((*flit));
      else if ((*flit)==Flavour(kf_u).Bar() || (*flit)==Flavour(kf_d).Bar())
	xpdf -= p_pdf->GetXPDF((*flit));
    }
    else if (m_bunch==Flavour(kf_p_plus).Bar()) {
      if ((*flit)==Flavour(kf_u) || (*flit)==Flavour(kf_d))
	xpdf += p_pdf->GetXPDF((*flit));
      else if ((*flit)==Flavour(kf_u).Bar() || (*flit)==Flavour(kf_d).Bar())
	xpdf -= p_pdf->GetXPDF((*flit));
    }
  }
  return xpdf;
}

 
