#include "PDF/Main/Continued_PDF.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PDF;
using namespace ATOOLS;

Continued_PDF::Continued_PDF(PDF_Base * pdf) :
  p_pdf(pdf), m_bunch(p_pdf->Bunch()), 
  m_xmin(p_pdf->XMin()), m_xmax(p_pdf->XMax()), m_Q02(p_pdf->Q2Min()),
  m_geta(1.), m_glambda(0.25)
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
    Flavour flav = (*flit); 
    m_xpdfmax[flav] = 0.; m_xmaxpdf[(*flit)] = 0.;
  }
  m_x0 = 0.25/(rpa->gen.Ecms()/2.);
  CalculateNorms();
  Scan();
  //Test();
}

Continued_PDF::~Continued_PDF() {
  m_pdfpartons.clear();
}

void Continued_PDF::CalculateNorms() {
  Sea_Kernel sea(p_pdf,m_bunch,&m_pdfpartons,m_Q02);
  Gauss_Integrator sintegrator(&sea);
  m_Snorm = sintegrator.Integrate(m_xmin,m_xmax,0.0001,1);
  Valence_Kernel val(p_pdf,m_bunch,&m_pdfpartons,m_Q02);
  Gauss_Integrator vintegrator(&val);
  m_Vnorm = vintegrator.Integrate(m_xmin,m_xmax,0.0001,1);
  // this assumes xg(x) = (1-x)^eta x^lambda
  m_gnorm = 
    exp(Gammln(m_geta+1.))*exp(Gammln(m_glambda+1.))/
    exp(Gammln(m_geta+m_glambda+2.));
  // this assumes xg(x) = (x+x_0)^lambda
  //m_gnorm = (pow(1.+m_x0,1.+m_glambda)-pow(m_x0,1.+m_glambda))/(1.+m_glambda);  
}


void Continued_PDF::Test()  {
  Continued_PDF_Test test(this);
  msg_Out()<<METHOD<<": sea = "<<m_Snorm<<", valence = "<<m_Vnorm<<", gluons = "<<m_gnorm<<".\n"
	   <<"   Q_0^2 for the transition is given by Q_0^2 = "<<m_Q02<<".\n";
  Gauss_Integrator testint(&test);
  for (size_t i=0;i<20;i++) {
    test.SetQ2(double(i)/5.);
    double sum = testint.Integrate(0.,1.,0.000001,1);
    msg_Out()<<"   --> sum at Q2 = "<<m_Q2<<" = "<<sum<<"\n";
  }
  exit(1);
}
  
void Continued_PDF::Scan()  {
  for (size_t i=0;i<2000;i++) {
    if (i==0 || i==1000) continue;
    double x = (i<1000?double(i)/1000.:0.001*double(i-1000)/1000.);
    for (size_t beam=0;beam<2;beam++) AllPartons(x,0.);
  }
  Calculate(m_xmaxpdf[Flavour(kf_gluon)],0.);
  msg_Out()<<METHOD<<" yields fmax(xmax = "<<m_xmaxpdf[Flavour(kf_gluon)]<<") = "
  	   <<m_xpdfmax[Flavour(kf_gluon)]<<"\n";
}

void Continued_PDF::Calculate(const double & x,const double & Q2) {
  m_x  = x; 
  m_Q2 = Q2;
  if (Q2<m_Q02) p_pdf->Calculate(m_x,m_Q02);
           else p_pdf->Calculate(m_x,m_Q2);
} 

double Continued_PDF::AllPartons(const double & x,const double & Q2) {
  Calculate(x,Q2);
  double val(0.), test(0.);
  for (std::list<ATOOLS::Flavour>::iterator flit=m_pdfpartons.begin();
       flit!=m_pdfpartons.end();flit++) {
    val += XPDF((*flit),true);
  }
  return val;
}

double Continued_PDF::XPDF(const Flavour & flav,const bool & defmax) {
  if (m_Q2>m_Q02) {
    if (!flav.IsDiQuark())                                       return p_pdf->GetXPDF(flav);
    if (m_bunch==Flavour(kf_p_plus) && !flav.IsAnti())           return XPDF(Flavour(kf_u)); 
    else if (m_bunch==Flavour(kf_p_plus).Bar() && flav.IsAnti()) return XPDF(Flavour(kf_u).Bar()); 
  }
  double seapart(0.), valpart(0.), total(0.);
  if (flav==Flavour(kf_gluon)) {
    seapart = p_pdf->GetXPDF(flav);
    valpart = GluonAtZero(m_x);
    total = seapart * m_Q2/m_Q02 + valpart * (1.-m_Q2/m_Q02);
  }
  else {
    if (m_bunch==Flavour(kf_p_plus)) {
      if (flav==Flavour(kf_u) || flav==Flavour(kf_d)) {
	seapart = p_pdf->GetXPDF(flav.Bar());
	valpart = p_pdf->GetXPDF(flav)-p_pdf->GetXPDF(flav.Bar());
      }
      else
	seapart = p_pdf->GetXPDF(flav);
    }
    else if (m_bunch==Flavour(kf_p_plus).Bar()) {
      if (flav==Flavour(kf_u).Bar() || flav==Flavour(kf_d).Bar()) {
	seapart = p_pdf->GetXPDF(flav.Bar());
	valpart = p_pdf->GetXPDF(flav)-p_pdf->GetXPDF(flav.Bar());
      }
      else
	seapart = p_pdf->GetXPDF(flav);
    }
    total = seapart * m_Q2/m_Q02 + valpart * (1.+m_Snorm*(1.-m_Q2/m_Q02));
  }
  if (defmax && total>m_xpdfmax[flav]) {
    m_xmaxpdf[flav] = m_x;
    m_xpdfmax[flav] = total;
  }
  return total;
}

double Continued_PDF::GluonAtZero(const double & x) {
  // this assumes xg(x) = (1-x)^eta x^lambda
  return pow(1-x,m_geta) * + pow(x,m_glambda) * (1.-m_Vnorm)/m_gnorm;
  // this assumes xg(x) = (x+x_0)^lambda
  //return pow(x+m_x0,m_glambda) * (1.-m_Vnorm)/m_gnorm;
}

double Continued_PDF_Test::operator()(double x) {
  p_pdf->AllPartons(x,m_Q2);
  return (p_pdf->XPDF(Flavour(kf_gluon)) +
	  p_pdf->XPDF(Flavour(kf_u)) + p_pdf->XPDF(Flavour(kf_u).Bar()) +
	  p_pdf->XPDF(Flavour(kf_d)) + p_pdf->XPDF(Flavour(kf_d).Bar()) +
	  p_pdf->XPDF(Flavour(kf_s)) + p_pdf->XPDF(Flavour(kf_s).Bar()) +
	  p_pdf->XPDF(Flavour(kf_d)) + p_pdf->XPDF(Flavour(kf_c).Bar()) +
	  p_pdf->XPDF(Flavour(kf_b)) + p_pdf->XPDF(Flavour(kf_b).Bar()));
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

