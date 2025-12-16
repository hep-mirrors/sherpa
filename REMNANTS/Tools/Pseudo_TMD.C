#include "REMNANTS/Tools/Pseudo_TMD.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace REMNANTS;
using namespace ATOOLS;

Pseudo_TMD::
Pseudo_TMD(PDF::PDF_Base * pdf,const double & Y0, const double & sigma0) :
  p_pdf(pdf), m_bunch(pdf->Bunch()), 
  m_xmin(p_pdf->XMin()), m_xmax(p_pdf->XMax()), m_Q02(p_pdf->Q2Min()),
  m_pt2min(0.), m_pt2max(1.e12), m_Y0(Y0), m_sigma02(sqr(sigma0)),
  m_x0(0.25/(rpa->gen.Ecms()/2.))
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
  Scan();
  //Test();
}

Pseudo_TMD::~Pseudo_TMD() {
  m_pdfpartons.clear();
}

void Pseudo_TMD::Scan()  {
  for (std::list<Flavour>::iterator flit=m_pdfpartons.begin();
       flit!=m_pdfpartons.end();flit++) {
    Flavour flav = (*flit); 
    m_xpdfmax[flav] = 0.;
  }
  for (size_t i=0;i<2000;i++) {
    if (i==0 || i==1000) continue;
    double x = (i<1000?double(i)/1000.:0.001*double(i-1000)/1000.);
    for (size_t beam=0;beam<2;beam++) AllPartons(x,0.,0.,0.);
  }
  for (std::list<ATOOLS::Flavour>::iterator flit=m_pdfpartons.begin();
       flit!=m_pdfpartons.end();flit++) {
    m_xpdfmax[(*flit)] = XPDF((*flit));
  }
}

void Pseudo_TMD::Calculate(const double & x,const double & Q2,
			   const double & pt2,const double & y) {
  m_x   = x; 
  m_Q2  = Q2;
  m_pt2 = pt2;
  m_Y   = y;
  if (Q2<m_Q02) p_pdf->Calculate(m_x,m_Q02);
           else p_pdf->Calculate(m_x,m_Q2);
} 

double Pseudo_TMD::AllPartons(const double & x,const double & Q2,
			      const double & pt2,const double & y) {
  Calculate(x,Q2,pt2,y);
  double val(0.), test(0.);
  for (std::list<ATOOLS::Flavour>::iterator flit=m_pdfpartons.begin();
       flit!=m_pdfpartons.end();flit++) {
    val += XTMD((*flit),true);
  }
  return val;
}

double Pseudo_TMD::XPDF(const Flavour & flav,const bool & defmax) {
  double xpdf   = p_pdf->GetXPDF(flav);
  if (defmax && xpdf>m_xpdfmax[flav]) m_xpdfmax[flav] = xpdf;
  return xpdf;
}

double Pseudo_TMD::XTMD(const Flavour & flav,const bool & defmax) {
  double xpdf   = XPDF(flav,defmax), sigma2 = Sigma2(m_Y);
  //msg_Out()<<"--> "<<METHOD<<"(Y = "<<m_Y<<" --> sigma^2 = "<<sigma2<<", xpdf = "<<xpdf<<")\n";
  return 1./(M_PI*sigma2) * exp(-m_pt2/(sigma2*dabs(xpdf)));
  return xpdf/(M_PI*sigma2) * exp(-m_pt2/sigma2);
}

const double Pseudo_TMD::Sigma2(const double & y) const {
  return m_sigma02;
  return m_sigma02 * sqr(dabs(y-m_Y0));
}

void Pseudo_TMD::Test()  {
  Pseudo_TMD_PT2Test testPT2(this);
  Gauss_Integrator   intPT2(&testPT2);
  double x,Q2,y;
  for (size_t i=1;i<5;i++) {
    for (size_t j=0;j<4;j++) {
      for (size_t k=0;k<6;k++) {
	x  = pow(0.1,i);
	Q2 = pow(10.,j);
	y  = double(k);
	testPT2.SetXQ2(x,Q2);
	testPT2.Sety(y);
	testPT2.SetFlav(Flavour(kf_gluon));
	double sum = intPT2.Integrate(0.,10.,0.0001,2);
	double tru = XPDF(Flavour(kf_gluon));
	msg_Out()<<"Int_pt2 TMD("<<x<<", "<<Q2<<", y = "<<y<<") = "
		 <<sum<<" vs. "<<tru<<", "
		 <<"delta = "<<((sum-tru)*100./tru)<<"%.\n";
      }
    }
  }
  exit(1);
}
  

double Pseudo_TMD_PT2Test::operator()(double pt) {
  //if (m_flav==Flavour(kf_none))
  //return p_pdf->AllPartons(m_x,m_Q2,sqr(pt),m_y);
  p_pdf->Calculate(m_x,m_Q2,sqr(pt),m_y);
  return 2.*M_PI*pt*p_pdf->XTMD(m_flav);
}
