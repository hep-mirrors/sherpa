#include "Continued_PDF.H"

#include "MathTools.H"

using namespace PDF;

Continued_PDF::Continued_PDF(PDF_Base *const pdf,const pcs::type scheme):
  p_pdf(pdf),
  m_scheme(scheme),
  m_epsilon(1.e-6)
{
  m_type=std::string("CPDF(")+p_pdf->Type()+std::string(")");
  m_xmin=0.0;
  m_xmax=p_pdf->XMax();
  m_q2min=0.0;
  m_q2max=ATOOLS::sqr(p_pdf->Q2Max());
  m_bunch=p_pdf->Bunch();
  m_partons=p_pdf->Partons();
}

Continued_PDF::~Continued_PDF()
{
  delete p_pdf;
}

void Continued_PDF::Output()
{ 
  return p_pdf->Output(); 
}

double Continued_PDF::ContinueLinear(const ATOOLS::Flavour flavour)
{
  double x=ATOOLS::Max(p_pdf->XMin(),m_x);
  double mu2=ATOOLS::Max(p_pdf->Q2Min(),m_mu2);
  p_pdf->Calculate(x,mu2);
  double xpdf=p_pdf->GetXPDF(flavour);
  return ATOOLS::Min(x/p_pdf->XMin(),mu2/p_pdf->Q2Min())*xpdf;
}

double Continued_PDF::ContinueConstant(const ATOOLS::Flavour flavour)
{
  double x=ATOOLS::Max(p_pdf->XMin(),m_x);
  double mu2=ATOOLS::Max(p_pdf->Q2Min(),m_mu2);
  p_pdf->Calculate(x,mu2);
  return p_pdf->GetXPDF(flavour);
}

void Continued_PDF::Calculate(double x,double mu2)
{
  if (mu2>p_pdf->Q2Min() && x>p_pdf->XMin()) {
    p_pdf->Calculate(x,mu2);
    m_continue=false;
  }
  else {
    m_continue=true;
    m_x=x;
    m_mu2=mu2;
  }
}

double Continued_PDF::GetXPDF(const ATOOLS::Flavour flavour)
{
  if (!m_continue) return p_pdf->GetXPDF(flavour);
  switch (m_scheme) {
  case pcs::constant: return ContinueConstant(flavour);
  case pcs::linear: return ContinueLinear(flavour);
  default: return ContinueLinear(flavour);
  }
  return 0.0;
}

PDF_Base *Continued_PDF::GetCopy()
{ 
  return p_pdf->GetCopy(); 
}

void Continued_PDF::AssignKeys(ATOOLS::Integration_Info *const info)
{
}

bool Continued_PDF::Collinear(const double kp2) const
{
  return true;
}

PDF_Base *Continued_PDF::GetBasicPDF() 
{
  return p_pdf;
}

double Continued_PDF::Cut(const std::string &type) 
{
  return PDF_Base::Cut(type);
}
