#include "Continued_PDF.H"

using namespace PDF;

Continued_PDF::Continued_PDF(PDF_Base *const pdf,const pcs::type scheme):
  p_pdf(pdf),
  m_scheme(scheme) {}

void Continued_PDF::Output()
{ 
  return p_pdf->Output(); 
}

void Continued_PDF::Calculate(double x,double z,double kperp2,double mu2)
{
  std::cout<<"continue";
  p_pdf->Calculate(x,z,kperp2,mu2);
}

double Continued_PDF::GetXPDF(const ATOOLS::Flavour flavour)
{
  std::cout<<"continued"<<std::endl;
  return p_pdf->GetXPDF(flavour);
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
