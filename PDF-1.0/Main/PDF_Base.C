#include "PDF_Base.H"

#include "Info_Key.H"

using namespace PDF;

PDF_Base::PDF_Base():
  m_type("none"),
  m_exponent(1.) {}

PDF_Base::~PDF_Base() {}

bool PDF_Base::Collinear(const double kp2) const
{
  return true;
}

PDF_Base *PDF_Base::GetBasicPDF() 
{
  return this;
}

double PDF_Base::Cut(const std::string &type)
{
  return ATOOLS::UNDEFINED_LOWER;
}
