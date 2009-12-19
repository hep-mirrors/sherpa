#include "PDF/Main/ISR_Base.H"

using namespace PDF;

ISR_Base::ISR_Base(PDF_Base *pdf):
  p_pdf(pdf),
  m_exponent(0.),
  m_xmax(1.),
  m_on((bool)pdf)
{
  if (pdf!=NULL) {
    m_exponent=p_pdf->Exponent(); 
    m_xmax=p_pdf->XMax();
  }
}

ISR_Base::~ISR_Base() 
{
  if (p_pdf!=NULL) delete p_pdf;
}
