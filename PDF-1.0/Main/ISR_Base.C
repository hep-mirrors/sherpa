#include "ISR_Base.H"

using namespace PDF;

ISR_Base::ISR_Base(PDF_Base *_p_pdf):
  p_pdf(_p_pdf),
  m_on((bool)_p_pdf), 
  m_kmr(_p_pdf->Type()==std::string("DUPDF")) {}

ISR_Base::~ISR_Base() {}
