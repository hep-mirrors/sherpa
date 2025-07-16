#include "PDF/Main/ISR_Base.H"

using namespace PDF;

ISR_Base::ISR_Base(PDF_Base *pdf)
    : p_pdf(pdf), m_type(isrtype::unknown), m_weight(0.), m_exponent(0.),
      m_xmax(1.), m_on((bool)pdf) {
  if (pdf != nullptr) {
    m_exponent = p_pdf->Exponent();
    m_xmax = p_pdf->XMax();
  }
}

ISR_Base::~ISR_Base() { delete p_pdf; }

std::ostream &PDF::operator<<(std::ostream &s, const PDF::isrtype::code type) {
  switch (type) {
  case isrtype::intact:
    s << "intact";
    break;
  case isrtype::lepton:
    s << "lepton";
    break;
  case isrtype::hadron:
    s << "hadron";
    break;
  case isrtype::unknown:
  default:
    s << "unknown";
    break;
  }
  return s;
}
