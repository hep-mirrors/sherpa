#include "Parton.H"

using namespace CS_SHOWER;
using namespace std;

std::ostream& CS_SHOWER::operator<<(std::ostream& str, const Parton &part) {
  str<<"  Parton : "<<part.m_flav<<" : "<<part.m_mom<<endl;
  if (part.m_pst==pst::IS)      str<<"     (Initial state parton)";
  else if (part.m_pst==pst::FS) str<<"     (Final state parton)  ";
  else                     str<<"                           ";
  if (part.m_cpartner==1)       str<<"  Colour partner: Next."<<endl;
  else if (part.m_cpartner==-1) str<<"  Colour partner: Prev."<<endl;
  else str<<endl;
  return str;
}
