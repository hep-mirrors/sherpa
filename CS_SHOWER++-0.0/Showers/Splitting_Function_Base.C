#include "Splitting_Function_Base.H"

using namespace CS_SHOWER;
using namespace std;

const double Splitting_Function_Base::s_Nc = 3.;
const double Splitting_Function_Base::s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double Splitting_Function_Base::s_CA = s_Nc;
const double Splitting_Function_Base::s_TR = 1./2.;
const double Splitting_Function_Base::s_Tf = s_TR*3.0;

ostream& CS_SHOWER::operator<<(std::ostream& str, const Splitting_Function_Base &base) {
  str<<"  "<<base.m_flavs[0]<<" -> "<<base.m_flavs[1]<<" + "<<base.m_flavs[2]
     <<" : "<<base.m_lastint<<endl;
  return str;
}
