#include "CSSHOWER++/Sudakov/CS_SingleSudakov_Base.H"

using namespace CSSHOWER;
using namespace std;

const double CS_SingleSudakov_Base::s_Nc = 3.;
const double CS_SingleSudakov_Base::s_CF = (s_Nc*s_Nc-1.)/(2.*s_Nc);
const double CS_SingleSudakov_Base::s_CA = s_Nc;
const double CS_SingleSudakov_Base::s_TR = 1./2.;
const double CS_SingleSudakov_Base::s_Tf = s_TR*3.0;

double  CS_SingleSudakov_Base::Lambda(const double a,const double b,const double c) {
  return a*a+b*b+c*c-2.*(a*b + a*c + b*c);
}
