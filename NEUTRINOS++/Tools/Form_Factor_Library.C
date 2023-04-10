#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include <iomanip>

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;

std::ostream & NEUTRINOS::operator<<(std::ostream & s,const ff_type::code & type) {
  if (type==ff_type::none)        s<<setw(12)<<"none";
  if (type==ff_type::dipole)      s<<setw(12)<<"dipole";
  if (type==ff_type::exponential) s<<setw(12)<<"exponential";
  if (type==ff_type::Gaussian)    s<<setw(12)<<"Gaussian";
  if (type==ff_type::unknown)     s<<setw(12)<<"unknown";
  return s;
}

double Dipole_Form_Factor::Calc(const double & q2)      { return m_norm / sqr(1.+m_arg*q2); }

double Exponential_Form_Factor::Calc(const double & q2) { return m_norm * exp(-m_arg*q2); }

double Gaussian_Form_Factor::Calc(const double & q2)    { return m_norm * exp(-sqr(m_arg*q2)); }
