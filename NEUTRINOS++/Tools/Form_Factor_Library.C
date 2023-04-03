#include "NEUTRINOS++/Tools/Form_Factor_Library.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace NEUTRINOS;
using namespace ATOOLS;
using namespace std;


double Dipole_Form_Factor::Calc(const double & q2)      { return m_norm / sqr(1.+m_arg*q2); }

double Exponential_Form_Factor::Calc(const double & q2) { return m_norm * exp(-m_arg*q2); }

double Gaussian_Form_Factor::Calc(const double & q2)    { return m_norm * exp(-sqr(m_arg*q2)); }
