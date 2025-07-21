#include "HADRON_RESCATTERING/XSecs/PiPi.H"
#include "HADRON_RESCATTERING/XSecs/HR_Parameters.H"
#include "ATOOLS/Org/Message.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

PiPi::PiPi() : m_offset(5.)
{
}

double PiPi::XStot(const Flavour & A,const Flavour & B,const double & s) {
  return 0.;
}

double PiPi::XSel(const double & s) {
  return 0.;
}



