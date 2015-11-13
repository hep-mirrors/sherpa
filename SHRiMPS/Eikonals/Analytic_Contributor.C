#include "SHRiMPS/Eikonals/Analytic_Contributor.H"
#include "ATOOLS/Math/MathTools.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Analytic_Contributor::
operator()(const double & b,const double & y) const {
  if (y<-m_Y || y>m_Y || b>p_ff->Bmax()) return 0.;
  return p_ff->FourierTransform(b)*exp(m_Delta*(m_Y+m_sign*y));
}

