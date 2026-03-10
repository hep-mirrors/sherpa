#include "HADRON_RESCATTERING/XSecs/HeraParam.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/MathTools.H"

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

HeraParam::HeraParam()
{
    InitParams();
    
}
HeraParam::~HeraParam()
{
    //dector init.
}

void HeraParam::InitParams() {
  // parameters are                         a,  b,    n,     c,   d,  
  m_params[HeraParameters::pPiMinus]    = { 0, 11.4, -0.4, 0.079, 0};
  m_params[HeraParameters::pPiPlus]     = { 0, 11.4, -0.4, 0.079, 0};
  m_params[HeraParameters::nPiPlus]     = { 0, 11.4, -0.4, 0.079, 0};
}


double HeraParam::xs_tot(HeraParameters::code tag,const double & s) {
//   if (m_params.find(tag)==m_params.end()) 
//   {
//     THROW(fatal_error,"Tag not found.");
//   }
  double a  = m_params[tag][0];
  double b  = m_params[tag][1];
  double n = m_params[tag][2];
  double c = m_params[tag][3];
  double d = m_params[tag][4];
  double p = sqrt(sqr(s)-4.*m1*s)/(2.*m2);

  double retVal = a + b*pow(p,n) + + c*pow(log(p),2) + d*log(p);
  return (retVal);
}


  