#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"


using namespace HADRON_RESCATTERING;
using namespace ATOOLS;
using namespace std;

HPR1R2::HPR1R2() :
  m_eta1(-0.4473), m_eta2(-0.5486)
{
  InitParams();
}

HPR1R2::~HPR1R2() {}

double HPR1R2::xs_tot(hpr1r2::code tag,const double & s) {
  if (m_params.find(tag)==m_params.end()) {
    THROW(fatal_error,"Tag not found.");
  }
  double H  = m_params[tag][0];
  double P  = m_params[tag][1];
  double R1 = m_params[tag][2];
  double R2 = m_params[tag][3];
  double s0 = m_params[tag][4];
  return ( P + H*sqr(log(s/s0)) +
	   R1*pow(s/s0,m_eta1)  +
	   R2*pow(s/s0,m_eta2) );
}

double HPR1R2::xs_el(hpr1r2::code tag,const double & s) {
  if (m_params.find(tag)==m_params.end()) {
    THROW(fatal_error,"Tag not found.");
  }
  double H  = m_params[tag][0];
  double P  = m_params[tag][1];
  double R1 = m_params[tag][2];
  double R2 = m_params[tag][3];
  double s0 = m_params[tag][4];
  double B  = m_params[tag][5];
  double xstot = ( P + H*sqr(log(s/s0)) +
		   R1*pow(s/s0,m_eta1)  +
		   R2*pow(s/s0,m_eta2) );
  double rho   = ( M_PI*H*log(s/s0)                        -
		   R1*pow(s/s0,m_eta1)*tan(m_eta1*M_PI/2.) +
		   ( R2*pow(s/s0,m_eta2)*
		     cos(m_eta2*M_PI/2.)/sin(m_eta2*M_PI/2.) ) )/ xstot;
  return sqr(xstot) * (1+sqr(rho)) / (16.*M_PI*B); 
}

void HPR1R2::InitParams() {
  // parameters are               H,     P,    R1,     R2,    s0,  B
  m_params[hpr1r2::pp]    = { 0.272, 34.41, 13.07, -7.394, 15.98, 13. };
  m_params[hpr1r2::ppbar] = { 0.272, 34.41, 13.07,  7.394, 15.98, 13. };
  m_params[hpr1r2::pn]    = { 0.272, 34.71, 12.52,  6.660, 15.98, 13. };
  //newly added
  // parameters are               H,     P,    R1,     R2,    s0,  B

  m_params[hpr1r2::pPiMinus]    = { 0.272, 18.75, 9.56,  1.767, 15.98, 13. };
  m_params[hpr1r2::pPiPlus]    = { 0.272, 18.75, 9.56,  -1.767, 15.98, 13. };
  m_params[hpr1r2::nPiMinus]    = { 0.272, 18.75, 9.56,  1.767, 15.98, 13. };
  m_params[hpr1r2::nPiPlus]    = { 0.272, 18.75, 9.56,  -1.767, 15.98, 13. };
  m_params[hpr1r2::pKMinus]    = { 0.272, 16.36, 4.29,  3.408, 15.98, 13. };
  m_params[hpr1r2::pKBarZero]    = { 0.272, 16.36, 4.29, 3.408, 15.98, 13. };
  m_params[hpr1r2::nKMinus]    = { 0.272, 16.31, 3.70,  1.826, 15.98, 13. };
  m_params[hpr1r2::nKBarZero]    = { 0.272, 16.31, 3.70,  1.826, 15.98, 13. };


  
}

  

  
