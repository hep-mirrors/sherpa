//#include <stdio.h>
//#include <math.h>
//#include <fstream> 
#include "Couplings_EW.H"
#include "MathTools.H"
#include "Run_Parameter.H"

using namespace AMEGIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;


Couplings_EW::Couplings_EW()
{
  _CKM    = 0;
}

Couplings_EW::~Couplings_EW()
{
  if (_CKM!=0) {
    for (short int i=0;i<3;i++) delete[] _CKM[i];
    delete[] _CKM;
  }
}

void Couplings_EW::Init()
{
  Data_Read dr(rpa.GetPath()+std::string("/")+rpa.me.ModelFile());

  vev           = dr.GetValue<double>("v");
  sinthetaW_MZ  = dr.GetValue<double>("sinTW^2");
  if (sinthetaW_MZ==NotDefined<double>())
    sinthetaW_MZ  = dr.GetValue<double>("SIN2_TW");

  double lambda = dr.GetValue<double>("lambda");
  double A      = dr.GetValue<double>("A");
  double rho    = dr.GetValue<double>("rho");
  double eta    = dr.GetValue<double>("eta");

  if (_CKM==0) {
    _CKM = new Complex*[3];
    for (short int i=0;i<3;i++) _CKM[i] = new Complex[3];
  }
  _CKM[0][0] = Complex(1.-sqr(lambda)/2.,0.);
  _CKM[0][1] = Complex(lambda,0.);
  _CKM[0][2] = Complex(rho*A*pow(lambda,3.),-eta*A*pow(lambda,3.));
  _CKM[1][0] = Complex(-lambda,0.);
  _CKM[1][1] = Complex(1.-sqr(lambda)/2.,0.);
  _CKM[1][2] = Complex(A*sqr(lambda),0.);
  _CKM[2][0] = Complex((1-rho)*A*pow(lambda,3.),-eta*A*pow(lambda,3.));
  _CKM[2][1] = Complex(-A*sqr(lambda),0.);
  _CKM[2][2] = Complex(1.,0.); 
}

double Couplings_EW::VEV()       {return vev;}
double Couplings_EW::SinThetaW() {return sqrt(sinthetaW_MZ);}
double Couplings_EW::CosThetaW() {return sqrt(1.-sinthetaW_MZ);}
Complex Couplings_EW::CKM(short int i,short int j) {return _CKM[i][j];}












