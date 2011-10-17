#include "HADRONS++/Main/Tools.H"
#include "ATOOLS/Org/Message.H"


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  general tools  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace HADRONS;
using namespace ATOOLS;

const double Tools::Vud = 0.97377;
const double Tools::Vus = 0.2257;
const double Tools::Vub = 4.31e-3;
const double Tools::Vcd = 0.230;
const double Tools::Vcs = 0.957;
const double Tools::Vcb = 41.6e-3;
const double Tools::Vtd = 7.4e-3;
const double Tools::Vts = Tools::Vtd/0.208;
const double Tools::Vtb = 1.0;

// 3 particle phase space function lambda
double Tools::Lambda( double a, double b, double c )
{
  double L = sqr(a-b-c)-4.*b*c;
  if (L>0.) return L;
  return 0.;
}

// standard Breit Wigner with given Mass * Width
Complex Tools::BreitWigner( double s, double Mass2, double MassWidth )
{
  return Mass2/Complex(Mass2-s,-1.*MassWidth );
}

// standard Breit Wigner with given Mass * Width
Complex Tools::BreitWignerFix( double s, double Mass2, double MassWidth )
{
  return Complex(Mass2,-1.*MassWidth)/Complex(Mass2-s,-1.*MassWidth );
}

// off shell mass * width (2 particle final state with same mass)
double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms )
{
  if (s>4.*ms && Mass2>4.*ms)
    return( sqrt(s)*Width*Mass2/s * pow( (s-4.*ms)/(Mass2-4.*ms), 1.5 ) );
  return 0.;	
}

// off shell mass * width (2 particle final state with different masses)
double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms1, double ms2 )
{
  double threshold = ms1+ms2+2.*sqrt(ms1*ms2);
  if (Mass2>threshold && s>threshold)
	  return( sqrt(s)*Width*Mass2/s * pow( Mass2/s*Lambda(s,ms1,ms2)/Lambda(Mass2,ms1,ms2), 1.5 ) );
  return 0;
}
