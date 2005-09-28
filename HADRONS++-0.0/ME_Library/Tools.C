#include "Tools.H"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  general tools  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace HADRONS;
using namespace ATOOLS;

double Tools::Lambda( double a, double b, double c )
{
  double L = sqr(a-b-c)-4.*b*c;
  if (L>0.) return L;
  return 0.;
}

Complex Tools::BreitWigner( double s, double Mass2, double MassWidth )
{
  return Mass2/Complex(Mass2-s,-1.*MassWidth );
}
Complex Tools::BreitWigner( double s, double Mass2, double Width, double ms, double lambda )
{
  double MassWidth = OffShellMassWidth( s, Mass2, Width, ms, lambda );
  return BreitWigner( s, Mass2, MassWidth );
}

Complex Tools::BreitWigner( double s, double Mass2, double Width, double ms1, double ms2, double lambda )
{
  double MassWidth = OffShellMassWidth( s, Mass2, Width, ms1, ms2, lambda );
  return BreitWigner( s, Mass2, MassWidth );
}

double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms, double lambda )
{
  if (s>4.*ms) {
	if (Mass2>4.*ms)
	  return( sqrt(s)*Width*pow(Mass2/s,lambda) * pow( (s-4.*ms)/(Mass2-4.*ms), 1.5 ) );
	return sqrt(Mass2)*Width;
  }
  return 0.;	
}

double Tools::OffShellMassWidth( double s, double Mass2, double Width, double ms1, double ms2, double lambda )
{
  double threshold = ms1+ms2+2.*sqrt(ms1*ms2);
  if (s>threshold) {
	if (Mass2>threshold)
	  return( sqrt(s)*Width*pow(Mass2/s,lambda) * pow( Mass2/s*Lambda(s,ms1,ms2)/Lambda(Mass2,ms1,ms2), 1.5 ) );
	return sqrt(Mass2)*Width;
  }
  return 0.;
}

