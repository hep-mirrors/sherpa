#include "Message.H"
#include "Tools.H"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  general tools  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

using namespace HADRONS;
using namespace ATOOLS;

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

// Breit Wigner with running width (2 particle final state with same mass) 
Complex Tools::BreitWigner( double s, double Mass2, double Width, double ms )
{
  double MassWidth = OffShellMassWidth( s, Mass2, Width, ms );
  return BreitWigner( s, Mass2, MassWidth );
}

// Breit Wigner with running width (2 particle final state with differenent masses)
Complex Tools::BreitWigner( double s, double Mass2, double Width, double ms1, double ms2 )
{
  double MassWidth = OffShellMassWidth( s, Mass2, Width, ms1, ms2 );
  return BreitWigner( s, Mass2, MassWidth );
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


Vec4D Tools::RealBosonPolarizationVector( Vec4D p, int lambda, double M2, bool & iszero )
{
  iszero = false;
  Vec4D eps;
  double ct (p.CosTheta()), st (p.SinTheta()), cp (p.CosPhi()), sp (p.SinPhi());
  if (IsEqual(p.PSpat(), 0.)) {      // decay from rest
    ct = 0.; st = 1.;
    cp = 1.; sp = 0.;
  }
  PRINT_INFO(ct<<" "<<st<<"    "<<cp<<" "<<sp);
  switch( lambda ) {
    case 1  : if ( !IsEqual(p.PSpat(),0.) ) 
                eps = Vec4D( p.PSpat(), p[0]/p.PSpat()*Vec3D(p) );
              else {
                eps = Vec4D( 0.,0.,0.,0. );
                iszero = true;
              }
              break;
    case 2  : eps = Vec4D( 0., ct*cp, ct*sp, -1.*st );
              break;
    case 3  : eps = Vec4D( 0., -1.*sp, cp, 0. );
              break;
    default : if( !IsEqual(p.Abs2(),M2) ) {
                eps = Vec4D( p );
                eps *= sqrt( (p.Abs2()-M2)/(p.Abs2()*M2) );
              }
              else {
                eps = Vec4D( 0.,0.,0.,0. );
              }
              break;
  }
  return eps;
}


void Tools::ComplexBosonPolarizationVector(  Vec4D p, int lambda, double M2, Vec4D * eps)
{
  double ct (p.CosTheta()), st (p.SinTheta()), cp (p.CosPhi()), sp (p.SinPhi());
  if (IsEqual(p.PSpat(), 0.)) {      // decay from rest
    ct = 0.; st = 1.;
    cp = 1.; sp = 0.;
  }
  switch( lambda ) {
    case 0: 
      eps[0] = 1.0/p.Mass() * Vec4D(p.PSpat(),p[0]*Vec3D(p)/p.PSpat());
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
    case 1:
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,-sp,cp,0.0);
      break;
    case 2:
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,sp,-cp,0.0);
      break; 
    case 3:
      eps[0] = sqrt( (p.Abs2()-M2)/(p.Abs2()*M2) ) * Vec4D(p);
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
  }
}


void Tools::ComplexBosonPolarizationVector( Vec4D p, int lambda, Vec4D * eps )
{
  double ct (p.CosTheta()), st (p.SinTheta()), cp (p.CosPhi()), sp (p.SinPhi());
  if (IsEqual(p.PSpat(), 0.)) {      // decay from rest
    ct = 0.; st = 1.;
    cp = 1.; sp = 0.;
  }
  switch( lambda ) {
    case 0: 
      eps[0] = 1.0/p.Mass() * Vec4D(p.PSpat(),p[0]*Vec3D(p)/p.PSpat());
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
    case 1:
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,-sp,cp,0.0);
      break;
    case 2:
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,sp,-cp,0.0);
      break; 
    case 3:
      eps[0] = Vec4D(0.0,0.0,0.0,0.0);
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
  }
}
