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

// standard Breit Wigner with given Mass * Width
Complex Tools::BreitWignerFix( double s, double Mass2, double MassWidth )
{
  return Complex(Mass2,-1.*MassWidth)/Complex(Mass2-s,-1.*MassWidth );
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
//   PRINT_INFO(ct<<" "<<st<<"    "<<cp<<" "<<sp);
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


ATOOLS::Vec4D* Tools::ComplexBosonPolarizationVector( Vec4D p, int lambda, double M2 )
{
  ATOOLS::Vec4D* eps = new ATOOLS::Vec4D[2];
  double ct (p.CosTheta()), st (p.SinTheta()), cp (p.CosPhi()), sp (p.SinPhi());
  if (IsEqual(p.PSpat(), 0.)) {      // decay from rest
    ct = 0.; st = 1.;
    cp = 1.; sp = 0.;
  }
  switch( lambda ) {
    case 0: // +
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,-sp,cp,0.0);
      break;
      case 1: // 0
      eps[0] = 1.0/p.Mass() * Vec4D(p.PSpat(),p[0]*Vec3D(p)/p.PSpat());
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
      case 2: // -
      eps[0] = 1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st);
      eps[1] = 1.0/sqrt(2.0) * Vec4D(0.0,sp,-cp,0.0);
      break; 
      case 3: // s
      eps[0] = sqrt( (p.Abs2()-M2)/(p.Abs2()*M2) ) * Vec4D(p);
      eps[1] = Vec4D(0.0,0.0,0.0,0.0);
      break;
  }
  return eps;
}


ATOOLS::Vec4D* Tools::ComplexBosonPolarizationVector( Vec4D p, int lambda )
{
  double M2 = p.Abs2();
  ATOOLS::Vec4D* eps = Tools::ComplexBosonPolarizationVector(p,lambda,M2);
  return eps;
}


ATOOLS::ComplexVec4D Tools::ComplexBosonPolarizationVectorC( Vec4D p, int lambda, double M2 )
{
  ATOOLS::ComplexVec4D eps;
  double ct (p.CosTheta()), st (p.SinTheta()), cp (p.CosPhi()), sp (p.SinPhi());
  if (IsEqual(p.PSpat(), 0.)) {      // decay from rest
    ct = 0.; st = 1.;
    cp = 1.; sp = 0.;
  }
  if( p.Abs2()<0.0 ) {
    msg.Error()<<METHOD<<": p^2 < 0. Please check definition of "
        <<"complex boson polarization vector for correctness (sign factors) "
        <<"and remove this warning."<<std::endl;
  }
  switch( lambda ) {
    case 0: // +
      return ComplexVec4D(  1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st),
                            1.0/sqrt(2.0) * Vec4D(0.0,-sp,cp,0.0) );
    case 1: // 0
      return ComplexVec4D(  1.0/p.Mass() * Vec4D(p.PSpat(),(p[0]/p.PSpat())*Vec3D(p)),
                            Vec4D(0.0,0.0,0.0,0.0) );
    case 2: // -
      return ComplexVec4D(  1.0/sqrt(2.0) * Vec4D(0.0,ct*cp,ct*sp,-st),
                            1.0/sqrt(2.0) * Vec4D(0.0,sp,-cp,0.0) );
    case 3: // s
      if( p.Abs2()<M2 ) {
        msg.Error()<<METHOD<<": p^2 < m^2. Please check definition of "
            <<"complex boson polarization vector for correctness (sign factors) "
            <<"and remove this warning."<<std::endl;
      }
      return ComplexVec4D(  sqrt( (p.Abs2()-M2)/(p.Abs2()*M2) ) * Vec4D(p),
                            Vec4D(0.0,0.0,0.0,0.0) );
  }
}


ATOOLS::ComplexVec4D Tools::ComplexBosonPolarizationVectorC( Vec4D p, int lambda )
{
  double M2 = p.Abs2();
  return  ComplexBosonPolarizationVectorC(p,lambda,M2);
}


ATOOLS::CMatrix Tools::ComplexSpin2BosonPolarizationVectorC( Vec4D p, int lambda, double M2 )
{
  CMatrix tensor = CMatrix(4);
  ComplexVec4D eps_plus = ComplexBosonPolarizationVectorC( p, 0, M2);
  ComplexVec4D eps_null = ComplexBosonPolarizationVectorC( p, 1, M2);
  ComplexVec4D eps_minus= ComplexBosonPolarizationVectorC( p, 2, M2);
  switch(lambda) {
    case 0:
      for(int mu = 0; mu < 4; mu++){ 
          for(int nu = 0; nu < 4; nu++){      
            Complex eps_munu = (eps_plus[mu]*eps_plus[nu]);    
            tensor[mu][nu] = eps_munu;
          }
      }
      break;
    case 1:
      for(int mu = 0; mu < 4; mu++){ 
          for(int nu = 0; nu < 4; nu++){      
            Complex eps_munu = (eps_plus[mu]*eps_null[nu]+eps_null[mu]*eps_plus[nu])/sqrt(2.0);    
            tensor[mu][nu] = eps_munu;
          }
      }
      break;
    case 2:
      for(int mu = 0; mu < 4; mu++){ 
          for(int nu = 0; nu < 4; nu++){      
            Complex eps_munu = (eps_plus[mu]*eps_minus[nu]+eps_minus[mu]*eps_plus[nu]
				-2.0*eps_null[mu]*eps_null[nu])/sqrt(6.0);    
            tensor[mu][nu] = eps_munu;
          }
      }
      break;
    case 3:
      for(int mu = 0; mu < 4; mu++){ 
          for(int nu = 0; nu < 4; nu++){      
            Complex eps_munu = (eps_minus[mu]*eps_null[nu]+eps_null[mu]*eps_minus[nu])/sqrt(2.0);    
            tensor[mu][nu] = eps_munu;
          }
      }
      break;
    case 4:
      for(int mu = 0; mu < 4; mu++){ 
          for(int nu = 0; nu < 4; nu++){      
            Complex eps_munu = (eps_minus[mu]*eps_minus[nu]);    
            tensor[mu][nu] = eps_munu;
          }
      }
      break;
    default:
      msg.Error()<<METHOD<<" Sorry, wrong index in ComplexSpin2BosonPolarizationVectorC"
      <<". Aborting."<<std::endl;
      abort();
  }
  return tensor;
}


ATOOLS::CMatrix Tools::ComplexSpin2BosonPolarizationVectorC( Vec4D p, int lambda )
{
  double M2 = p.Abs2();
  return  ComplexSpin2BosonPolarizationVectorC( p, lambda, M2);
}



Vec4D Tools::Cross( Vec4D a, Vec4D b, Vec4D c )
{
  return cross(a,b,c);
}
