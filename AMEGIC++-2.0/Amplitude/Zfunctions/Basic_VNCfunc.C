#include "Basic_Func.H"
#include "String_Generator.H"
#include "Run_Parameter.H"
#include "Matrix.H"

using namespace AMEGIC;
using namespace AMATOOLS;

Kabbala Basic_VNCfunc::VNC(const int& z)
{
  int a = iabs(arg[4*z].numb);
  int b = iabs(arg[4*z+2].numb);

  //Called from Z-function
  int sign = 1;

  if (a<BS->GetNmomenta()) sign *= BS->Sign(a);
  if (b<BS->GetNmomenta()) sign *= BS->Sign(b);

  return (sign>0) ?
    sgen->Get_VNCnumber(a,b,VNCcalc(a,b))
    :
    -sgen->Get_VNCnumber(a,b,VNCcalc(a,b));
}

Kabbala Basic_VNCfunc::VNC(const int& a,const int &b)
{
  //Called from Multiple photon-function
  int sign = Sign(a)*Sign(b)*ps[iabs(a)].direction*ps[iabs(b)].direction; 

  if (ps[iabs(a)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(a)].numb);
  if (ps[iabs(b)].numb<BS->GetNmomenta()) sign *= BS->Sign(ps[iabs(b)].numb);

  return (sign>0) ?
    sgen->Get_VNCnumber(ps[iabs(a)].numb,ps[iabs(b)].numb,VNCcalc(ps[iabs(a)].numb,ps[iabs(b)].numb))
    :
    -sgen->Get_VNCnumber(ps[iabs(a)].numb,ps[iabs(b)].numb,VNCcalc(ps[iabs(a)].numb,ps[iabs(b)].numb));
}

Complex Basic_VNCfunc::VNCcalc(const int& a,const int &b)
{ 
  /*
  if ((a==2) && (b==3)) return Complex(0.,0.);
  if ((a==3) && (b==4)) return Complex(0.,0.);
  */

  //cout<<"VNC: "<<a<<":"<<b<<endl;
  Matrix<4>* Theta = AORGTOOLS::rpa.me.GetTheta();
  return (Vec4D(BS->Momentum(a)[0],-1.*Vec3D(BS->Momentum(a)))*((*Theta)*BS->Momentum(b)))/2.;
}
  


















