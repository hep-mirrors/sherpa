#include "PHASIC++/Channels/Transverse_Kinematics.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace ATOOLS;

LN_Pair PHASIC::GetLN
(const Vec4D &pi,const Vec4D &pk,const int mode)
{
  double mi2(pi.Abs2()), mk2(pk.Abs2());
  double eps(pi*pk), kap(eps*eps-mi2*mk2);
  if (kap<0.0) return LN_Pair();
  kap=Sign(eps)*sqrt(kap);
  Vec4D l(((eps+kap)*pi-mi2*pk)/(2.0*kap));
  Vec4D n(((eps+kap)*pk-mk2*pi)/(2.0*kap));
  return LN_Pair(l,n,mode);
}

Vec4D PHASIC::LT(const Vec4D &a,const Vec4D &b,const Vec4D &c)
{
  double t(a[1]*b[2]*c[3]+a[2]*b[3]*c[1]+a[3]*b[1]*c[2]
	   -a[1]*b[3]*c[2]-a[3]*b[2]*c[1]-a[2]*b[1]*c[3]);
  double x(-a[0]*b[2]*c[3]-a[2]*b[3]*c[0]-a[3]*b[0]*c[2]
	   +a[0]*b[3]*c[2]+a[3]*b[2]*c[0]+a[2]*b[0]*c[3]);
  double y(-a[1]*b[0]*c[3]-a[0]*b[3]*c[1]-a[3]*b[1]*c[0]
	   +a[1]*b[3]*c[0]+a[3]*b[0]*c[1]+a[0]*b[1]*c[3]);
  double z(-a[1]*b[2]*c[0]-a[2]*b[0]*c[1]-a[0]*b[1]*c[2]
	   +a[1]*b[0]*c[2]+a[0]*b[2]*c[1]+a[2]*b[1]*c[0]);
  return Vec4D(t,-x,-y,-z);
}

double PHASIC::ComputePhi(Vec4D pijt,Vec4D pkt,Vec4D pi)
{
  Vec4D n_perp(0.0,cross(Vec3D(pijt),Vec3D(pkt)));
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_IODebugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    Poincare zrot(pijt,Vec4D::ZVEC);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(LT(pijt,pkt,n_perp));
  l_perp*=1.0/sqrt(dabs(l_perp.Abs2()));
  double cp(-pi*n_perp), sp(-pi*l_perp), phi(atan(sp/cp));
  return cp<0.0?phi+M_PI:(sp>0.0?phi:phi+2.0*M_PI);
}
