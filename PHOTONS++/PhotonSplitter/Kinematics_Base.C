#include "PHOTONS++/PhotonSplitter/Kinematics_Base.H"

#include "PHOTONS++/PhotonSplitter/Sudakov.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHOTONS;
using namespace PHASIC;
using namespace ATOOLS;

double Kinematics_FF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double mi2,const double mj2,const double mk2) const
{
  // Q2 = |(pi+pj+pk)|^2
  double pipj=(Q2-mi2-mj2-mk2)*y;
  return pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2; // transverse momentum ordering 
}

double Kinematics_FF::GetY(const double Q2,const double kt2,const double z,
			   const double mi2,const double mj2,const double mk2) const
{
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+mk2) return -1.0;
  return (kt2/(z*(1.0-z))+(1.0-z)/z*mi2+z/(1.0-z)*mj2)/(Q2-mi2-mj2-mk2);
}

double Kinematics_FF::GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double mk2) const
{
  return (q2-mi2-mj2)/(Q2-mi2-mj2-mk2);
}

double Kinematics_FF::GetVirt(const double &Q2,const double &y,const double &z,
           const double mi2,const double mj2,const double mk2) const
{
  // this is the virtuality q2, for t=qbar2=q2-m2 subtract off mij2
  double pipj=(Q2-mi2-mj2-mk2)*y;
  return pipj + mi2 + mj2;
}

bool Kinematics_FF::MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double mk2, const double mij2)
{
  Vec4D Q(pij+pk), rpij(pij);
  Vec4D n_perp(0.0,cross(Vec3D(pij),Vec3D(pk)));
  Poincare cms(pij+pk);
  cms.Boost(rpij);
  if (n_perp.PSpat2()<=rpa->gen.SqrtAccu()) {
    msg_Debugging()<<"Set fixed n_perp\n";
    n_perp=Vec4D(0.0,1.0,1.0,0.0);
    Poincare zrot(rpij,Vec4D::ZVEC);
    zrot.RotateBack(n_perp);
  }
  n_perp*=1.0/n_perp.PSpat();
  Vec4D l_perp(0.0,cross(Vec3D(rpij),Vec3D(n_perp)));
  l_perp*=1.0/l_perp.PSpat();
  double Q2(Q.Abs2()), sij(y*(Q2-mk2)+(1.0-y)*(mi2+mj2));
  double po(sqr(Q2-mij2-mk2)-4.0*mij2*mk2);
  double pn(sqr(Q2-sij-mk2)-4.0*sij*mk2);
  if (po<0.0 || pn<0.0) {
    msg_Debugging()<<METHOD<<"(): Kinematics does not fit."<<std::endl;
    return false;
  }
  po=sqrt(po);
  pn=sqrt(pn);
  double ecm(Q2-sij-mk2), gam(0.5*(ecm+pn));
  double zt(ecm/pn*(z-mk2/gam*(sij+mi2-mj2)/ecm));
  double ktt(sij*zt*(1.0-zt)-(1.0-zt)*mi2-zt*mj2);
  if (ktt<0.0 || gam<=0.0) {
    msg_Debugging()<<METHOD<<"(): Invalid kinematics."<<std::endl;
    return false;
  }
  ktt=sqrt(ktt);
  
  pk=pn/po*(pk-(Q2-mij2+mk2)/(2.0*Q2)*Q)+(Q2-sij+mk2)/(2.0*Q2)*Q;
  pj=Q-pk;
  pi=ktt*sin(phi)*l_perp;
  cms.BoostBack(pi);
  pi+=ktt*cos(phi)*n_perp+zt/pn*(gam*pj-sij*pk)+
    (mi2+ktt*ktt)/zt/pn*(pk-mk2/gam*pj);
  pj=Q-pk-pi;
  
  return true;
}

//////////////// begin Kinematics_FI ///////////////

double Kinematics_FI::GetKT2(const double &Q2,const double &y,const double &z,
			     const double mi2,const double mj2,const double ma2) const
{
  // Q2 = |(pi+pj-pk)|^2
  double pipj=(ma2-Q2-mi2-mj2)*y;
  return pipj*z*(1.0-z)-sqr(1.0-z)*mi2-sqr(z)*mj2; // transverse momentum ordering 
}

double Kinematics_FI::GetVirt(const double &Q2,const double &y,const double &z,
           const double mi2,const double mj2,const double ma2) const
{
  // this is the virtuality q2, for t=qbar2=q2-m2 subtract off mij2
  double pipj=(ma2-Q2-mi2-mj2)*y;
  return pipj + mi2 + mj2;
}

double Kinematics_FI::GetY(const double Q2,const double kt2,const double z,
			   const double mi2,const double mj2,const double ma2) const
{
  // to do
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+ma2) return -1.0;
  return 0;
}

double Kinematics_FI::GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double ma2) const
{
  // to do
  return 0;
}

bool Kinematics_FI::MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double ma2, const double mij2)
{
  // to do
  return 0;
}

/////////////// begin Kinematics_IF /////////////////

double Kinematics_IF::GetKT2(const double &Q2,const double &y,const double &z,
			     const double mi2,const double mj2,const double ma2) const
{
  return 0; // to do...
}

double Kinematics_IF::GetVirt(const double &Q2,const double &y,const double &z,
           const double mi2,const double mj2,const double ma2) const
{
  return 0; // to do...
}

double Kinematics_IF::GetY(const double Q2,const double kt2,const double z,
			   const double mi2,const double mj2,const double ma2) const
{
  // not used 
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+ma2) return -1.0;
  return 0;
}

double Kinematics_IF::GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double ma2) const
{
  // not used 
  return 0;
}

bool Kinematics_IF::MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double ma2, const double mij2)
{
  // not used
  return 0;
}

/////////////// begin Kinematics_II /////////////////

double Kinematics_II::GetKT2(const double &Q2,const double &y,const double &z,
			     const double mi2,const double mj2,const double ma2) const
{
  return 0; // to do...
}

double Kinematics_II::GetVirt(const double &Q2,const double &y,const double &z,
           const double mi2,const double mj2,const double ma2) const
{
  return 0; // to do...
}

double Kinematics_II::GetY(const double Q2,const double kt2,const double z,
			   const double mi2,const double mj2,const double ma2) const
{
  // not used 
  if (z<=0.0 || z>=1.0 || Q2<=mi2+mj2+ma2) return -1.0;
  return 0;
}

double Kinematics_II::GetYVirt(const double Q2, const double q2, const double mi2, const double mj2,
      const double ma2) const
{
  // not used 
  return 0;
}

bool Kinematics_II::MakeKinematics(const double z, const double y, const double phi, 
      Vec4D &pij,Vec4D &pk, Vec4D &pi, Vec4D &pj, const double mi2, const double mj2,
      const double ma2, const double mij2)
{
  // NOT CORRECT YET 
  // NOT USED YET 
  return 0;
}