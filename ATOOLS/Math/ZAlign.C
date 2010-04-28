#include "ATOOLS/Math/ZAlign.H"

using namespace ATOOLS;

ZAlign::ZAlign(const Vec4D &pa,const Vec4D &pb,
	       const double &ma2,const double &mb2):
  m_pao(pa), m_pb(pb)
{
  Vec4D Q(pa+pb);
  double Q2=Q.Abs2(), papb=0.5*(Q2-ma2-mb2), ea=0.0;
  if (IsZero(mb2)) {
    ea=0.5*(papb+ma2*sqr(pb[3])/papb)/pb[0];
    m_pan=Vec4D(ea,0.0,0.0,sqrt(ea*ea-ma2));
    Vec4D ph(m_pan[0],0.0,0.0,-m_pan[3]);
    if (dabs((ph+pb).Abs2()-Q2)<
	dabs((m_pan+pb).Abs2()-Q2)) m_pan[3]=-m_pan[3];
  }
  else {
    double kr((pb.PMinus()*Q.PPlus())/(pb.PPlus()*Q.PMinus()));
    double eps(pb[0]*papb), kap(pb[3]*sqrt(eps*eps-ma2*mb2));
    if (dabs((eps+kap)/(eps-kap)-kr)<
	dabs((eps-kap)/(eps+kap)-kr)) kap=-kap;
    m_pan=Vec4D(eps+kap,0.0,0.0,eps-kap);
  }
  Vec4D pao(m_pao), pan(m_pan);
  m_cmso=Poincare(pao+pb);
  m_cmsn=Poincare(pan+pb);
  m_cmso.Boost(pao);
  m_cmsn.Boost(pan);
  m_zrot=Poincare(pao,pan);
}

void ZAlign::Align(Vec4D &p) const
{
  m_cmso.Boost(p);
  m_zrot.Rotate(p);
  m_cmsn.BoostBack(p);
}
