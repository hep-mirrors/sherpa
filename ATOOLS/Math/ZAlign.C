#include "ATOOLS/Math/ZAlign.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

ZAlign::ZAlign(const Vec4D &pa,const Vec4D &pb,
	       const double &ma2,const double &mb2):
  m_pao(pa), m_pb(pb), m_stat(1)
{
  Vec4D Q(pa+pb);
  double Q2=Q.Abs2(), papb=0.5*(Q2-ma2-mb2);
  if (!IsEqual(papb,pa*pb,1.0e-3))
    msg_Error()<<METHOD<<"(): p_a*p_b = "<<papb
	       <<" vs. "<<pa*pb<<", rel. diff. "
	       <<papb/(pa*pb)-1.0<<std::endl;
  if (!IsZero(pb[1],1.0e-3) || !IsZero(pb[2],1.0e-3))
    msg_Error()<<METHOD<<"(): p_b not aligned -> "
	       <<pb<<std::endl;
  if (IsZero(mb2)) {
    double ea(0.5*(papb+ma2*sqr(pb[3])/papb)/pb[0]);
    if (ea*ea<ma2) m_stat=-1;
    m_pan=Vec4D(ea,0.0,0.0,sqrt(ea*ea-ma2));
    Vec4D ph(m_pan[0],0.0,0.0,-m_pan[3]);
    if (dabs((ph+pb).Abs2()-Q2)<
	dabs((m_pan+pb).Abs2()-Q2)) m_pan[3]=-m_pan[3];
  }
  else {
    double kap(papb*papb-ma2*mb2);
    if (kap<0.0) m_stat=-1;
    else kap=sqrt(kap);
    double kr((pb.PMinus()*Q.PPlus())/(pb.PPlus()*Q.PMinus()));
    double yo(log(kr)), yn(log((papb+mb2+kap)/(papb+mb2-kap)));
    if (dabs(yn+yo)<dabs(yn-yo)) kap=-kap;
    m_pan=Vec4D((pb[0]*papb+pb[3]*kap)/mb2,0.0,
		0.0,(pb[3]*papb+pb[0]*kap)/mb2);
  }
  if (!IsEqual(Q2,(m_pan+m_pb).Abs2(),1.0e-3))
    msg_Error()<<METHOD<<"(): Q = "<<sqrt(Q2)<<" vs. "
	       <<(m_pan+m_pb).Mass()<<", rel. diff. "
	       <<sqrt(Q2/(m_pan+m_pb).Abs2())-1.0<<std::endl;
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

Vec4D ZAlign::Align(const Vec4D &p) const
{
  Vec4D cp(p);
  m_cmso.Boost(cp);
  m_zrot.Rotate(cp);
  m_cmsn.BoostBack(cp);
  return cp;
}
