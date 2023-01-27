#include "ATOOLS/Math/Poincare.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

Poincare::Poincare(const Vec4D &v,const double &rsq):
  m_type(1), m_l(v), m_rsq(rsq>0.?rsq:v.Mass()) {}

Poincare::Poincare(const Vec4D &v1,const Vec4D &v2,int mode):
  m_type(mode?3:2), m_l(1.,0.,0.,0.), m_rsq(1.)
{
  if (m_type==3) {
    m_l=v1;
    m_t=v2;
    return;
  }
  Vec4D b(0.,Vec3D(v2)/v2.PSpat());
  m_l=Vec4D(0.,Vec3D(v1)/v1.PSpat());
  m_t=b+m_l*(m_l*b);
  double mt(m_t.PSpat2());
  if (mt!=0.) m_t*=1./sqrt(mt);
  int l[3]{1,2,3};
  double ml[4]={0.,dabs(m_l[1]),dabs(m_l[2]),dabs(m_l[3])};
  if (ml[l[2]]>ml[l[1]]) std::swap<int>(l[1],l[2]);
  if (ml[l[1]]>ml[l[0]]) std::swap<int>(l[0],l[1]);
  if (ml[l[2]]>ml[l[1]]) std::swap<int>(l[1],l[2]);
  double tdp(m_t[l[1]]*m_l[l[1]]+m_t[l[2]]*m_l[l[2]]);
  if (tdp!=0.) m_t[l[0]]=-tdp/m_l[l[0]];
  if (m_t.PSpat2()==0.) m_t[l[1]]=1.;
  m_omct=m_l.SmallOMCT(b);
  m_st=-m_t*b;
}

void Poincare::Boost(Vec4D &v) const
{
  double lv(m_l[1]*v[1]+m_l[2]*v[2]+m_l[3]*v[3]);
  double v0((m_l[0]*v[0]-lv)/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_l[0]));
  v=Vec4D(v0,Vec3D(v)-c1*Vec3D(m_l));
}

void Poincare::BoostBack(Vec4D &v) const
{
  double lv(m_l[1]*v[1]+m_l[2]*v[2]+m_l[3]*v[3]);
  double v0((m_l[0]*v[0]+lv)/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_l[0]));
  v=Vec4D(v0,Vec3D(v)+c1*Vec3D(m_l));
}

void Poincare::BoostBack(Vec4C &v) const
{
  Complex lv(m_l[1]*v[1]+m_l[2]*v[2]+m_l[3]*v[3]);
  Complex v0((m_l[0]*v[0]+lv)/m_rsq);
  Complex c1((v[0]+v0)/(m_rsq+m_l[0]));
  v=Vec4<Complex>(v0,Vec3<Complex>(v)+c1*Vec3<Complex>(m_l));
}

void Poincare::Rotate(Vec4D &v) const
{
  double vx(-m_l*v), vy(-m_t*v);
  v-=(m_omct*vx+m_st*vy)*m_l;
  v-=(-m_st*vx+m_omct*vy)*m_t;
}

void Poincare::RotateBack(Vec4D &v) const
{
  double vx(-m_l*v), vy(-m_t*v);
  v-=(m_omct*vx-m_st*vy)*m_l;
  v-=(m_st*vx+m_omct*vy)*m_t;
}

void Poincare::Lambda(Vec4D &v) const
{
  double m2=v.Abs2();
  v=v-2.0*(v*(m_l+m_t))/(m_l+m_t).Abs2()*(m_l+m_t)
    +2.0*(v*m_l)/m_l.Abs2()*m_t;
  v[0]=Sign(v[0])*sqrt(v.PSpat2()+m2);
}

void Poincare::LambdaBack(Vec4D &v) const
{
  double m2=v.Abs2();
  v=v-2.0*(v*(m_l+m_t))/(m_l+m_t).Abs2()*(m_l+m_t)
    +2.0*(v*m_t)/m_t.Abs2()*m_l;
  v[0]=Sign(v[0])*sqrt(v.PSpat2()+m2);
}

void Poincare::Invert() 
{
  if (m_type==3) { std::swap<Vec4D>(m_l,m_t); return; }
  if (m_type==2) { m_st=-m_st; return; }
  for (int i(1);i<4;++i) m_l[i]=-m_l[i];
}

Vec4D Poincare_Sequence::operator*(const Vec4D &p) const
{
  Vec4D np(p);
  for (const_iterator pit(begin());pit!=end();++pit) np=(*pit)*np;
  return np;
}

void Poincare_Sequence::Invert()
{
  Poincare_Sequence copy(*this);
  reverse_iterator cit(copy.rbegin());
  for (iterator pit(begin());pit!=end();++pit,++cit) {
    cit->Invert();
    *pit=*cit;
  }
}
