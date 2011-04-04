#include "ATOOLS/Math/Poincare.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;

Poincare::Poincare(): 
  m_status(0), m_beta(1.,0.,0.,0.), m_rsq(1.), m_usen(false)
{
}
  
Poincare::Poincare(Vec4D v1, Vec4D v2,int mode):
  m_status(mode?3:2), m_beta(1.,0.,0.,0.), m_rsq(1.), m_usen(false)
{
  if (m_status==3) {
    m_beta=v1;
    m_n=v2;
    return;
  }
  ATOOLS::Vec3D b(v1), a(v2);
  m_ct=a*b/(a.Abs()*b.Abs());
  m_n=Vec4D(0.0,cross(a,b));
  m_usen=!IsZero(m_n.PSpat());
  if (!m_usen) {
    m_n=Vec4D(0.0,cross(a,Vec3D::XVEC));
    m_usen=!IsZero(m_n.PSpat());
    if (!m_usen) {
      m_n=Vec4D(0.0,cross(a,Vec3D::YVEC));
      m_usen=!IsZero(m_n.PSpat());
      if (!m_usen) {
	m_n=Vec4D(0.0,cross(a,Vec3D::ZVEC));
	m_usen=!IsZero(m_n.PSpat());
      }
    }
  }
  m_usen&=m_ct<1.0;
  m_nsq=m_n.PSpat2();
  m_nabs=sqrt(m_nsq);
  if(sqr(m_ct)>1.0) {
    if (!IsEqual(sqr(m_ct),1.0)) {
      int prc(msg_Error().precision());
      msg_Error().precision(14);
      msg_Error()<<METHOD<<"(): cos\\theta = "<<m_ct
		 <<". Set to +-1. "<<(m_usen?"Will":"Won't")
		 <<" rotate."<<std::endl;
      msg_Error().precision(prc);
    }
    m_ct=m_ct>0.0?1.0:-1.0;
  }
  m_st=-sqrt(1.0-sqr(m_ct));
  if (m_usen) {
    a=1.0/a.Abs()*a;
    b=1.0/b.Abs()*b;
    Vec3D n(m_n), at(n*(n*a)/m_nsq), ap(a-at);
    Vec3D ta(at+m_ct*ap-m_st/m_nabs*cross(n,ap));
    if (!ATOOLS::IsZero((ta-b).Sqr(),1.0e-6)) 
      msg_Error()<<METHOD<<"(): Inaccurate rotation {\n"
 		 <<"  a    = "<<a<<"\n"
		 <<"  b    = "<<b<<"\n"
		 <<"  a'   = "<<ta<<" -> rel. dev. ("
		 <<ta[1]/b[1]-1.0<<","<<ta[2]/b[2]-1.0<<","
		 <<ta[3]/b[3]-1.0<<")\n"
		 <<"  m_ct = "<<m_ct<<"\n"
		 <<"  m_st = "<<m_st<<"\n"
		 <<"  m_n  = "<<Vec3D(m_n)
		 <<"\n}"<<std::endl;
  }
}
  
Poincare::Poincare(Vec4D v): 
  m_status(1), m_beta(v), m_rsq(v.Mass()), m_usen(false)
{
  if (IsZero(m_rsq) || IsEqual(m_rsq,-m_beta[0])) m_status=0;
}
   
void Poincare::Boost(Vec4D& v) const
{
  if (m_status==0) return;
  double v0((m_beta[0]*v[0]-Vec3D(m_beta)*Vec3D(v))/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_beta[0]));
  v=Vec4D(v0,Vec3D(v)-c1*Vec3D(m_beta)); 
}
  
void Poincare::BoostBack(Vec4D& v) const
{
  if (m_status==0) return;
  double v0((m_beta[0]*v[0]+Vec3D(m_beta)*Vec3D(v))/m_rsq);
  double c1((v[0]+v0)/(m_rsq+m_beta[0]));
  v=Vec4D(v0,Vec3D(v)+c1*Vec3D(m_beta));  
}
  
  
void Poincare::Rotate(Vec4D& v) const
{
  if (m_usen) {
    Vec3D n(m_n), a(v), at(n*(n*a)/m_nsq), ap(a-at);
    v=Vec4D(v[0],at+m_ct*ap+m_st/m_nabs*cross(n,ap));
  }
}
  
void Poincare::RotateBack(Vec4D& v) const
{
  if (m_usen) {
    Vec3D n(-1.0*m_n), a(v), at(n*(n*a)/m_nsq), ap(a-at);
    v=Vec4D(v[0],at+m_ct*ap+m_st/m_nabs*cross(n,ap));
  }
}
  
void Poincare::Lambda(Vec4D& v) const
{
  if (m_status!=3) THROW(fatal_error,"Invalid function call");
  v=v-2.0*(v*(m_beta+m_n))/(m_beta+m_n).Abs2()*(m_beta+m_n)
    +2.0*(v*m_beta)/m_beta.Abs2()*m_n;
}

void Poincare::LambdaBack(Vec4D& v) const
{
  if (m_status!=3) THROW(fatal_error,"Invalid function call");
  v=v-2.0*(v*(m_beta+m_n))/(m_beta+m_n).Abs2()*(m_beta+m_n)
    +2.0*(v*m_n)/m_n.Abs2()*m_beta;
}

bool Poincare::CheckBoost() 
{
  if (m_beta.Abs2()<=0.) return false;
  return true;
}

bool Poincare::CheckRotation() 
{
  if (m_usen) {
    if ((!(m_ct>=0) && !(m_ct<0)) || (!(m_st>=0) && !(m_st<0)) ||
	m_n.Nan()) return false;
  }
  return true;
}

bool Poincare::CheckLambda() 
{
  if (!IsEqual(m_beta.Abs2(),m_n.Abs2())) return false;
  return true;
}

void Poincare::Invert() 
{
  if (m_status==3) {
    std::swap<Vec4D>(m_beta,m_n);
    return;
  }
  for (short unsigned int i(1);i<4;++i) {
    m_beta[i]=-m_beta[i];
    m_n[i]=-m_n[i];
  }
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
