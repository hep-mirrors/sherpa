#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace ATOOLS;

template<> const Vec4D Vec4D::XVEC=Vec4D(1.,1.,0.,0.);
template<> const Vec4D Vec4D::YVEC=Vec4D(1.,0.,1.,0.);
template<> const Vec4D Vec4D::ZVEC=Vec4D(1.,0.,0.,1.);

template<> const Vec3D Vec3D::XVEC=Vec3D(1.,0.,0.);
template<> const Vec3D Vec3D::YVEC=Vec3D(0.,1.,0.);
template<> const Vec3D Vec3D::ZVEC=Vec3D(0.,0.,1.);

template<> double Vec4D::CosPhi() const {
  return Max(Min(m_x[1]/PPerp(),1.0),-1.0);
}
template<> double Vec4D::SinPhi() const {
  if (PPerp() == 0.0) return 0.0;
  return Max(Min(m_x[2]/PPerp(),1.0),-1.0);
}
template<> double Vec4D::Phi() const {
  if(m_x[2]>0.) return acos(CosPhi());
  else return -acos(CosPhi());
}
template<> double Vec4D::CosTheta() const {
  return Max(Min(m_x[3]/PSpat(),1.0),-1.0);
}
template<> double Vec4D::SinTheta() const { 
  return Max(Min(sqrt(PPerp2()/PSpat2()),1.0),-1.0);
}
template<> double Vec4D::Theta() const {
  return acos(CosTheta());
}
template<> double Vec4D::Eta() const {
  double pt2=PPerp2();
  double pp =P();
  double pz =dabs(m_x[3]);
  double sn =Sign(m_x[3]);
  if (pt2<1.e-10*pp*pp) {
    return sn*20.;
  }
  return sn*0.5*log(sqr(pp+pz)/pt2);
}
template<> double Vec4D::CosTheta(const Vec4D& ref) const {
  Vec3D pref=Vec3D(ref), p=Vec3D(*this);
  return Max(Min(pref*p/(pref.Abs()*p.Abs()),1.0),-1.0);
}
template<> double Vec4D::Theta(const Vec4D& ref) const {
  return acos(CosTheta(ref));
}
template<> double Vec4D::Eta(const Vec4D& ref) const {
  double cos=CosTheta(ref);
  return 0.5*log(sqr(1.0+cos)/(1.0-cos*cos));
}
template<> double Vec4D::CosDPhi(const Vec4D& ref) const {
  Vec3D pref=Vec3D(ref[1],ref[2],0.0), p=Vec3D(m_x[1],m_x[2],0.0);
  return Max(Min(pref*p/(pref.Abs()*p.Abs()),1.0),-1.0);
}
template<> double Vec4D::DPhi(const Vec4D& ref) const {
  return acos(CosDPhi(ref));
}
template<> double Vec4D::DEta(const Vec4D& ref) const {
  return Eta()-ref.Eta();
}
template<> double Vec4D::DY(const Vec4D& ref) const {
  return Y()-ref.Y();
}
template<> double Vec4D::DR(const Vec4D& ref) const {
  return sqrt(sqr(DPhi(ref))+sqr(DEta(ref)));
}

std::istream& ATOOLS::operator>>(std::istream& s,Vec4D& vec)
{
  std::string out;
  s>>out;
  if (out.length()==0 || out[0]!='(' || out[out.length()-1]!=')')
    THROW(critical_error,"String to vector translation failed.");
  out=out.substr(0,out.length()-1).substr(1);
  for (short unsigned int i=0;i<4;++i) {
    size_t pos=out.find(",");
    vec[i]=ToType<double>(out.substr(0,pos));
    if (pos!=std::string::npos) out=out.substr(pos+1);
    else out="";
  }
  if (out.length()>0)
    THROW(critical_error,"Vector is not a four vector.");
  return s;
}

std::istream& ATOOLS::operator>>(std::istream& s,Vec3D& vec)
{
  std::string out;
  s>>out;
  if (out.length()==0 || out[0]!='(' || out[out.length()-1]!=')')
    THROW(critical_error,"String to vector translation failed.");
  out=out.substr(0,out.length()-1).substr(1);
  for (short unsigned int i=0;i<3;++i) {
    size_t pos=out.find(",");
    vec[i]=ToType<double>(out.substr(0,pos));
    if (pos!=std::string::npos) out=out.substr(pos+1);
    else out="";
  }
  if (out.length()>0)
    THROW(critical_error,"Vector is not a three vector.");
  return s;
}

bool ATOOLS::IsEqual(const Vec4D& v1, const Vec4D& v2, const double crit)
{
  double maxp=Max(dabs(v1[0]),Max(dabs(v1[1]),Max(dabs(v1[2]),dabs(v1[3])))); 
  double q(IsZero(maxp)?1.0:1.0/maxp);
  for(short int i=0;i<4;i++) 
   if (dabs(q*(v1[i]-v2[i]))>crit &&
       !(dabs(v1[i])<=crit && 
	 dabs(v2[i])<=crit)) return false;
  return true;
}

bool ATOOLS::IsEqual(const Vec3D& v1, const Vec3D& v2, const double crit)
{
  double maxp=Max(dabs(v1[1]),Max(dabs(v1[2]),dabs(v1[3]))); 
  double q=1.;
  if (!IsZero(maxp)) q=1./maxp;
  for(short int i=1;i<4;i++) {
    if (dabs(q*(v1[i]-v2[i]))>crit &&
	!(dabs(v1[i])<=crit && 
	  dabs(v2[i])<=crit)) return false;
  }
  return true;
}
