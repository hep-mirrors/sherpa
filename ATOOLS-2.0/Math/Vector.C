/*
 * vec.C: Definitions for
 * 3 dim. euklidean vectors and 4 dim. Minkowski vectors
 * Poincare transformations
 */


#include "Vector.H"
#include "Run_Parameter.H"
#include "Exception.H"
#include "MyStrStream.H"

using namespace ATOOLS;

double Vec4D::s_accu=1.0e-12;

// "reference" to tags for computational constructors
const Tag::Tsum Tag::sum={};
const Tag::Tdiff Tag::diff={};
const Tag::Tsmul Tag::smul={};
const Tag::Tcross Tag::cross={};

const Vec4D Vec4D::XVEC=Vec4D(1.,1.,0.,0.);
const Vec4D Vec4D::YVEC=Vec4D(1.,0.,1.,0.);
const Vec4D Vec4D::ZVEC=Vec4D(1.,0.,0.,1.);

const Vec3D Vec3D::XVEC=Vec3D(1.,0.,0.);
const Vec3D Vec3D::YVEC=Vec3D(0.,1.,0.);
const Vec3D Vec3D::ZVEC=Vec3D(0.,0.,1.);

Vec3D::Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tcross) 
{
  m_x[0]=v1[2]*v2[3]-v1[3]*v2[2];
  m_x[1]=v1[3]*v2[1]-v1[1]*v2[3];
  m_x[2]=v1[1]*v2[2]-v1[2]*v2[1];
}

std::ostream& ATOOLS::operator<< (std::ostream& s, const Vec4D& vec)
{
  return s<<'('<<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
}

std::ostream& ATOOLS::operator<< (std::ostream& s, const Vec3D& vec)
{
  return s<<'('<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
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

bool ATOOLS::operator==(const Vec4D& v1, const Vec4D& v2) 
{
  double maxp=Max(dabs(v1[0]),Max(dabs(v1[1]),Max(dabs(v1[2]),dabs(v1[3])))); 
  double q(IsZero(maxp)?1.0:1.0/maxp);
  for(short int i=0;i<4;i++) 
    if (dabs(q*(v1[i]-v2[i]))>Vec4D::Accu()) return false;
  return true;
}

bool Vec4D::Nan() const
{
  for(short unsigned int i(0);i<4;++i) 
    if (!(m_x[i]>=0.0) && !(m_x[i]<=0.0)) return true;
  return false;
}

bool Vec4D::IsZero() const
{
  for(short unsigned int i(0);i<4;++i) 
    if (!ATOOLS::IsZero(m_x[i])) return false;
  return true;
}
void Vec4D::ResetAccu()                
{ 
  s_accu=rpa.gen.Accu(); 
}

const double Vec4D::PPerp(const Vec4D &ref) const 
{ 
  Vec3D perp=1./ATOOLS::Max(Vec3D(ref).Abs(),1.e-12)*Vec3D(ref);
  perp=Vec3D(*this)-perp*(perp*Vec3D(*this));
  return perp.Abs();
}

const double Vec4D::PPerp2(const Vec4D &ref) const 
{ 
  return sqr(PPerp(ref)); 
}
  
const double Vec4D::Theta(const Vec4D &ref) const 
{ 
  return acos(CosTheta(ref));
}
  
const double Vec4D::Theta() const 
{
  return acos(CosTheta());
}

const double Vec4D::CosDPhi(const Vec4D &ref) const 
{ 
  Vec3D pref=Vec3D(ref[1],ref[2],0.0), p=Vec3D(m_x[1],m_x[2],0.0);
  return Max(Min(pref*p/(pref.Abs()*p.Abs()),1.0),-1.0);
}

const double Vec4D::CosTheta(const Vec4D &ref) const 
{ 
  Vec3D pref=Vec3D(ref), p=Vec3D(*this);
  return Max(Min(pref*p/(pref.Abs()*p.Abs()),1.0),-1.0);
}

const double Vec4D::Eta(const Vec4D &ref) const 
{ 
  double cos=CosTheta(ref);
  return (2.0*cos>M_PI?-0.5:0.5)*log(sqr(1.0+cos)/(1.0-cos*cos));
}
  
bool ATOOLS::operator==(const Vec3D& v1, const Vec3D& v2) 
{
  double maxp=Max(v1[1],Max(v1[2],v1[3])); 
  double q=1.;
  if (!IsZero(maxp)) q=1./maxp;
  for(short int i=1;i<4;i++) {
    if (!IsZero(q*(v1[i]-v2[i]))) return false;
  }
  return true;
}

