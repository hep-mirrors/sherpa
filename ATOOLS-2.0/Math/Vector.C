/*
 * vec.C: Definitions for
 * 3 dim. euklidean vectors and 4 dim. Minkowski vectors
 * Poincare transformations
 */


#include "Vector.H"
#include "Run_Parameter.H"



using namespace ATOOLS;

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

bool ATOOLS::operator==(const Vec4D& v1, const Vec4D& v2) 
{
  double maxp=Max(v1[0],Max(v1[1],Max(v1[2],v1[3]))); 
  double q=1.;
  if (!IsZero(maxp)) q=1./maxp;
  for(short int i=0;i<4;i++) {
    if (dabs(q*(v1[i]-v2[i]))>ATOOLS::rpa.gen.Accu()) return false;
  }
  return true;
}
