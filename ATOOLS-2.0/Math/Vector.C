/*
 * vec.C: Definitions for
 * 3 dim. euklidean vectors and 4 dim. Minkowski vectors
 * Poincare transformations
 */


#include "Vector.H"


using namespace AMATOOLS;

// "reference" to tags for computational constructors
const Tag::Tsum Tag::sum={};
const Tag::Tdiff Tag::diff={};
const Tag::Tsmul Tag::smul={};
const Tag::Tcross Tag::cross={};

const Vec4D XVEC=Vec4D(1.,1.,0.,0.);
const Vec4D YVEC=Vec4D(1.,0.,1.,0.);
const Vec4D ZVEC=Vec4D(1.,0.,0.,1.);

std::ostream& AMATOOLS::operator<< (std::ostream& s, const Vec4D& vec)
{
  return s<<'('<<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
}

std::ostream& AMATOOLS::operator<< (std::ostream& s, const Vec3D& vec)
{
  return s<<'('<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
}



