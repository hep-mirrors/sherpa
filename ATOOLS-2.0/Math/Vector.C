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

const vec4d xvec=vec4d(1.,1.,0.,0.);
const vec4d yvec=vec4d(1.,0.,1.,0.);
const vec4d zvec=vec4d(1.,0.,0.,1.);

std::ostream& AMATOOLS::operator<< (std::ostream& s, const vec4d& vec)
{
  return s<<'('<<vec[0]<<','<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
}

std::ostream& AMATOOLS::operator<< (std::ostream& s, const vec3d& vec)
{
  return s<<'('<<vec[1]<<','<<vec[2]<<','<<vec[3]<<')';
}



