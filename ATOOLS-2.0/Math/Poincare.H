#ifndef poincare_h
#define poincare_h
#include "Vector.H"

namespace AMATOOLS {
  class Poincare {
    double mat[4][4];               
    vec4d beta;
    int status;  // 0 no operation, 1 boost, 2 rot, 3 boost+rot
  public:
    Poincare();                     // standard constructor Unity
    Poincare(vec4d v1, vec4d v2);   // rotations constructor
    Poincare(vec4d v);              // boost constuctor 
 
    void Boost(vec4d& v);          // boosts vectors in CMS (of vector given at initialization)
    void BoostBack(vec4d& v);      // boost back to LAB Frame

    void Rotate(vec4d& v);         // rotate
    void RotateBack(vec4d& v);     // rotate back

    vec4d operator*(const vec4d& vin) {
      vec4d v(vin);
      if ((status)&(1)) Boost(v);
      if ((status)&(2)) Rotate(v);
      return v;
    }
  };
}

#endif
