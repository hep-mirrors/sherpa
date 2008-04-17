#ifndef Photon_Definition_h
#define Photon_Definition_h

#include "Object_Definition_Base.H"

namespace ANALYSIS {
  class Photon_Definition : public Object_Definition_Base {
  public:
    Photon_Definition();
    void Evaluate(const ATOOLS::Particle_List & plist,double weight, int ncount);
    void Evaluate(const ATOOLS::Blob_List & blobs, double weight, int ncount);
  };
}
#endif
