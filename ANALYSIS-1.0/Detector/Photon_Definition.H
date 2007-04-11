#ifndef Photon_Definition_h
#define Photon_Definition_h

#include "Object_Definition_Base.H"

namespace ANALYSIS {
  class Photon_Definition : public Object_Definition_Base {
  public:
    Photon_Definition(const std::string="ET_UP");

    void FillMCTruthList(ATOOLS::Particle_List *,Primitive_Detector *);
    void FillSimpleDetectorList(ATOOLS::Particle_List *,Primitive_Detector *);
  };
}
#endif
