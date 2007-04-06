#ifndef Analysis_Object_Base_h
#define Analysis_Object_Base_h

#include "Blob_List.H"
#include "Particle_List.H"

namespace ANALYSIS {
  class Analysis_Object_Base {
  private:
  public:
    Analysis_Object_Base() {}
    virtual ~Analysis_Object_Base()=0;
    virtual void Evaluate(const ATOOLS::Particle_List & plist,double weight, int ncount)=0;
    virtual void Evaluate(const ATOOLS::Blob_List & blobs, double weight, int ncount)=0;
  };
}
#endif
