#ifndef Had_Calorimeter_h
#define Had_Calorimeter_h

#include "Detector_Element.H"
#include "Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Had_Calorimeter : public Detector_Element {
  private:
  public:
    Had_Calorimeter(Primitive_Analysis *);
     virtual ~Had_Calorimeter();

    void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);
  };
}

#endif
