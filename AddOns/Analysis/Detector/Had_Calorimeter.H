#ifndef Analysis_Detector_Had_Calorimeter_H
#define Analysis_Detector_Had_Calorimeter_H

#include "AddOns/Analysis/Detector/Detector_Element.H"
#include "ATOOLS/Math/Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Had_Calorimeter : public Detector_Element {
  private:
  public:
    Had_Calorimeter(Primitive_Analysis *);
     virtual ~Had_Calorimeter();

    void Evaluate(const ATOOLS::Blob_List &,double,double) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);
  };
}

#endif
