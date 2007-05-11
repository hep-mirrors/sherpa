#ifndef Tracker_h
#define Tracker_h

#include "Detector_Element.H"
#include "Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Tracker : public Detector_Element {
  private:
    double m_etamin, m_etamax;
  public:
    Tracker(Primitive_Analysis *);
    ~Tracker();

    void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);
    void SetAcceptance(const double etamin,const double etamax);

    void Print(ATOOLS::kf::code=ATOOLS::kf::none);
  };
}

#endif
