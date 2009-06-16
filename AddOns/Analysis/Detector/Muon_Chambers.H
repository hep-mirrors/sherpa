#ifndef Analysis_Detector_Muon_Chambers_H
#define Analysis_Detector_Muon_Chambers_H

#include "AddOns/Analysis/Detector/Detector_Element.H"
#include "ATOOLS/Math/Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Muon_Chambers : public Detector_Element {
  private:
    double m_etamin, m_etamax;
  public:
    Muon_Chambers(Primitive_Analysis *);
    ~Muon_Chambers();

    void Evaluate(const ATOOLS::Blob_List &,double,double) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);
    void SetAcceptance(const double etamin,const double etamax);
  };
}

#endif
