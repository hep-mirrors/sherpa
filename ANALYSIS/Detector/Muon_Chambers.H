#ifndef Muon_Chambers_h
#define Muon_Chambers_h

#include "Detector_Element.H"
#include "Histogram.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Muon_Chambers : public Detector_Element {
  private:
    double m_etamin, m_etamax;
  public:
    Muon_Chambers(Primitive_Analysis *);
    ~Muon_Chambers();

    void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    Analysis_Object * GetCopy() const;

    void Reset();
    bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part);
    void SetAcceptance(const double etamin,const double etamax);
  };
}

#endif
