#ifndef Analysis_Detector_Detector_Segment_H
#define Analysis_Detector_Detector_Segment_H

#include "AddOns/Analysis/Detector/Cell.H"
#include "AddOns/Analysis/Main/Analysis_Object.H"
#include <set>

namespace ANALYSIS {
  class Detector_Segment : public Analysis_Object {
    double     m_etamin,m_etamax;
    long int   m_neta,m_nphi;

    Etastrip * p_first, * p_last;
  public:
    Detector_Segment(const double,const double,
		     const long int,const long int);
    ~Detector_Segment();

    void              Evaluate(const ATOOLS::Blob_List &,double,double) {}
    Analysis_Object * GetCopy() const;

    bool              InEtaRange(const double eta) const { return (eta>=m_etamin && eta<=m_etamax); }
    void              Dimensions(double &,double &) const;
    double            Etamin()                     const { return m_etamin; }
    double            Etamax()                     const { return m_etamax; }
    Etastrip        * GetFirst()                         { return p_first; }
    Etastrip        * GetLast()                          { return p_last; }

    void Reset();
    Cell * AddDeposit(const double,const double,const double);
    Cell * AddParticle(const double,const double,ATOOLS::Particle *,double=-1.);
    Cell * LocateCell(const double,const double);
  };

  class DS_Order {
  public:
    bool operator()(const Detector_Segment * a,const Detector_Segment * b) const {
      return (a->Etamin()<b->Etamin())?0:1;
    }
  };

  typedef std::set<Detector_Segment *,DS_Order> DS_List;
}

#endif
