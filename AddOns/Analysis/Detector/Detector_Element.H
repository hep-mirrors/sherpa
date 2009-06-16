#ifndef Analysis_Detector_Detector_Element_H
#define Analysis_Detector_Detector_Element_H

#include "AddOns/Analysis/Detector/Detector_Segment.H"
#include "AddOns/Analysis/Main/Analysis_Object.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Detector_Element : public Analysis_Object {
  protected:
    std::string m_name;
    bool m_istested;

    DS_List             m_segments;
    std::list<Cell *>   m_hitcells;
    std::list<Track *>  m_tracks;
  public:
    Detector_Element(Primitive_Analysis *,const std::string name);
    Detector_Element(Detector_Segment *);
    virtual ~Detector_Element();

    virtual void AddDetectorSegment(Detector_Segment * seg);
    virtual DS_List            * GetDetectorSegments() { return (&m_segments); }
    virtual std::list<Cell *>  * GetHitCells() { return (&m_hitcells); }
    virtual Cell               * GetCell(const double eta,const double phi);
    virtual void GetTracks(std::list<Track *> &,const double eta=0,const double phi=0,
			   const double R2=0,kf_code=kf_none);

    virtual void Evaluate(const ATOOLS::Blob_List &,double,int);

    virtual Analysis_Object *GetCopy() const=0;
    virtual void Reset()=0;
    virtual bool Fill(const double E,const double eta,const double phi,ATOOLS::Particle * part)=0;


    virtual std::string Name();
    virtual bool        IsTested() const { return m_istested; }

    virtual void PrintHits() const;
  };
}

#endif
