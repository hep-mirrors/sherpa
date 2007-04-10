#ifndef Primitive_Detector_H
#define Primitive_Detector_H

#include "Primitive_Detector_Element.H"
#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Primitive_Detector: public Primitive_Observable_Base {
    std::string                m_inlistname;
    String_DetectorElement_Map m_elements;
  public:
    Primitive_Detector(const std::string &inlist,const std::string &outlist);
    ~Primitive_Detector();

    Primitive_Observable_Base* Copy() const;
    
    void Evaluate(const ATOOLS::Blob_List &,double=1.,int=1);
    
    void Reset();

    void Add(Primitive_Detector_Element *);
    void Fill(const ATOOLS::Blob_List *);
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    Primitive_Detector_Element * GetElement(std::string name);

    std::string Name() const;
    void        Print();

    void AddSelector(const double etmin,const double etamin,
		     const double etamax,const double rmin,
		     const int bjets);

  };

  inline std::string Primitive_Detector::Name() const { return m_name; }
  inline void Primitive_Detector::Reset() {
    for (String_DetectorElement_Iter sdeiter=m_elements.begin();
	 sdeiter!=m_elements.end();sdeiter++) {
      sdeiter->second->Reset();
    }
  }
}

#endif
