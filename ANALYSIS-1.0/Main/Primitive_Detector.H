#ifndef Primitive_Detector_H
#define Primitive_Detector_H

#include "Blob_List.H"
#include "Particle_List.H"
#include "Vector.H"
#include "Message.H"

namespace ANALYSIS {
  class Primitive_Detector_Element {
  protected:
    int         m_nx, m_ny;
    std::string m_name;
    double   ** p_cells;
  public:
    Primitive_Detector_Element();
    Primitive_Detector_Element(std::string name);
    Primitive_Detector_Element(int neta,int nphi,std::string name);
    virtual ~Primitive_Detector_Element();
    virtual void Fill(const ATOOLS::Particle_List *) = NULL;
    virtual void Reset()                             = NULL;
    double       Cell(const int i,const int j);
    std::string  Name() const;
  };

  inline double Primitive_Detector_Element::Cell(const int i,const int j) {
    if (i>-1&&j>-1&&i<m_nx&&j<m_ny) return p_cells[i][j];
    else ATOOLS::msg.Error()<<"Error in Primitive_Detector_Element "<<m_name<<std::endl
			    <<"   GetCell("<<i<<","<<j<<") out of bounds, return 0."<<std::endl;
    return 0.;
  }
  inline std::string Primitive_Detector_Element::Name() const { return m_name; }



  typedef std::map<std::string,Primitive_Detector_Element *> String_DetectorElement_Map;
  typedef String_DetectorElement_Map::iterator               String_DetectorElement_Iter;


  class Primitive_Detector {
    std::string                m_name;
    String_DetectorElement_Map m_elements;
  public:
    Primitive_Detector();
    ~Primitive_Detector();

    void Add(Primitive_Detector_Element *);
    void Fill(const ATOOLS::Blob_List *);
    void Fill(const ATOOLS::Particle_List *);
    Primitive_Detector_Element * GetElement(std::string name);

    std::string Name() const;
    void        Print();
    void        Reset();
  };

  inline std::string Primitive_Detector::Name() const { return m_name; }
  inline void Primitive_Detector::Reset() {
    for (String_DetectorElement_Iter sdeiter=m_elements.begin();
	 sdeiter!=m_elements.end();sdeiter++) {
      sdeiter->second->Reset();
    }
  };
}

#endif
