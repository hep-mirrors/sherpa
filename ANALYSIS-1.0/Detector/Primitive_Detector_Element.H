#ifndef Primitive_Detector_Element_H
#define Primitive_Detector_Element_H

#include "Particle_Qualifier.H"
#include "Particle_List.H"
#include "Getter_Function.H"

namespace ANALYSIS {
  class Primitive_Detector_Element {
  protected:
    int         m_neta, m_nphi;
    double      m_etamin,m_etamax,m_phimin,m_phimax;
    double      m_delta_eta, m_delta_phi;
    std::string m_name;
    double   ** p_cells;

    ATOOLS::Particle_Qualifier_Base * p_qualifier;

    virtual void FillParticleInDetectorElement(const ATOOLS::Particle *) =0;
    void         MatchCell(const double,const double,int &,int &) const;
    void         MatchCell(const ATOOLS::Vec4D,int &,int &) const;
    void         PseudoRapidityNAzimuthalAngle(const ATOOLS::Vec4D &,double &,double &); 
  public:
    Primitive_Detector_Element(const int,const int,
			       const double=-5.,const double=5.,
			       const double=0.,const double=2.*M_PI);
    virtual ~Primitive_Detector_Element();
    virtual Primitive_Detector_Element * Copy() const=0;
    virtual void Fill(const ATOOLS::Particle_List *) =0;
    virtual void Extract(ATOOLS::Particle_List *)    =0;
    virtual void Reset()                             =0;
 
    void          GetDimensions(int &,int &,double &,double &,double &, double &); 
    double        Cell(const int,const int) const;
    void          PseudoRapidityNAzimuthalAngle(const int i,const int j,
						double & eta, double & phi);
    ATOOLS::Vec4D ReconstructMasslessFourMom(const int i,const int j);
    void          SetName(std::string);
    std::string   Name() const;
  };


  typedef std::map<std::string,Primitive_Detector_Element *> String_DetectorElement_Map;
  typedef String_DetectorElement_Map::iterator               String_DetectorElement_Iter;

  typedef ATOOLS::Getter_Function<Primitive_Detector_Element,std::string> 
  Primitive_Detector_Element_Getter;
}

#endif
