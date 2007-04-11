#ifndef Primitive_Detector_Element_H
#define Primitive_Detector_Element_H

#include "Analysis_Object.H"
#include "Particle_Qualifier.H"
#include "Particle_List.H"
#include "Getter_Function.H"

namespace ANALYSIS {
  struct pde {
    enum code {
      cells=0,
      tracks=1
    };
  };

  class Primitive_Detector_Element : public Analysis_Object {
  protected:
    pde::code   m_mode;
    long int    m_neta, m_nphi;
    double      m_etamin,m_etamax,m_phimin,m_phimax;
    double      m_delta_eta, m_delta_phi;
    std::string m_name;
    std::vector<std::vector< double> >            m_cells;
    std::vector<std::vector< double> >::iterator  m_etastrip;
    std::map<std::pair<long int,long int>,double> m_tracks;

    ATOOLS::Particle_Qualifier_Base * p_qualifier;

    virtual bool FillParticleInDetectorElement(const ATOOLS::Particle *) =0;
    void         MatchCell(const double,const double,long int &,long int &) const;
    void         MatchCell(const ATOOLS::Vec4D,long int &,long int &) const;
    void         PseudoRapidityNAzimuthalAngle(const ATOOLS::Vec4D &,double &,double &); 
    virtual void AddNoise() {}
  public:
    Primitive_Detector_Element(const pde::code,const long int,const long int,
			       const double=-5.,const double=5.,
			       const double=0.,const double=2.*M_PI);
    virtual ~Primitive_Detector_Element();

    virtual void Evaluate(const ATOOLS::Blob_List &,double,int) {}
    virtual Analysis_Object * GetCopy() { return NULL; }
    virtual void Fill(const ATOOLS::Particle_List *)=0;
    virtual void Extract(ATOOLS::Particle_List *)   =0;
    virtual void Reset();
 
    void          GetNumbersOfCells(long int &,long int &);
    void          GetDimensions(long int &,long int &,double &,double &,double &, double &); 
    double        Cell(const long int,const long int) const;
    void          PseudoRapidityNAzimuthalAngle(const long int i,const long int j,
						double & eta, double & phi);
    ATOOLS::Vec4D ReconstructMasslessFourMom(const long int i,const long int j);
    void          SetName(std::string);
    std::string   Name() const;
  };


  typedef std::map<std::string,Primitive_Detector_Element *> String_DetectorElement_Map;
  typedef String_DetectorElement_Map::iterator               String_DetectorElement_Iter;
}

#endif
