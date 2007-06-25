#ifndef Jet_Maker_h
#define Jet_Maker_h

#include "Object_Definition_Base.H"
#include "Simple_Cone.H"
#include "Particle.H"
#include "BTagging.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Jet_Maker : public Object_Definition_Base {
    BTagging    * p_btagging;
    Simple_Cone * p_simplecone;
    double        m_inv, m_spread;
    int           m_jetmode, m_correction;

    void IdentifyBs();
    void CorrectEnergies();
  public:
    Jet_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~Jet_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);
    void SetSimpleCone(const double Ecut,const double Estart,const double R);
    void SetECorrection(const double,const int=1);
    void SetBTagging(const double,const double,const double);
  };
}

#endif
