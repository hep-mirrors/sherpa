#ifndef Jet_Maker_h
#define Jet_Maker_h

#include "Object_Definition_Base.H"
#include "Simple_Cone.H"
#include "Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class Jet_Maker : public Object_Definition_Base {
    Simple_Cone * p_simplecone;
    double        m_inv;
    int           m_jetmode;

    void ReconstructSimpleCones();
    void CorrectEnergies();
  public:
    Jet_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~Jet_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);
    void SetSimpleCone(const double Ecut,const double Estart,const double R);
    void SetECorrection(const double);
  };
}

#endif
