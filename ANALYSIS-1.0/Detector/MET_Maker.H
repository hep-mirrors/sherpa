#ifndef MET_Maker_h
#define MET_Maker_h

#include "Object_Definition_Base.H"
#include "Particle.H"

namespace ANALYSIS {
  class Primitive_Analysis;

  class MET_Maker : public Object_Definition_Base {
  private:
    double m_cutHad, m_cutEM;
  public:
    MET_Maker(Primitive_Analysis *,const std::string="ET_UP");
    ~MET_Maker();

    void ReconstructObjects(ATOOLS::Particle_List *,ATOOLS::Vec4D & METvector);
    void SetCuts(const double cutHad,const double cutEM);
  };
}

#endif
