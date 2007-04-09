#ifndef Primitive_Hadron_Calorimeter_H
#define Primitive_Hadron_Calorimeter_H

#include "Primitive_Detector.H"
#include "Particle_List.H"
#include "Vector.H"
#include "Message.H"

namespace ATOOLS { class Particle_Qualifier_Base; }

namespace ANALYSIS {
  class Primitive_Hadron_Calorimeter : public Primitive_Detector_Element {
  private:
    ATOOLS::Particle_Qualifier_Base *p_qualifier;

    double SmearEnergy(const ATOOLS::Particle *);
    void   AddNoise();
  public:
    Primitive_Hadron_Calorimeter(const int,const int,
				 const double=-5.,const double=5.,
				 const double=0.,const double=2.*M_PI);
    ~Primitive_Hadron_Calorimeter();
    void Fill(const ATOOLS::Particle_List *);
    void Extract(ATOOLS::Particle_List *);
    void Reset();
    void Print(std::ostream & =std::cout);

    Primitive_Detector_Element * Copy() const ;
  };
}

#endif
