#ifndef Four_Particle_Observables_H
#define Four_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class Four_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    std::vector<Flavour> m_flavs;
    void virtual Evaluate(Vec4D,Vec4D,Vec4D,Vec4D,double=1.) = 0; 
    void virtual Evaluate(double,double=1.); 
  public:
    Four_Particle_Observable_Base(std::vector<Flavour> &,int,double,double,int,
				  std::string=std::string(""));
    void virtual Evaluate(int,Vec4D *,Flavour *,double=1.);
    void virtual Evaluate(const Particle_List &,double=1.);
  };

  class Four_Particle_PlaneAngle : public ATOOLS::Four_Particle_Observable_Base {  
    void Evaluate(Vec4D,Vec4D,Vec4D,Vec4D,double=1.); 
  public:
    Four_Particle_PlaneAngle(std::vector<Flavour> &,int,double,double,int,
			     std::string=std::string("NRAngle:"));
  };
}

#endif
