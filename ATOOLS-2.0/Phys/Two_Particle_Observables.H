#ifndef Two_Particle_Observables_H
#define Two_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class Two_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    Flavour      m_flav1,m_flav2;
    void virtual Evaluate(Vec4D,Vec4D,double=1.) = 0; 
    void virtual Evaluate(double,double=1.); 
  public:
    Two_Particle_Observable_Base(Flavour &,Flavour &,int,double,double,int,
				 std::string=std::string(""));
    void virtual Evaluate(int,Vec4D *,Flavour *,double=1.);
    void virtual Evaluate(const Particle_List &,double=1.);
  };

  class Two_Particle_Mass : public ATOOLS::Two_Particle_Observable_Base {  
    void Evaluate(Vec4D,Vec4D,double=1.); 
  public:
    Two_Particle_Mass(Flavour &,Flavour &,
		      int,double,double,int,std::string=std::string("Mass:"));
  };

  class Two_Particle_PT : public ATOOLS::Two_Particle_Observable_Base {  
    void Evaluate(Vec4D,Vec4D,double=1.); 
  public:
    Two_Particle_PT(Flavour &,Flavour &,
		    int,double,double,int,std::string=std::string("P_T:"));
  };
}

#endif
