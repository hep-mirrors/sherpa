#ifndef Two_Particle_Observables_H
#define Two_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class Two_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    Flavour      m_flav1,m_flav2;
    void virtual Evaluate(const Vec4D &,const Vec4D &,double=1.,int =1) = 0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    Two_Particle_Observable_Base(const Flavour &,const Flavour &,int,double,double,int,
				 const std::string & =std::string(""));
    void virtual Evaluate(int,const Vec4D *,const Flavour *,double=1., int=1);
    void virtual Evaluate(const Particle_List &,double=1., int=1);
  };

  class Two_Particle_Mass : public ATOOLS::Two_Particle_Observable_Base {  
    void Evaluate(const Vec4D & ,const Vec4D &,double=1.,int=1); 
  public:
    Two_Particle_Mass(const Flavour &,const Flavour &,
		      int,double,double,int,const std::string & =std::string("Mass:"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_PT : public ATOOLS::Two_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,const Vec4D &,double=1., int=1); 
  public:
    Two_Particle_PT(const Flavour &,const Flavour &,
		    int,double,double,int,const std::string & =std::string("P_T:"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Eta : public ATOOLS::Two_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,const Vec4D &,double=1., int=1); 
  public:
    Two_Particle_Eta(const Flavour &,const Flavour &,
		    int,double,double,int,const std::string & =std::string("P_T:"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
