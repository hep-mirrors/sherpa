#ifndef One_Particle_Observables_H
#define One_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class One_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    Flavour     m_flav;
    void virtual Evaluate(Vec4D,double=1.) = 0; 
    void virtual Evaluate(double,double=1.); 
  public:
    One_Particle_Observable_Base(Flavour &,int,double,double,int,
				 std::string=std::string(""));
    void virtual Evaluate(int,Vec4D *,Flavour *,double=1.);
    void virtual Evaluate(const Particle_List &,double=1.);
  };

  class One_Particle_ET : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(Vec4D,double=1.); 
  public:
    One_Particle_ET(Flavour &,int,double,double,int,std::string=std::string("E_T:"));
  };

  class One_Particle_PT : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(Vec4D,double=1.); 
  public:
    One_Particle_PT(Flavour &,int,double,double,int,std::string=std::string("p_T:"));
  };

  class One_Particle_Eta : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(Vec4D,double=1.); 
  public:
    One_Particle_Eta(Flavour &,int,double,double,int,std::string=std::string("eta:"));
  };

  class One_Particle_E : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(Vec4D,double=1.); 
  public:
    One_Particle_E(Flavour &,int,double,double,int,std::string=std::string("E:"));
  };
}

#endif
