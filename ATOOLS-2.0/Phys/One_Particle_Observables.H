#ifndef One_Particle_Observables_H
#define One_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class One_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    Flavour     m_flav;
    void virtual Evaluate(const Vec4D &,double=1.,int=1) = 0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    One_Particle_Observable_Base(const Flavour &,int,double,double,int,
				 const std::string & =std::string(""));
    void virtual Evaluate(int,const Vec4D *,const Flavour *,double=1.,int ncount=1);
    void virtual Evaluate(const Particle_List &,double=1.,int ncount=1);
  };

  class One_Particle_ET : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,double=1.,int=1); 
  public:
    One_Particle_ET(const Flavour &,int,double,double,int,const std::string & =std::string("E_T:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_PT : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,double=1., int=1); 
  public:
    One_Particle_PT(const Flavour &,int,double,double,int,const std::string & =std::string("p_T:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Eta : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,double=1.,int=1); 
  public:
    One_Particle_Eta(const Flavour &,int,double,double,int,const std::string & =std::string("eta:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_E : public ATOOLS::One_Particle_Observable_Base {  
    void Evaluate(const Vec4D &,double=1., int=1); 
  public:
    One_Particle_E(const Flavour &,int,double,double,int,const std::string & =std::string("E:"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
