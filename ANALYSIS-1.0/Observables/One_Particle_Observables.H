#ifndef One_Particle_Observables_H
#define One_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class One_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour m_flav;
    void virtual Evaluate(const ATOOLS::Vec4D &,double=1.,int=1) = 0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    One_Particle_Observable_Base(const ATOOLS::Flavour &,int,double,double,int,
				 const std::string & =std::string(""));
    void virtual Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double=1.,int ncount=1);
    void virtual Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  };

  class One_Particle_ET : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,double=1.,int=1); 
  public:
    One_Particle_ET(const ATOOLS::Flavour &,int,double,double,int,const std::string & =std::string("E_T:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_PT : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    One_Particle_PT(const ATOOLS::Flavour &,int,double,double,int,const std::string & =std::string("p_T:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Eta : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,double=1.,int=1); 
  public:
    One_Particle_Eta(const ATOOLS::Flavour &,int,double,double,int,const std::string & =std::string("eta:"));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_E : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    One_Particle_E(const ATOOLS::Flavour &,int,double,double,int,const std::string & =std::string("E:"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
