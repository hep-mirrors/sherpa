#ifndef Two_Particle_Observables_H
#define Two_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Two_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour      m_flav1,m_flav2;
    void virtual Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1.,int =1) = 0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    Two_Particle_Observable_Base(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
				 int,double,double,int,const std::string & =std::string(""),const std::string & =std::string("FinalState"));
    void virtual Evaluate(int,const ATOOLS::Vec4D *,const ATOOLS::Flavour *,double=1., int=1);
    void virtual Evaluate(const ATOOLS::Particle_List &,double=1., int=1);
  };

  class Two_Particle_Mass : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & ,const ATOOLS::Vec4D &,double=1.,int=1); 
  public:
    Two_Particle_Mass(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		      int,double,double,int,const std::string & =std::string("Mass:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_PT(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		    int,double,double,int,const std::string & =std::string("P_T:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Scalar_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_Scalar_PT(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
			   int,double,double,int,const std::string & =std::string("SP_T:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

#ifdef ROOT_SUPPORT
  class Two_Particle_Angles : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_Angles(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
			int,double,double,int,const std::string & =std::string("SP_T:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };
#endif

  class Two_Particle_Eta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_Eta(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		    int,double,double,int,const std::string & =std::string("P_T:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

class Two_Particle_DEta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_DEta(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		      int,double,double,int,const std::string & =std::string("DEta:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DPhi : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_DPhi(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		      int,double,double,int,const std::string & =std::string("DPhi:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DR : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
  public:
    Two_Particle_DR(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
		    int,double,double,int,const std::string & =std::string("DR:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
