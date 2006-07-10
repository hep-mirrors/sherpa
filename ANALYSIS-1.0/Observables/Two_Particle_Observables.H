#ifndef Two_Particle_Observables_H
#define Two_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Two_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour      m_flav1,m_flav2;
    void virtual Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
			  double weight, int ncount) = 0;
  public:
    Two_Particle_Observable_Base(const String_Matrix & parameters, const std::string & obsname);
    Two_Particle_Observable_Base(const Two_Particle_Observable_Base*);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };

  class Two_Particle_Mass : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_Mass(const String_Matrix &);
    Two_Particle_Mass(const Two_Particle_Mass * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_PT(const String_Matrix &);
    Two_Particle_PT(const Two_Particle_PT * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Scalar_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_Scalar_PT(const String_Matrix &);
    Two_Particle_Scalar_PT(const Two_Particle_Scalar_PT * old);
    Primitive_Observable_Base * Copy() const;
  };

// #ifdef USING__ROOT
//   class Two_Particle_Angles : public Two_Particle_Observable_Base {  
//     void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., int=1); 
//   public:
//     Two_Particle_Angles(const String_Matrix &,const Primitive_Observable_Defaults *);
//     Two_Particle_Angles(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
// 			int,double,double,int,const std::string & =std::string("SP_T:"),const std::string & =std::string("FinalState"));
//     Two_Particle_Angles(const Two_Particle_Angles * old);
//     Primitive_Observable_Base * Copy() const;
//   };
// #endif

  class Two_Particle_Eta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_Eta(const String_Matrix &);
    Two_Particle_Eta(const Two_Particle_Eta * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DEta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_DEta(const String_Matrix &);
    Two_Particle_DEta(const Two_Particle_DEta * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Y : public Two_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2,
		  double weight, int ncount);
  public:
    Two_Particle_Y(const String_Matrix &);
    Two_Particle_Y(const Two_Particle_Y * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DY : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_DY(const String_Matrix &);
    Two_Particle_DY(const Two_Particle_DY * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DPhi : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_DPhi(const String_Matrix &);
    Two_Particle_DPhi(const Two_Particle_DPhi * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Angle : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_Angle(const String_Matrix &);
    Two_Particle_Angle(const Two_Particle_Angle * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DR : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_DR(const String_Matrix &);
    Two_Particle_DR(const Two_Particle_DR * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_CMS_Angle : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_CMS_Angle(const String_Matrix &);
    Two_Particle_CMS_Angle(const Two_Particle_CMS_Angle * old);
    Primitive_Observable_Base * Copy() const;
  };
   
  class Two_Particle_Mass2 : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, int ncount); 
  public:
    Two_Particle_Mass2(const String_Matrix &);
    Two_Particle_Mass2(const Two_Particle_Mass2 * old);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
