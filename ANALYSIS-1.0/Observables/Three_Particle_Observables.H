#ifndef Three_Particle_Observables_H
#define Three_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Three_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour      m_flav1,m_flav2,m_flav3;
    void virtual Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
			  const ATOOLS::Vec4D & mom3, double weight, int ncount) = 0; 
  public:
    Three_Particle_Observable_Base(const String_Matrix & parameters, const std::string & obsname );
    Three_Particle_Observable_Base(const Three_Particle_Observable_Base*);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };
  
  class Three_Particle_Y : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_Y(const String_Matrix &);
    Three_Particle_Y(const Three_Particle_Y * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_DEta : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DEta(const String_Matrix &);
    Three_Particle_DEta(const Three_Particle_DEta * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_DPhi : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DPhi(const String_Matrix &);
    Three_Particle_DPhi(const Three_Particle_DPhi * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_DR : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DR(const String_Matrix &);
    Three_Particle_DR(const Three_Particle_DR * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_3Mass2 : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_3Mass2(const String_Matrix &);
    Three_Particle_3Mass2(const Three_Particle_3Mass2 * old);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_3Mass : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_3Mass(const String_Matrix &);
    Three_Particle_3Mass(const Three_Particle_3Mass * old);
    Primitive_Observable_Base * Copy() const;
  };

}
#endif
