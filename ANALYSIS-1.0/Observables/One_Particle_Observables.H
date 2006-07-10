#ifndef One_Particle_Observables_H
#define One_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class One_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour m_flav;
    void virtual Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount) = 0; 
    void virtual Evaluate(double value, double weight, int ncount); 
  public:
    One_Particle_Observable_Base(const String_Matrix & parameters, const std::string & obsname);
    One_Particle_Observable_Base(const One_Particle_Observable_Base*);
    void virtual Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  const ATOOLS::Flavour * flavs,
			  double weight, int ncount);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };

  class One_Particle_ET : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_ET(const String_Matrix &);
    One_Particle_ET(const One_Particle_ET * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_PT : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_PT(const String_Matrix &);
    One_Particle_PT(const One_Particle_PT * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Eta : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_Eta(const String_Matrix &);
    One_Particle_Eta(const One_Particle_Eta * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_E : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_E(const String_Matrix &);
    One_Particle_E(const One_Particle_E * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_BeamAngle : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_BeamAngle(const String_Matrix &);
    One_Particle_BeamAngle(const One_Particle_BeamAngle * old);
    Primitive_Observable_Base * Copy() const;
  };

// ---------------------------------------------------------------------------

  class One_Particle_EVis : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_EVis(const String_Matrix &);
    One_Particle_EVis(const One_Particle_EVis * old);
    Primitive_Observable_Base * Copy() const;
    void Evaluate(int nout, const ATOOLS::Vec4D * moms,
		  const ATOOLS::Flavour * flavs,
		  double weight, int ncount);
    void Evaluate(const ATOOLS::Particle_List & plist, 
		  double weight, int ncount);
  };
   
// ---------------------------------------------------------------------------
   
  class One_Particle_Px : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_Px(const String_Matrix &);
    One_Particle_Px(const One_Particle_Px * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Py : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_Py(const String_Matrix &);
    One_Particle_Py(const One_Particle_Py * old);
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Pz : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    One_Particle_Pz(const String_Matrix &);
    One_Particle_Pz(const One_Particle_Pz * old);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
