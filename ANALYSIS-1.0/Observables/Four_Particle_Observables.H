#ifndef Four_Particle_Observables_H
#define Four_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Four_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    bool f_special;
    std::vector<ATOOLS::Flavour> m_flavs;
    void virtual Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
			  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
			  double=1.,int=1) =0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    Four_Particle_Observable_Base(const String_Matrix & parameters, const std::string & obsname );
    Four_Particle_Observable_Base(const Four_Particle_Observable_Base*);
    void virtual Evaluate(int,const ATOOLS::Vec4D *,
			  const ATOOLS::Flavour *,double=1.,int=1);
    void virtual Evaluate(const ATOOLS::Particle_List &,
			  double=1.,int=1);
  };

  // ======================================================================

  class Di_Mass : public Primitive_Observable_Base {
    public:
      Di_Mass(const String_Matrix &);
      Di_Mass(const Di_Mass * old) : Primitive_Observable_Base(*old) { }
      Primitive_Observable_Base * Copy() const { return new Di_Mass(this); }
      void Evaluate(const ATOOLS::Blob_List & blobs,double value, int);
  };

  //==========================================================================

  class Four_Particle_PlaneAngle : public Four_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  double=1.,int=1); 
  public:
    Four_Particle_PlaneAngle(const String_Matrix &);
    Four_Particle_PlaneAngle(const Four_Particle_PlaneAngle * old);
    Primitive_Observable_Base* Copy() const;
  };

  class Four_Particle_PT : public Four_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  double=1.,int=1);
  public:
    Four_Particle_PT(const String_Matrix &);
    Four_Particle_PT(const Four_Particle_PT * old);
    Primitive_Observable_Base* Copy() const;
  };

  class Two_Partonpair_PTdiff : public Four_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
		  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
		  double=1.,int=1);
  public:
    Two_Partonpair_PTdiff(const String_Matrix &);
    Two_Partonpair_PTdiff(const Two_Partonpair_PTdiff * old);
    Primitive_Observable_Base* Copy() const;
  };

  class Two_Partonpair_Theta : public Four_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
		  const ATOOLS::Vec4D&, const ATOOLS::Vec4D&,
		  double=1.,int=1);
  public:
    Two_Partonpair_Theta(const String_Matrix &);
    Two_Partonpair_Theta(const Two_Partonpair_Theta * old);
    Primitive_Observable_Base* Copy() const;
  };

}

#endif
