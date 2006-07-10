#ifndef Six_Particle_Observables_H
#define Six_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Six_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    std::vector<ATOOLS::Flavour> m_flavs;
    void virtual Evaluate(
        const ATOOLS::Vec4D &, 
        const ATOOLS::Vec4D &, 
        const ATOOLS::Vec4D &, 
        const ATOOLS::Vec4D &, 
        const ATOOLS::Vec4D &, 
        const ATOOLS::Vec4D &,
        double=1.,int=1) =0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    Six_Particle_Observable_Base( const String_Matrix & parameters, const std::string & obsname );
    Six_Particle_Observable_Base(const Six_Particle_Observable_Base*);
    void virtual Evaluate(int,const ATOOLS::Vec4D *,
			  const ATOOLS::Flavour *,double=1.,int=1);
    void virtual Evaluate(const ATOOLS::Particle_List &,
			  double=1.,int=1);
  };

  //==========================================================================

  class Six_Particle_PlaneAngle : public Six_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  double=1.,int=1); 
  public:
    Six_Particle_PlaneAngle(const String_Matrix & parameters);
    Six_Particle_PlaneAngle(const Six_Particle_PlaneAngle * old);
    Primitive_Observable_Base* Copy() const;
  };

}

#endif
