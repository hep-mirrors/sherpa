#ifndef Four_Particle_Observables_H
#define Four_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ATOOLS {
  class Four_Particle_Observable_Base : public ATOOLS::Primitive_Observable_Base {  
  protected:
    std::vector<Flavour> m_flavs;
    void virtual Evaluate(const Vec4D &, const Vec4D &, 
			  const Vec4D &, const Vec4D &,double=1.,int=1) =0; 
    void virtual Evaluate(double,double=1.,int=1); 
  public:
    Four_Particle_Observable_Base(const std::vector<Flavour> &,int,double,double,int,
				  const std::string & =std::string(""));
    void virtual Evaluate(int,const Vec4D *,const Flavour *,double=1.,int=1);
    void virtual Evaluate(const Particle_List &,double=1.,int=1);
  };

  class Four_Particle_PlaneAngle : public ATOOLS::Four_Particle_Observable_Base {  
    void Evaluate(const Vec4D &, const Vec4D &, 
		  const Vec4D &, const Vec4D &,double=1.,int=1); 
  public:
    Four_Particle_PlaneAngle(const std::vector<Flavour> &,int,double,double,int,
			     const std::string & =std::string("NRAngle:"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
