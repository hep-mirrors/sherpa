#ifndef Six_Particle_Observables_H
#define Six_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Six_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    bool f_special;
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
    Six_Particle_Observable_Base(
        const std::vector<ATOOLS::Flavour> &,
        int,double,double,int,
        const std::string& listname,
        const std::string& name);
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
    Six_Particle_PlaneAngle(const std::vector<ATOOLS::Flavour>&,
			     int,double,double,int,
			     const std::string& listname);//=std::string("NRAngle:"));
    Primitive_Observable_Base* Copy() const;
  };

  class Six_Particle_PT : public Six_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &, 
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
		  double=1.,int=1);
  public:
    Six_Particle_PT(const std::vector<ATOOLS::Flavour>&,
		     int,double,double,int,
		     const std::string& listname);
    Primitive_Observable_Base* Copy() const;
  };

}

#endif
