#ifndef Shape_Observables_EE_H
#define Shape_Observables_EE_H

#include "Event_Shapes_EE.H"
#include "Primitive_Observable_Base.H"
#include "Vector.H"
#include "Matrix.H"

namespace ANALYSIS {
  class Event_Shapes_Observable_Base : public Primitive_Observable_Base {  
  protected:
    std::string m_key,m_inlistname;
  public:
    Event_Shapes_Observable_Base(int,double,double,int,const std::string & =std::string(""));
  };

  class Thrust : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    Thrust(int,double,double,int,const std::string & =std::string("Thrust"));
    Primitive_Observable_Base * Copy() const;
  };

  class Major : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    Major(int,double,double,int,const std::string & =std::string("Major"));
    Primitive_Observable_Base * Copy() const;
  };

  class Minor : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    Minor(int,double,double,int,const std::string & =std::string("Minor"));
    Primitive_Observable_Base * Copy() const;
  };

  class Oblateness : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    Oblateness(int,double,double,int,const std::string & =std::string("Oblateness"));
    Primitive_Observable_Base * Copy() const;
  };

  class PT_In_Thrust : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    PT_In_Thrust(int,double,double,int,const std::string & _inlistname,
		 const std::string & =std::string("PT_In_Thrust"));
    Primitive_Observable_Base * Copy() const;
  };

  class PT_Out_Thrust : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List &,double=1.,int ncount=1);
  public:
    PT_Out_Thrust(int,double,double,int,const std::string & _inlistname,
		  const std::string & =std::string("PT_Out_Thrust"));
    Primitive_Observable_Base * Copy() const;
  };
}

#endif
