#ifndef Analysis_Observables_Shape_Observables_EE_H
#define Analysis_Observables_Shape_Observables_EE_H

#include "AddOns/Analysis/Observables/Event_Shapes_EE.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Matrix.H"

namespace ANALYSIS {
  class Event_Shapes_Observable_Base : public Primitive_Observable_Base {  
  protected:
    std::string m_key;
  public:
    Event_Shapes_Observable_Base(int,double,double,int,const std::string & ="");
  };

  class Thrust : public Event_Shapes_Observable_Base {
  public:
    Thrust(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="Thrust");
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class Major : public Event_Shapes_Observable_Base {
  public:
    Major(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="Major");
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class Minor : public Event_Shapes_Observable_Base {
  public:
    Minor(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="Minor");
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class Oblateness : public Event_Shapes_Observable_Base {
  public:
    Oblateness(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="Oblateness");
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    Primitive_Observable_Base * Copy() const;
  };

  class PT_In_Thrust : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List & pl,double weight,double ncount);
  public:
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    PT_In_Thrust(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="PT_In_Thrust");
    Primitive_Observable_Base * Copy() const;
  };

  class PT_Out_Thrust : public Event_Shapes_Observable_Base {
  protected:
    void Evaluate(const ATOOLS::Particle_List & pl,double weight,double ncount);
  public:
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
    PT_Out_Thrust(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="PT_Out_Thrust");
    Primitive_Observable_Base * Copy() const;
  };

  class Eta_Thrust : public Event_Shapes_Observable_Base {
  public:
    Eta_Thrust(int type, double xmin, double xmax, int nbin, 
	      const std::string & listname, const std::string & name ="EtaThrust");
    void Evaluate(const ATOOLS::Blob_List & ,double weight, double ncount);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
