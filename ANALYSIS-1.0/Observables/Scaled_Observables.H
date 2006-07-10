#ifndef Scaled_Observables_H
#define Scaled_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Scaled_Observable_Base : public Primitive_Observable_Base {  
  protected:
    double m_ecms;
    void virtual Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount) = 0; 
    void virtual Evaluate(double value, double weight, int ncount); 
  public:
  Scaled_Observable_Base(const String_Matrix & parameters, const std::string & obsname);
  Scaled_Observable_Base(const Scaled_Observable_Base*);
//     Scaled_Observable_Base(const String_Matrix &,const Primitive_Observable_Defaults *,
// 			   const std::string, const double);
//     Scaled_Observable_Base(int type, double xmin, double xmax, int nbins, 
// 			   const std::string & listname,
// 			   const std::string & name, const double ecms);
    void virtual Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  double weight, int ncount);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };

  class Scaled_Momentum : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Scaled_Momentum(const String_Matrix &);
    Scaled_Momentum(const Scaled_Momentum *);
    Primitive_Observable_Base * Copy() const;
  };

  class Log_Scaled_Momentum : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Log_Scaled_Momentum(const String_Matrix &);
    Log_Scaled_Momentum(const Log_Scaled_Momentum *);
    Primitive_Observable_Base * Copy() const;
  };

  class Scaled_Energy : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Scaled_Energy(const String_Matrix &);
    Scaled_Energy(const Scaled_Energy *);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
