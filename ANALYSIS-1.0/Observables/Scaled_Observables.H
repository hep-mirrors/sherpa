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
    Scaled_Observable_Base(int type, double xmin, double xmax, int nbins, 
			   const std::string & listname,
			   const std::string & name, double ecms);
    void virtual Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  double weight, int ncount);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };

  class Scaled_Momentum : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Scaled_Momentum(int type,double xmin,double xmax,int nbins,
		    const std::string & listname, double ecms);
    Primitive_Observable_Base * Copy() const;
  };

  class Log_Scaled_Momentum : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Log_Scaled_Momentum(int type,double xmin,double xmax,int nbins,
		    const std::string & listname, double ecms);
    Primitive_Observable_Base * Copy() const;
  };

  class Scaled_Energy : public Scaled_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, int ncount); 
  public:
    Scaled_Energy(int type,double xmin,double xmax,int nbins,
		  const std::string & listname, double ecms);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
