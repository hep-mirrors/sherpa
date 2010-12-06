#ifndef Analysis_Observables_One_Particle_Observables_H
#define Analysis_Observables_One_Particle_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {
  class One_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour m_flav;
    void virtual Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount) = 0; 
    void virtual Evaluate(double value, double weight, double ncount); 
  public:
    One_Particle_Observable_Base(const ATOOLS::Flavour & flav,
				 int type, double xmin, double xmax, int nbins, 
				 const std::string & listname,
				 const std::string & name);
    void virtual Evaluate(int nout, const ATOOLS::Vec4D * moms,
			  const ATOOLS::Flavour * flavs,
			  double weight, double ncount);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, double ncount);
    void virtual EvaluateNLOcontrib(const ATOOLS::Vec4D & mom,
				    double weight, double ncount);
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
  };

  class One_Particle_ET : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_ET(const ATOOLS::Flavour & flav,
		    int type,double xmin,double xmax,int nbins,
		    const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_PT : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_PT(const ATOOLS::Flavour & flav,
		    int type,double xmin,double xmax,int nbins,
		    const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Eta : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_Eta(const ATOOLS::Flavour & flav,
		     int type,double xmin,double xmax,int nbins,
		     const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Y : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_Y(const ATOOLS::Flavour & flav,
		   int type,double xmin,double xmax,int nbins,
		   const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_Phi : public One_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_Phi(const ATOOLS::Flavour & flav,
		     int type,double xmin,double xmax,int nbins,
		     const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_E : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_E(const ATOOLS::Flavour & flav,
		   int type,double xmin,double xmax,int nbins,
		   const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_M : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_M(const ATOOLS::Flavour & flav,
		   int type,double xmin,double xmax,int nbins,
		   const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

  class One_Particle_BeamAngle : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_BeamAngle(const ATOOLS::Flavour & flav,
		   int type,double xmin,double xmax,int nbins,
		   const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
  };

// ---------------------------------------------------------------------------

  class One_Particle_EVis : public One_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom, double weight, double ncount); 
  public:
    One_Particle_EVis(const ATOOLS::Flavour & flav,
		      int type,double xmin,double xmax,int nbins,
		      const std::string & listname=std::string(""));
    Primitive_Observable_Base * Copy() const;
    void Evaluate(int nout, const ATOOLS::Vec4D * moms,
		  const ATOOLS::Flavour * flavs,
		  double weight, double ncount);
    void Evaluate(const ATOOLS::Particle_List & plist, 
		  double weight, double ncount);
  };
}

#endif
