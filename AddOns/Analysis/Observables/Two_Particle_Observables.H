#ifndef Analysis_Observables_Two_Particle_Observables_H
#define Analysis_Observables_Two_Particle_Observables_H

#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Two_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour      m_flav1,m_flav2;
    void virtual Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
			  double weight, double ncount) = 0; 
    //    void virtual Evaluate(double value, double weight, double ncount); 
  public:
    Two_Particle_Observable_Base(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
				 int type, double xmin, double xmax, int nbins, 
				 const std::string & listname,
				 const std::string & name);
//     void virtual Evaluate(int nout, const ATOOLS::Vec4D * moms,
// 			  const ATOOLS::Flavour * flavs,
// 			  double weight, double ncount);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, double ncount);
    void virtual EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
			    double weight, double ncount); 
    void EvaluateNLOcontrib(double weight, double ncount);
    void EvaluateNLOevt();
  };

  class Two_Particle_Mass : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Mass(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_PT(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		    int type, double xmin, double xmax, int nbins, 
		    const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Scalar_PT : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Scalar_PT(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			   int type, double xmin, double xmax, int nbins, 
			   const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

#ifdef USING__ROOT
  class Two_Particle_Angles : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D &,const ATOOLS::Vec4D &,double=1., double=1.); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D &, const ATOOLS::Vec4D &,
      double=1., double=1.);
  public:
    Two_Particle_Angles(const ATOOLS::Flavour &,const ATOOLS::Flavour &,
			int,double,double,int,const std::string & =std::string("SP_T:"),const std::string & =std::string("FinalState"));
    Primitive_Observable_Base * Copy() const;
  };
#endif

  class Two_Particle_Eta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Eta(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		     int type, double xmin, double xmax, int nbins, 
		     const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DEta : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_DEta(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Y : public Two_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2,
                  double weight, double ncount);
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Y(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		   int type, double xmin, double xmax, int nbins,
		   const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DY : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_DY(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Phi : public Two_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Phi(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		     int type, double xmin, double xmax, int nbins, 
		     const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DPhi : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_DPhi(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_Angle : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Angle(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_DR : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_DR(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		    int type, double xmin, double xmax, int nbins, 
		    const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_CMS_Angle : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_CMS_Angle(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			   int type, double xmin, double xmax, int nbins, 
			   const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };
   
  class Two_Particle_Mass2 : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_Mass2(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      int type, double xmin, double xmax, int nbins, 
		      const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Two_Particle_MT2 : public Two_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
    void EvaluateNLOcontrib(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  double weight, double ncount); 
  public:
    Two_Particle_MT2(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		     int type, double xmin, double xmax, int nbins, 
		     const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

}

#endif
