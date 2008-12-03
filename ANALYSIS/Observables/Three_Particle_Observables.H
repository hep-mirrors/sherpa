#ifndef Three_Particle_Observables_H
#define Three_Particle_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {
  class Three_Particle_Observable_Base : public Primitive_Observable_Base {  
  protected:
    ATOOLS::Flavour      m_flav1,m_flav2,m_flav3;
    void virtual Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
			  const ATOOLS::Vec4D & mom3, double weight, int ncount) = 0; 
  public:
    Three_Particle_Observable_Base(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
				   const ATOOLS::Flavour & flav3, 
				   int type, double xmin, double xmax, int nbins, 
				   const std::string & listname,
				   const std::string & name);
    void virtual Evaluate(const ATOOLS::Particle_List & plist, 
			  double weight, int ncount);
  };
  
  class Three_Particle_PT : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D& mom1, const ATOOLS::Vec4D& mom2, 
                  const ATOOLS::Vec4D& mom3, double weight, int ncount); 
  public:
    Three_Particle_PT(const ATOOLS::Flavour& flav1,
		      const ATOOLS::Flavour& flav2,
		      const ATOOLS::Flavour& flav3, int type, double xmin, 
		      double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base* Copy() const;
  };

  class Three_Particle_Y : public Three_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D& mom1, const ATOOLS::Vec4D& mom2,
                  const ATOOLS::Vec4D& mom3, double weight, int ncount);
  public:
    Three_Particle_Y(const ATOOLS::Flavour& flav1, const ATOOLS::Flavour& flav2,
		     const ATOOLS::Flavour & flav3,
		     int type, double xmin, double xmax, int nbins,
		     const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_DEta : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DEta(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			const ATOOLS::Flavour & flav3,
			int type, double xmin, double xmax, int nbins, 
			const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };
  
  class Three_Particle_DPhi : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DPhi(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			const ATOOLS::Flavour & flav3, int type, double xmin, 
			double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_DR : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_DR(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
		      const ATOOLS::Flavour & flav3, int type, double xmin, 
		      double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_3Mass2 : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_3Mass2(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			const ATOOLS::Flavour & flav3, int type, double xmin, 
			double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_3Mass : public Three_Particle_Observable_Base {  
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
		  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_3Mass(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			 const ATOOLS::Flavour & flav3, int type, double xmin, 
			 double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };

  class Three_Particle_3EnergyCMS : public Three_Particle_Observable_Base {
    void Evaluate(const ATOOLS::Vec4D & mom1, const ATOOLS::Vec4D & mom2, 
                  const ATOOLS::Vec4D & mom3, double weight, int ncount); 
  public:
    Three_Particle_3EnergyCMS(const ATOOLS::Flavour & flav1, const ATOOLS::Flavour & flav2,
			      const ATOOLS::Flavour & flav3, int type, double xmin, 
			      double xmax, int nbins, const std::string & listname);
    Primitive_Observable_Base * Copy() const;
  };
   
}
#endif
