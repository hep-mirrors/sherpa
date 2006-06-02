#ifndef Primitive_Observable_h
#define Primitive_Observable_h

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class Multiplicity : public Primitive_Observable_Base {  
  public:
    Multiplicity(const std::string &,const std::string &,
		 const std::string & =std::string(""));
    Multiplicity(const std::string &,int,double,double,int,const std::string & =std::string(""));
    Multiplicity(int type,double xmin,double xmax,int nbins,
		 const std::string & listname=std::string(""));
    Multiplicity(ATOOLS::Histogram_Base * histo,
		 const std::string & listname=std::string(""));

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class Multiplicity

  class Hadron_Multiplicities : public Primitive_Observable_Base {  
  public:
    Hadron_Multiplicities(int type,double xmin,double xmax,int nbins,
			  const std::string & listname=std::string(""));
    Hadron_Multiplicities(ATOOLS::Histogram_Base * histo,
			  const std::string & listname=std::string(""));
    Hadron_Multiplicities(const std::string &,const std::string &,
			  const std::string & =std::string(""));
    Hadron_Multiplicities(const std::string &,int,double,double,int,
			  const std::string & =std::string(""));

    void Evaluate(const ATOOLS::Particle_List & pl, double weight, int ncount);
    Primitive_Observable_Base * Copy() const;

  };// end of class Hadron_Multiplicities
}// end of namespace ANALYSIS

#endif
