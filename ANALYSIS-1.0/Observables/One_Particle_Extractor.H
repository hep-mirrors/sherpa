#ifndef One_Particle_Extractor_H
#define One_Particle_Extractor_H

#include "Primitive_Observable_Base.H"
#include "Variable.H"

namespace ANALYSIS {

  class One_Particle_Extractor: public Primitive_Observable_Base {  
  private:

    ATOOLS::Variable_Base<double> *p_variable;

    std::string     m_reflist, m_outlist;
    ATOOLS::Flavour m_flavour;

    size_t m_item;

  public:

    One_Particle_Extractor(const std::string &type,const ATOOLS::Flavour flav,
			   const size_t item,const double min,const double max,
			   const std::string &inlist,const std::string &reflist,
			   const std::string &outlist);
    
    ~One_Particle_Extractor();

    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_Particle_Extractor

}// end of namespace ANALYSIS

#endif
