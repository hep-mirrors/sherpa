#ifndef All_Particle_Selector_H
#define All_Particle_Selector_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class PT_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

  public:

    PT_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class PT_Selector

  class ET_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

  public:

    ET_Selector(const double min,const double max,
		const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class ET_Selector

  class Eta_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

  public:

    Eta_Selector(const double min,const double max,
		 const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class Eta_Selector

  class Y_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

  public:

    Y_Selector(const double min,const double max,
		 const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class Y_Selector

  class Phi_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

  public:

    Phi_Selector(const double min,const double max,
		 const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class Phi_Selector

  class DPhi_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_reflist, m_outlist;
    size_t      m_item;

    ATOOLS::Flavour m_flavour;

  public:

    DPhi_Selector(const double min,const double max,
		  const ATOOLS::Flavour flav,
		  const size_t item,const std::string &reflist,
		  const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class DPhi_Selector

}// end of namespace ANALYSIS

#endif
