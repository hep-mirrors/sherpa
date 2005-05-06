#ifndef One_Particle_Selector_H
#define One_Particle_Selector_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class One_PT_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_outlist;
    ATOOLS::Flavour m_flavour;
    
    size_t       m_item;
    unsigned int m_mode;
    
  public:

    One_PT_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_PT_Selector

  class One_ET_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_outlist;
    ATOOLS::Flavour m_flavour;

    size_t m_item;
    unsigned int m_mode;
  public:

    One_ET_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_ET_Selector

  class One_Eta_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_outlist;
    ATOOLS::Flavour m_flavour;

    size_t m_item;
    unsigned int m_mode;
  public:

    One_Eta_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_Eta_Selector

  class One_Y_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_outlist;
    ATOOLS::Flavour m_flavour;

    size_t m_item;
    unsigned int m_mode;
  public:

    One_Y_Selector(const ATOOLS::Flavour flav,
		    const size_t item,const int mode,const double min,const double max,
		    const std::string &inlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_Y_Selector

  class One_DPhi_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_reflist, m_outlist;
    ATOOLS::Flavour m_flavour, m_refflavour;

    size_t m_item, m_refitem;

  public:

    One_DPhi_Selector(const ATOOLS::Flavour flav,const size_t item,
		      const ATOOLS::Flavour ref,const size_t refitem,
		      const double min,const double max,
		      const std::string &inlist,const std::string &reflist,
		      const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_DPhi_Selector

  class One_ETFrac_Selector: public Primitive_Observable_Base {  
  private:

    std::string     m_reflist, m_outlist;
    ATOOLS::Flavour m_flavour, m_refflavour;

    size_t m_item, m_refitem;

  public:

    One_ETFrac_Selector(const ATOOLS::Flavour flav,const size_t item,
		       const ATOOLS::Flavour ref,const size_t refitem,
		       const double min,const double max,
		       const std::string &inlist,const std::string &reflist,
		       const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class One_ETFrac_Selector

}// end of namespace ANALYSIS

#endif
