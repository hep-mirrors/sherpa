#ifndef MI_Observables_H
#define MI_Observables_H

#include "Primitive_Observable_Base.H"
#include "../../AMISIC++-1.0/Main/Amisic_Histogram.H"

namespace ANALYSIS {


  class MI_Statistics: public Primitive_Observable_Base {  
  public:

    MI_Statistics(const std::string &listname,int mode=0);

    void Evaluate(const ATOOLS::Blob_List &blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

  };// end of class MI_Statistics

  class Forward_Backward_Eta_Correlation: public Primitive_Observable_Base {  
  private:

    AMISIC::Amisic_Histogram<double> m_etafw, m_etafwsq, m_etafwbw;

  public:

    Forward_Backward_Eta_Correlation(const int type,
				     const double detamin,const double detamax,
				     const int nbins,const std::string &jetlist,
				     const std::string &listname);

    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Forward_Backward_Eta_Correlation

  class Leading_PT_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

    size_t m_item;

  public:

    Leading_PT_Selector(const size_t item,const double ptmin,const double ptmax,
			const std::string &jetlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class Leading_PT_Selector

  class Leading_ET_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_outlist;

    size_t m_item;

  public:

    Leading_ET_Selector(const size_t item,const double ptmin,const double ptmax,
			const std::string &jetlist,const std::string &outlist);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;
    
    void EndEvaluation(double scale);

  };// end of class Leading_ET_Selector

  class Transverse_Region_Selector: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist, m_outlist;

    ATOOLS::Vec4D m_leadingjet;

  public:

    Transverse_Region_Selector(const double phimin,const double phimax,
			       const std::string &jetlist,const std::string &inlist,
			       const std::string &outlist);

    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Transverse_Region_Selector

  class Multiplicity_vs_JetPT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;
    
  public:

    Multiplicity_vs_JetPT(const int type,
			  const double ptmin,const double ptmax,
			  const int nbins,const std::string &jetlist,
			  const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_JetPT

  class Multiplicity_vs_PT: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

  public:

    Multiplicity_vs_PT(const int type,
		       const double ptmin,const double ptmax,
		       const int nbins,const std::string &jetlist,
		       const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;

  };// end of class Multiplicity_vs_PT

  class Multiplicity_vs_DPhi: public Primitive_Observable_Base {  
  private:

    std::string m_jetlist;

    AMISIC::Amisic_Histogram<double> m_histogram;
    
  public:

    Multiplicity_vs_DPhi(const int type,
			 const double dphimin,const double dphimax,
			 const int nbins,const std::string &jetlist,
			 const std::string &listname);
    
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    
    Primitive_Observable_Base *Copy() const;

    void EndEvaluation(double scale);

  };// end of class Multiplicity_vs_DPhi

}// end of namespace ANALYSIS

#endif
