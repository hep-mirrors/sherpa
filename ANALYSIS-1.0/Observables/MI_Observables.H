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
				     const int nbins,const std::string &listname);
    void Evaluate(const ATOOLS::Blob_List &bloblist,
		  double weight=1.,int ncount=1);
    Primitive_Observable_Base *Copy() const;
    void EndEvaluation(double scale);
  };// end of class Forward_Backward_Eta_Correlation

  class Transversal_PT_Jet: public Primitive_Observable_Base {  
  private:
    double      m_ptmin, m_ptmax,m_phiconesize;
    std::string m_jetlistname;
  public:
    Transversal_EtaPhi(const int type,
		       const double jetptmin,const double jetptmax,
		       const double ptmin,const double ptmax,const double phiconesize,
		       const int nbins,const std::string &listname);
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    Primitive_Observable_Base *Copy() const;
  };

  class Azimuthal_PT: public Primitive_Observable_Base {  
  public:
    Azimuthal_PT(const int type,
		 const double philocation,const double phiconesize,
		 const int nbins,const std::string &listname);
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    Primitive_Observable_Base *Copy() const;
    void EndEvaluation(double scale);
  };

  class Transversal_PT_Particle: public Primitive_Observable_Base {  
  public:
    Transversal_PT(const int type,
		   const double jetptmin,const double jetptmax,
		   const double ptmin,const double ptmax,const double phiconesize,
		   const int nbins,const std::string &listname);
    void Evaluate(const ATOOLS::Particle_List &particlelist,
		  double weight=1.,int ncount=1);
    Primitive_Observable_Base *Copy() const;
    void EndEvaluation(double scale);
  };

}// end of namespace ANALYSIS

#endif
