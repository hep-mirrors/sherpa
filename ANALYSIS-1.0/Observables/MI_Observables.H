#ifndef MI_Observables_H
#define MI_Observables_H

#include "Primitive_Observable_Base.H"

namespace ANALYSIS {

  class MI_Statistics: public Primitive_Observable_Base {  
  public:

    MI_Statistics(const std::string &listname,int mode=0);

    void Evaluate(const ATOOLS::Blob_List &blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

  };// end of class MI_Statistics

  class Forward_Backward_Eta_Correlation: public Primitive_Observable_Base {  
  private:

    ATOOLS::Histogram m_etafw, m_etabw;

  public:

    Forward_Backward_Eta_Correlation(const int type,
				     const double detamin,const double detamax,
				     const int nbins,const std::string &listname);

    void Evaluate(const ATOOLS::Particle_List &particles,
		  double weight=1.,int ncount=1);

    Primitive_Observable_Base *Copy() const;

  };// end of class Forward_Backward_Eta_Correlation

}// end of namespace ANALYSIS

#endif
