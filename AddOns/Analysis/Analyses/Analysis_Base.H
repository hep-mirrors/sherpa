#ifndef AddOns_Analysis_Analyses_Analysis_Base_H
#define AddOns_Analysis_Analyses_Analysis_Base_H

#include "AddOns/Analysis/Observables/Normalized_Observable.H"

namespace ANALYSIS {

  class Analysis_Base: public Primitive_Observable_Base {  
  protected:

    std::vector<Normalized_Observable*> m_dists;
    std::vector<ATOOLS::Histogram*>     m_histos;

    void AddZeroPoint(const double &ntrial,const int &mode);

  public:

    Analysis_Base(const std::string &listname);

    ~Analysis_Base();
    
    virtual void Evaluate(double weight,double ncount,int mode) = 0;

    void FillHisto(const size_t &i,const double &x,const double &weight,
		   const double &ntrial,const int &mode);
    void FillDist(const size_t &i,const double &x,const double &y,
		  const double &weight,const double &ntrial,const int &mode);

    void Evaluate(const ATOOLS::Particle_List &pl,
		  double weight, double ncount);
    void EvaluateNLOcontrib(double weight,double ncount);
    void EvaluateNLOevt();

    Analysis_Object &operator+=(const Analysis_Object &obj);

    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    void Output(const std::string & pname);

  };// end of class Analysis_Base

}// end of namespace ANALYSIS

#define AddZero(ntrial,mode) { AddZeroPoint(ntrial,mode); return; }

#endif
