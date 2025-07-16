#ifndef Analysis_Observables_PSM_Observables_H
#define Analysis_Observables_PSM_Observables_H

#include "AddOns/Analysis/Observables/Primitive_Observable_Base.H"

namespace ANALYSIS {



  class PSM_Observable : public Primitive_Observable_Base {
  protected:
    std::vector<int> m_pnb;
  public:
    PSM_Observable(unsigned int type,double xmin,double xmax,int nbins,
		   int p0, int p1, int p2, int p3, 
		   const std::string & =std::string("PSM"));

    void Evaluate(const ATOOLS::Particle_List &,double, double);
    void Evaluate(const ATOOLS::Blob_List & blobs,double value, double);
    
    void EndEvaluation(double scale=1.0);
    void Restore(double scale=1.0);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);
    void Reset();
    Primitive_Observable_Base * Copy() const;

  };

}
#endif
