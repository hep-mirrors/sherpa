#ifndef ISR_Statistics_H
#define ISR_Statistics_H

#include "Primitive_Observable_Base.H"
#include <iostream>

namespace ANALYSIS {

  class ISR_Statistics : public Primitive_Observable_Base {  
  protected:

    std::vector<ATOOLS::Histogram*> m_observables;

  public:

    ISR_Statistics(const std::string & inlistname, int mode=0);

    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base * Copy() const;

    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);

  };// end of class ISR_Statistics

#ifdef USING__ROOT
  class Sprime_Y_Distribution : public Primitive_Observable_Base {  
  private:

    double m_ymin, m_ymax;
    size_t m_ybins;

  public:

    Sprime_Y_Distribution(const double spmin,const double spmax,const size_t spbins,
			  const double ymin,const double ymax,const size_t ybins);

    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base * Copy() const;

  };// end of class Sprime_Y_Distribution
#endif

}// end of namespace ANALYSIS

#endif
