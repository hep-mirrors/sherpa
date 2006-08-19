#ifndef Shower_Statistics_H
#define Shower_Statistics_H

#include "Primitive_Observable_Base.H"
#include <iostream>

namespace ANALYSIS {

  class Shower_Statistics : public Primitive_Observable_Base {  
  protected:

    std::vector<ATOOLS::Histogram*> m_observables;

  public:

    Shower_Statistics(const std::string & inlistname, int mode=0);

    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base * Copy() const;

    void EndEvaluation(double);
    virtual void Output(const std::string & pname);

    Primitive_Observable_Base & operator+=(const Primitive_Observable_Base &);

  };

}

#endif
