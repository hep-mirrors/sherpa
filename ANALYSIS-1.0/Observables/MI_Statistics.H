#ifndef MI_Statistics_H
#define MI_Statistics_H

#include "Primitive_Observable_Base.H"
#include <iostream>

namespace ANALYSIS {

  class MI_Statistics : public Primitive_Observable_Base {  
  public:

    MI_Statistics(const std::string & inlistname, int mode=0);

    void Evaluate(const ATOOLS::Blob_List & blobs,double weight=1.,int ncount=1);

    Primitive_Observable_Base * Copy() const;

  };

}

#endif
