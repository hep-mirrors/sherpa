#define Sample_Analysis_H
#ifdef  Sample_Analysis_H

#include "IO_HepEvt.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"

namespace SHERPA {
  class Sample_Analysis {
    IO_HepEvt  convert;
    int status;
    APHYTOOLS::Primitive_Analysis * ana;
    std::vector<APHYTOOLS::Primitive_Observable_Base *> obs;
  public:
    Sample_Analysis();
    void Init();
    void AfterME(APHYTOOLS::Blob_List *);
    void AfterPartonShower(APHYTOOLS::Blob_List *);
    void AfterHadronization(APHYTOOLS::Blob_List *);
    void Finish();
  };


}

#endif
