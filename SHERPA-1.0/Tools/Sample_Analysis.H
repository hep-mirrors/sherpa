#ifndef  Sample_Analysis_H
#define Sample_Analysis_H

#include "IO_HepEvt.H"
#include "Parton_List.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"

namespace SHERPA {
  class Sample_Analysis {
    bool                            hepevt;
    bool                            InitIO;
    int                             status,hepevtmode;
    int                             nhep;
    IO_HepEvt                     * convert;
    APHYTOOLS::Primitive_Analysis * ana;
  public:
    Sample_Analysis(IO_HepEvt *&,bool);
    ~Sample_Analysis();
    void Init();
    void DoAnalysis(APHYTOOLS::Parton_List *) {}
    void AfterME(APHYTOOLS::Blob_List *);
    void AfterPartonShower(APHYTOOLS::Blob_List *);
    void AfterHadronization(APHYTOOLS::Blob_List *);
    void Finish();
  };   
}

#endif
