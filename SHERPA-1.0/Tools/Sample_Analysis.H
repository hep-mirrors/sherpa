#ifndef  Sample_Analysis_H
#define Sample_Analysis_H

#include "IO_HepEvt.H"
#include "Parton_List.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class Sample_Analysis;

  class Analysis_Phase : public Event_Phase_Handler {
    Sample_Analysis * p_analysis;
    int m_mode;
    int m_status;
  public:
    Analysis_Phase(Sample_Analysis *,int);
    bool Treat(APHYTOOLS::Blob_List *, double &);
    void CleanUp();
  };

  class Sample_Analysis  {
    bool                            hepevt;
    bool                            InitIO;
    int                             status,hepevtmode;
    int                             nhep;
    IO_HepEvt                     * convert;
    APHYTOOLS::Primitive_Analysis * ana;
    std::vector<APHYTOOLS::Primitive_Observable_Base *> obs;
  public:
    Sample_Analysis(IO_HepEvt *&,bool);
    ~Sample_Analysis();
    void Init();

    void AfterME(APHYTOOLS::Blob_List *,double);
    void AfterPartonShower(APHYTOOLS::Blob_List *,double);
    void AfterHadronization(APHYTOOLS::Blob_List *,double );
    void Finish();
  };   
}

#endif
