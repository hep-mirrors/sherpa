#ifndef  Sample_Analysis_H
#define Sample_Analysis_H

#include "Particle_List.H"
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
    bool Treat(ATOOLS::Blob_List *, double &);
    void CleanUp();
  };

  class Sample_Analysis  {
    bool                            hepevt;
    bool                            InitIO;
    int                             status,hepevtmode;
    int                             nhep;
    std::string m_path, m_file;
    ATOOLS::Primitive_Analysis * ana;
    std::vector<ATOOLS::Primitive_Observable_Base *> obs;
  public:
    Sample_Analysis(std::string,std::string);
    ~Sample_Analysis();
    void Init();

    void AfterME(ATOOLS::Blob_List *,double);
    void AfterPartonShower(ATOOLS::Blob_List *,double);
    void AfterHadronization(ATOOLS::Blob_List *,double );
    void Finish(std::string);
  };   
}

#endif
