#ifndef  Sample_Analysis_H
#define Sample_Analysis_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class Sample_Analysis  {
    bool        hepevt;
    bool        InitIO;
    int         status,hepevtmode;
    int         nhep;
    std::string m_path, m_file;
    ATOOLS::Primitive_Analysis * ana;
    std::vector<ATOOLS::Primitive_Observable_Base *> obs;
    // ----- test only
  public:
    static int    m_njet;
    static double m_pt_W;
    Sample_Analysis(const std::string & ,const std::string &);
    ~Sample_Analysis();
    void Init();

    void AfterME(ATOOLS::Blob_List *,double);
    void AfterPartonShower(ATOOLS::Blob_List *,double);
    void AfterHadronization(ATOOLS::Blob_List *,double );
    void Finish(const std::string);
  };   
}

#endif
