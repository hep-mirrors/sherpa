#ifndef  Sample_Analysis_H
#define Sample_Analysis_H

#include "Particle_List.H"
#include "Blob_List.H"
#include "Primitive_Analysis.H"
#include "Event_Phase_Handler.H"

namespace SHERPA {
  class PHard_Observable : public ATOOLS::Primitive_Observable_Base {  
  public:
    PHard_Observable(int _type,double _xmin,double _xmax,int _nbins, 
		     std::string _name);
    void Evaluate(double); 
    void Evaluate(int,ATOOLS::Vec4D *,ATOOLS::Flavour *,double w);
    void Evaluate(const ATOOLS::Particle_List &,double w);
    void Evaluate(const ATOOLS::Blob_List &,double w); 
  };
  
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
