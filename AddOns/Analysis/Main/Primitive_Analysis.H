#ifndef Analysis_Main_Primitive_Analysis_H
#define Analysis_Main_Primitive_Analysis_H

#include <vector>

#include "ATOOLS/Phys/Particle_List.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Org/File_IO_Base.H"

namespace ATOOLS {
  class Particle_Qualifier_Base;
}

namespace ANALYSIS {

  const static std::string finalstate_list("FinalState");

  enum code {

    splitt_jetseeds = 2,
    splitt_process  = 4,
    splitt_all      = 6,
    fill_histos     = 8,
    fill_helper     = 16,
    fill_all        = 24,
    
    output_this     = 128,
    output_jet      = 256,
    output_process  = 512,
    
    do_me           = 1024,
    do_mi           = 2048,
    do_shower       = 4096,
    do_hadron       = 8192,
    splitt_phase    = 16384,
    splitt_extra    = 32768,
    do_menlo        = 65536
  };

  class Analysis_Object;
  class Analysis_Handler;
  class Primitive_Analysis;

  typedef std::vector<Analysis_Object*>                 Object_List;
  typedef std::map<std::string,ATOOLS::Particle_List*>  PL_Container;
  typedef std::map<std::string,Primitive_Analysis *>    Analysis_List;

  class Primitive_Analysis: public ATOOLS::File_IO_Base {
  private:

    int                    m_mode;
    long                   m_nevt;
    std::string            m_name, m_maxjettag;

    Object_List            m_objects;
    PL_Container           m_pls;

    std::set<Primitive_Analysis*> m_called;

    ATOOLS::String_BlobDataBase_Map m_datacontainer;
    Analysis_List          m_subanalyses;

    // reference to the event record
    const ATOOLS::Blob_List    * p_blobs;

    Primitive_Analysis   * p_partner;

    bool  m_active, m_splitjetconts;

    Analysis_Handler *p_ana;

    std::string JetID(std::string name,std::string max) const;
    
    void Init();
    bool SelectBlob(const ATOOLS::Blob *blob); 
    void CreateFinalStateParticleList();

    Primitive_Analysis * GetSubAnalysis(const ATOOLS::Blob_List *const bl,
					const std::string & key, int mode);
    void CallSubAnalysis(const ATOOLS::Blob_List * const bl, double value);

    void PrintStatus();

  public :

    Primitive_Analysis(Analysis_Handler *const ana,const int);
    Primitive_Analysis(Analysis_Handler *const ana,const std::string, const int);

    ~Primitive_Analysis();
    
    void AddSubAnalysis(const std::string &,Primitive_Analysis *);
    void DoAnalysis(const ATOOLS::Blob_List * const , const double);
    bool DoAnalysisNLO(const ATOOLS::Blob_List * const , const double);
    void FinishAnalysis(const std::string &);
    void RestoreAnalysis();
    
    void AddObject(Analysis_Object *obj);
    Analysis_Object *GetObject(const std::string & key);

    void AddParticleList(const std::string &,ATOOLS::Particle_List *);
    ATOOLS::Particle_List * GetParticleList(const std::string &,
					    const bool=false);

    ATOOLS::Blob_Data_Base * operator[](const std::string name);
    void AddData(const std::string name,ATOOLS::Blob_Data_Base * data); 

    void ClearAllData();

    void SetPartner(Primitive_Analysis * const);
    void SetBlobType(const std::string &);
    int  NumberOfObjects() { return m_objects.size(); }

    void Test(const int mode=0);
    
    // inline methods
    inline void SetMaxJetTag(const std::string &tag) { m_maxjettag=tag; }
    inline void SetSplitJetConts(bool tag) { m_splitjetconts=tag; }

    inline const std::string &Name() const { return m_name; }
    inline int                Mode() const { return m_mode; }

    inline const ATOOLS::String_BlobDataBase_Map &GetData()
    { return m_datacontainer; }

    inline Analysis_Handler *AnalysisHandler() const { return p_ana; }

  };// end of class Primitive_Analysis

}// end of namespace ANALYSIS

#endif
