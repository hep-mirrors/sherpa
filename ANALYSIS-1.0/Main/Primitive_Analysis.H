#ifndef Primitive_Analysis_H
#define Primitive_Analysis_H

#include <vector>

#include "Particle_List.H"
#include "Blob_List.H"
#include "Blob.H"
#include "File_IO_Base.H"

namespace ATOOLS {
  class Particle_Qualifier_Base;
}

namespace ANALYSIS {
  enum code {
    unweighted      = 0,
    weighted        = 1,
    splitt_jetseeds = 2,
    splitt_process  = 4,
    splitt_all      = 6,
    fill_histos     = 8,
    fill_helper     = 16,
    fill_all        = 24,
    
    norm_to_xs      = 32,
    norm_to_one     = 64,
    norm            = 96,
    
    output_this     = 128,
    output_jet      = 256,
    output_process  = 512,
    
    do_me           = 1024,
    do_mi           = 2048,
    do_shower       = 4096,
    do_hadron       = 8192,
    splitt_phase    = 16384,
    splitt_extra    = 32768
  };

  class Primitive_Observable_Base;
  class Primitive_Analysis;

  typedef std::vector<Primitive_Observable_Base *>      Observable_List;
  typedef std::map<std::string,ATOOLS::Particle_List *> PL_Container;
  typedef std::map<std::string,Primitive_Analysis *>    Analysis_List;

  struct Weight_Statistics {
    double sum_weight;
    unsigned long nevt;
    double sum_weight_one;
    unsigned long nevt_one;
    Weight_Statistics() :
      sum_weight(0.),nevt(0),sum_weight_one(0.),nevt_one(0) {}
  };


  class Primitive_Analysis: public ATOOLS::File_IO_Base {
  private:
    int                    m_mode;
    long                   m_nevt;
    std::string            m_name;

    Observable_List        m_observables;

    PL_Container           m_pls;
    ATOOLS::String_BlobDataBase_Map m_datacontainer;
    Analysis_List          m_subanalyses;

    // reference to the event record
    const ATOOLS::Blob_List    * p_blobs;

    Primitive_Analysis   * p_partner;

    Weight_Statistics      m_stats;

    ATOOLS::Particle_Qualifier_Base * p_bfinder;
    bool  m_active;
    
    void Init();
    bool SelectBlob(const ATOOLS::Blob *blob); 
    void CreateFinalStateParticleList(bool markb=false);
    void CreateChargedParticleList();
    void CreateIntermediateHadronsList();

    Primitive_Analysis * GetSubAnalysis(const std::string & key, int mode);
    void CallSubAnalysis(const ATOOLS::Blob_List * const bl, double value);

    void PrintStatus();
  public :
    Primitive_Analysis(const std::string, const int);
    Primitive_Analysis(const int);
    ~Primitive_Analysis();

    void AddObservable(Primitive_Observable_Base *);
    void AddSubAnalysis(const std::string &,Primitive_Analysis *);
    // most simple analysis, assume the observables a know how to handle a single value
    void DoAnalysis(const ATOOLS::Particle_List &, const double);
    // advanced analysis, working with a bloblist, thus event generation properties can be exploited
    void DoAnalysis(const ATOOLS::Blob_List * const , const double);
    // finalise and write out results of analysis
    void FinishAnalysis(const std::string &, long=0, double =1.);
    


    Primitive_Observable_Base * GetObservable(const std::string &);

    ATOOLS::Particle_List * GetParticleList(const std::string &);
    void AddParticleList(const std::string &,ATOOLS::Particle_List *);

    ATOOLS::Blob_Data_Base * operator[](const std::string name);
    void AddData(const std::string name,ATOOLS::Blob_Data_Base * data); 
    void ClearAllData();

    void SetPartner(Primitive_Analysis * const);
    void SetBlobType(const std::string &);
    int  NumberOfObservables() { return m_observables.size(); }
    
    // inline methods
    inline const std::string &Name() const { return m_name; }
    inline const int          Mode() const { return m_mode; }
  };
}

#endif
