#ifndef Primitive_Analysis_H
#define Primitive_Analysis_H

#include <vector>

#include "Particle_List.H"
#include "Blob_List.H"

namespace ATOOLS {
  namespace ANALYSIS {
    enum code {
      unweighted = 0,
      weighted   = 1,
      splitt_jetseeds=2,
      splitt_process =4,
      splitt_all = 6,
      fill_histos= 8,
      fill_helper=16,
      fill_all   =24,

      norm       =96,
      norm_to_xs =32,
      norm_to_one=64,

      output_this=128,
      output_jet =256,
      output_process=512,

      do_me = 1024,
      do_shower = 2048,
      do_hadron = 4096,
      splitt_phase = 8192
    };
  }

  class Primitive_Observable_Base;

  typedef std::vector<Primitive_Observable_Base *> Observable_List;
  typedef std::map<std::string,Particle_List *>    PL_Container;

  class Primitive_Analysis;
  typedef std::map<std::string,Primitive_Analysis *> Analysis_List;

  class Primitive_Analysis {
  private:
    int   m_mode;
    long  m_nevt;
    std::string       m_name;
    Observable_List   m_observables;

    PL_Container      m_pls;
    Data_Container    m_datacontainer;
    Analysis_List     m_subanalyses;

    // reference to the event record
    Blob_List *       p_blobs;

    Primitive_Analysis * p_partner;

    void Init();
    void CreateFinalStateParticleList();
    void CreateChargedParticleList();

    Primitive_Analysis * GetSubAnalysis(const std::string & key, int mode);
    void CallSubAnalysis(Blob_List * const bl, double value);

    void PrintStatus();
  public :
    Primitive_Analysis(const std::string, const int);
    Primitive_Analysis(const int);
    ~Primitive_Analysis();

    void AddObservable(Primitive_Observable_Base *);
    void AddSubAnalysis(const std::string &,Primitive_Analysis *);
    // most simple analysis, assume the observables a know how to handle a single value
    void DoAnalysis(const Particle_List &, double);

    // advanced analysis, working with a bloblist, thus event generation properties can be exploit
 
    void DoAnalysis(Blob_List * const , double);
    //void DoAnalysis(const Blob_List &, double);
 
   // finialise and write out results of analysis
    void FinishAnalysis(const std::string &, long=0, double =1.);
    


    Primitive_Observable_Base * GetObservable(const std::string &);

    Particle_List * GetParticleList(const std::string &);
    void AddParticleList(const std::string &,Particle_List *);

    Blob_Data_Base * operator[](const std::string name) 
    {
      Data_Container::const_iterator cit=m_datacontainer.find(name);
      if (cit==m_datacontainer.end()) return 0;
      return cit->second;
    } 
    void AddData(const std::string name, Blob_Data_Base * data); 
    void ClearAllData();

    void SetPartner(Primitive_Analysis * const);
    void SetBlobType(const std::string &);
  };
}

#endif
