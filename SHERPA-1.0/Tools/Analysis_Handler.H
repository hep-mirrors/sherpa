#ifndef Analysis_Handler_H
#define Analysis_Handler_H

#include "Exception.H"
#include "Primitive_Analysis.H"
#include "Final_Selector.H"
#include "Particle_List.H"
#include "Blob_List.H"
#include <iostream>
#include <map>

namespace SHERPA {
  class Observable_Data {
  public:
    Observable_Data(std::string=std::string(""));
    std::string              type;
    std::vector<int>         ints;
    std::vector<double>      numbers;
    std::vector<std::string> keywords;
    int  Specify();
    void Output();
  };

  class Analysis_Handler: public ATOOLS::Terminator_Object  {
    std::string  m_phase, m_outputpath;
    std::string  m_prefix;

    ANALYSIS::Primitive_Analysis             * p_analysis;
    std::vector<Observable_Data *>             m_obsdata;
    std::multimap<ATOOLS::Flavour,std::string> m_subsamples;

    int  m_qualifier;
    int  m_mode;

    void ReadInFinalSelectors(std::ifstream *,ANALYSIS::Final_Selector *&);
    void ReadInObservables(std::ifstream *);
    void SetUpObservables();
    // This is the place to add specific observables ....
    void SetUpSubSamples();
    void PrepareTerminate();
  public:
    Analysis_Handler(std::ifstream *,std::string=std::string("Hadrons"),const std::string& =std::string(""));
    ~Analysis_Handler();
    void DoAnalysis(ATOOLS::Blob_List * const, double);
    void Clear()        { p_analysis->ClearAllData(); }
    void Finish();

    void SetOutputPath(const std::string &);
    std::string Phase() { return m_phase; }
    bool On()           { return (p_analysis->NumberOfObservables()>0); }
  };   

  typedef std::map<std::string,Analysis_Handler *> AnalysesMap;
  typedef AnalysesMap::iterator                    AnalysesIter;
}

#endif
