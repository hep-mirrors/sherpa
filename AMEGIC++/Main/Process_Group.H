#ifndef AMEGIC_Main_Process_Group_H
#define AMEGIC_Main_Process_Group_H

#include "AMEGIC++/Main/Amegic_Base.H"
#include "PHASIC++/Process/Process_Group.H"

namespace MODEL { class Model_Base; }

namespace AMEGIC {

  class Topology;
  class Process_Base;
  class Amplitude_Handler;
  class Helicity;

  class Process_Group: public Amegic_Base, public PHASIC::Process_Group {
    
  protected :
    
    Topology * p_top;
    ATOOLS::Vec4D* p_testmoms;
    
    MODEL::Model_Base* p_model;
    std::string m_mfname;

    std::vector<AMEGIC::Process_Base *> m_umprocs, m_errprocs;
    std::list<std::string>              *p_channellibnames;
    std::set<std::string>               m_whitelist;

  public :

    // constructor
    Process_Group();

    // destructor
    ~Process_Group();

    // member functions
    bool FillIntegrator(PHASIC::Phase_Space_Handler *const psh);
    void EndOptimize();

    PHASIC::Process_Base *GetProcess(const PHASIC::Process_Info &pi) const;
    bool Initialize(PHASIC::Process_Base *const proc);

    AMEGIC::Process_Base *Partner() const;

    bool SetUpIntegrator();
    int  InitAmplitude(MODEL::Model_Base *model,Topology *top);

    Amplitude_Handler *GetAmplitudeHandler(); 
    Helicity *GetHelicity();
    bool NewLibs();
    std::string PSLibName();
    void SetPrintGraphs(std::string gpath);

    int PerformTests();

    void WriteMappingFile();
    void Minimize();

    void PrintProcessSummary(int=0); 

    // inline functions
    inline void SetMFname(const std::string &name) { m_mfname=name; }

    void FillAlphaHistogram(ATOOLS::Histogram*,double);
  };
}

#endif


