#ifndef SHERPA_PerturbativePhysics_Matrix_Element_Handler_H
#define SHERPA_PerturbativePhysics_Matrix_Element_Handler_H

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#include <map> 

namespace ATOOLS { 
  class Jet_Finder_Base;
  class Cluster_Definitions_Base;
  class Cluster_Amplitude;
#ifdef USING__GZIP
  class igzstream;
  class ogzstream;
#endif
}

namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF { class ISR_Handler; }

namespace PHASIC { 
  class Single_Process;
  class Selector_Key;
  struct nlo_type;
}

namespace SHERPA {

  typedef std::vector<PHASIC::Process_Base*> Process_Vector;

  class Shower_Handler;

  class Matrix_Element_Handler {
  public:


    typedef std::map<std::string,std::pair<int,double> >      MPDV_Map;
    typedef std::map<std::string,std::pair<int,std::string> > MPSV_Map;

  private :

    PHASIC::ME_Generators m_gens;
    Process_Vector     m_procs;

    PHASIC::Process_Base *p_proc;

    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;

    std::string m_path, m_file, m_processfile, m_selectorfile, m_respath;
    int         m_eventmode, m_seedmode;

    PHASIC::Weight_Info m_evtinfo;

    Shower_Handler *p_shower;

    double m_totalxs;
    size_t m_nlomode, m_ranidx;

#ifdef USING__GZIP
    ATOOLS::igzstream *p_ranin;
    ATOOLS::ogzstream *p_ranout;
#else
    std::ifstream *p_ranin;
    std::ofstream *p_ranout;
#endif

    void BuildDecays
    (PHASIC::Subprocess_Info &ACFS,const std::vector<std::string> &dectags);
    void BuildSingleProcessList
    (PHASIC::Process_Info &pi, const std::string &ini,
     const std::string &fin,const std::vector<std::string> &dectags,
     MPDV_Map &vmaxerr,MPDV_Map &vmaxeps,MPDV_Map &vefac,MPSV_Map &vefunc,
     MPSV_Map &vycut,MPSV_Map &vscale,MPSV_Map &vcoupl,MPSV_Map &vsfile,
     const std::string &gycut,const std::string &selfile);

    std::string MakeString(const std::vector<std::string> &in,
			   const size_t &first);

    void GetMaxCouplings(PHASIC::Process_Base *const proc,
			 size_t &oqcd,size_t &oew);

    size_t ExtractFlavours(PHASIC::Subprocess_Info &infos,std::string buffer);

    template <typename Type> Type ExtractMPvalue(const std::string& str);
    template <typename Type>
    void AddMPvalue(std::string lstr,std::string rstr,const Type &val,
		    std::map<std::string,std::pair<int,Type> >& dv,
		    const int nfs,const int &priority);
    template <typename Type>
    bool GetMPvalue(std::map<std::string,std::pair<int,Type> >& dv,
		    const int nfs,const std::string &pnid,Type &rv);
    template <typename Type>
    void ExtractMPvalues(std::string& str,std::map
			 <std::string,std::pair<int,Type> >& dv,
			 const int &priority);

    std::vector<PHASIC::Process_Base*> InitializeProcess(const PHASIC::Process_Info &pi);

    void BuildProcesses();

    void SetRandomSeed();

    bool GenerateUnweightedEvent();
    bool GenerateWeightedEvent(); 

  public :

    Matrix_Element_Handler(const std::string &path,const std::string &file,
                           const std::string &processfile,
                           const std::string &selectorfile);

    ~Matrix_Element_Handler();

    bool InitializeProcesses(MODEL::Model_Base *model,
			     BEAM::Beam_Spectra_Handler *beam,
			     PDF::ISR_Handler *isr);

    bool CalculateTotalXSecs();
    bool GenerateOneEvent();

    // inline functions
    inline PHASIC::Process_Base *Process() const { return p_proc; }
    inline Process_Vector        AllProcesses() const { return m_procs; }

    inline int EventGenerationMode() const { return m_eventmode; }

    inline PHASIC::Weight_Info WeightInfo() const { return m_evtinfo; }

    inline BEAM::Beam_Spectra_Handler *GetBeam() const { return p_beam; }
    inline PDF::ISR_Handler           *GetISR() const  { return p_isr;  }

    inline void SetShowerHandler(Shower_Handler *const sh) { p_shower=sh; }

    inline const std::string &RPath() const { return m_respath; }

    inline double TotalXS() const { return m_totalxs; }

  };// end of class Matrix_Element_Handler

  typedef std::map<std::string,Matrix_Element_Handler *> MEHandlersMap;
  typedef MEHandlersMap::iterator                        MEHandlerIter;

}// end of namespace SHERPA

#endif

