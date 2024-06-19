#ifndef SHERPA_PerturbativePhysics_Matrix_Element_Handler_H
#define SHERPA_PerturbativePhysics_Matrix_Element_Handler_H

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "YFS/Main/YFS_Handler.H"
#include "REMNANTS/Main/Remnant_Handler.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Phys/Weight_Info.H"

#include <map>

namespace ATOOLS {
#ifdef USING__GZIP
  class igzstream;
  class ogzstream;
#endif
}

namespace BEAM { class Beam_Spectra_Handler; }
namespace PDF {
  class ISR_Handler;
  class NLOMC_Base;
}

namespace PHASIC {
  class Single_Process;
  class Selector_Key;
}

namespace YFS { class YFS_Handler; }

namespace SHERPA {

  class Shower_Handler;

  class Matrix_Element_Handler {
  public:

    typedef std::vector<PHASIC::NLOTypeStringProcessMap_Map*> ProcMap_Vector;

    typedef std::map<std::string,std::pair<int,int> >         MPIV_Map;
    typedef std::map<std::string,std::pair<int,double> >      MPDV_Map;
    typedef std::map<std::string,std::pair<int,std::string> > MPSV_Map;

    struct Processblock_Info {
      MPDV_Map m_vmaxerr, m_vmaxeps, m_vefac, m_vrsefac;
      MPSV_Map m_veobs, m_vefunc, m_vycut, m_vscale, m_vcoupl, m_vkfac,
        m_vnlomode, m_vnlopart, m_vnlocpl,
	m_vnlosubv, m_vasscontribs, m_vmegen, m_vrsmegen, m_vloopgen, m_vint,
        m_vrsint, m_vgpath, m_vcpl, m_vmaxcpl, m_vmincpl, m_vacpl,
	m_vmaxacpl, m_vminacpl, m_vaddname, m_vspecial;
      MPIV_Map m_vnmaxq, m_vnminq, m_vamegicmhv,
	m_vntchan, m_vmtchan, m_vitmin, m_vrsitmin;
      ATOOLS::Scoped_Settings m_selectors;

      std::string m_gycut;
      int m_cutcore;
      size_t m_sort;
      inline Processblock_Info(): m_cutcore(0), m_sort(3) {}
    };// end of struct Processblock_Info

  private :

    PHASIC::ME_Generators  m_gens;
    PHASIC::Process_Vector m_procs;

    PHASIC::Process_Base *p_proc;

    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;
    YFS::YFS_Handler           *p_yfs;
    REMNANTS::Remnant_Handler  *p_remnants;
    MODEL::Model_Base          *p_model;

    std::string m_respath;
    int         m_eventmode, m_seedmode, m_hasnlo;
    bool        m_pilotrunenabled;
    int         m_nloadd, m_ewaddmode, m_qcdaddmode;

    ATOOLS::Weight_Info m_evtinfo;

    Shower_Handler    *p_shower;
    PDF::NLOMC_Base   *p_nlomc;

    double m_sum, m_ovwth, m_weightfactor;
    size_t m_ranidx, m_fosettings;
    ATOOLS::nlo_mode::code m_nlomode;

#ifdef USING__GZIP
    ATOOLS::igzstream *p_ranin;
    ATOOLS::ogzstream *p_ranout;
#else
    std::ifstream *p_ranin;
    std::ofstream *p_ranout;
#endif

    ProcMap_Vector m_pmaps;

    void RegisterDefaults();
    void RegisterMainProcessDefaults(ATOOLS::Scoped_Settings&);

    void LimitCouplings
    (MPSV_Map &pbi,const size_t &nfs,const std::string &pnid,
     std::vector<double> &mincpl,std::vector<double> &maxcpl,const int mode);
    void BuildDecays
    (PHASIC::Subprocess_Info &ACFS,const std::vector<std::string> &dectags);

    struct Single_Process_List_Args {
      PHASIC::Process_Info pi;
      Processblock_Info pbi;
      std::string ini;
      std::string fin;
      std::vector<std::string> dectags;
    };
    void BuildSingleProcessList(Single_Process_List_Args&);

    std::string MakeOrderString(ATOOLS::Scoped_Settings&&) const;
    std::string MakeString(const std::vector<std::string>&) const;

    template <typename Type> Type ExtractMPvalue(const std::string& str);
    template <typename Type>
    void AddMPvalue(std::string lstr,std::string rstr,const Type &val,
		    std::map<std::string,std::pair<int,Type> >& dv,
		    const int nfs,const int &priority);
    template <typename Type>
    bool GetMPvalue(std::map<std::string,std::pair<int,Type> >& dv,
		    const int nfs,const std::string &pnid,Type &rv);
    template <typename Type>
    void ExtractMPvalues(std::string str,
                         std::pair<size_t, size_t> multirange,
                         const int &priority,
                         std::map<std::string,std::pair<int,Type> >&) const;

    std::vector<PHASIC::Process_Base*> InitializeProcess(
        PHASIC::Process_Info, PHASIC::NLOTypeStringProcessMap_Map*&);
    std::vector<PHASIC::Process_Base*> InitializeSingleProcess
    (const PHASIC::Process_Info &pi,PHASIC::NLOTypeStringProcessMap_Map *&pmap);
    void CheckInitialStateOrdering(const PHASIC::Process_Info&);

    void BuildProcesses();
    void ReadFinalStateMultiIndependentProcessSettings(
      const std::string& procname, ATOOLS::Scoped_Settings, Single_Process_List_Args&);
    void ReadFinalStateMultiSpecificProcessSettings(
      ATOOLS::Scoped_Settings,
      Single_Process_List_Args&,
      std::string range="-") const;

    void InitNLOMC();

    // event generation helper functions
    void SetRandomSeed();
    bool GenerateOneTrialEvent();

  public :

    Matrix_Element_Handler(MODEL::Model_Base*);

    ~Matrix_Element_Handler();

    int InitializeProcesses(BEAM::Beam_Spectra_Handler*,
                            PDF::ISR_Handler*,
                            YFS::YFS_Handler*);

    int InitializeTheReweighting(ATOOLS::Variations_Mode mode);


    bool CalculateTotalXSecs();
    bool GenerateOneEvent();

    // inline functions
    inline PHASIC::Process_Base *Process() const { return p_proc; }

    inline PHASIC::Process_Vector AllProcesses() const { return m_procs; }

    inline PHASIC::ME_Generators *Generators() { return &m_gens; }

    inline ProcMap_Vector ProcMaps() const { return m_pmaps;}

    inline int EventGenerationMode() const { return m_eventmode; }

    inline ATOOLS::Weight_Info WeightInfo() const { return m_evtinfo; }

    inline void SetRemnantHandler(REMNANTS::Remnant_Handler *const rh) {
      p_remnants = rh;
    }
    inline void SetShowerHandler(Shower_Handler *const sh) {
      p_shower=sh;
    }
    inline BEAM::Beam_Spectra_Handler *GetBeam() const  { return p_beam; }
    inline PDF::ISR_Handler           *GetISR() const   { return p_isr;  }
    inline YFS::YFS_Handler           *GetYFS() const   { return p_yfs;  }
    inline REMNANTS::Remnant_Handler  *Remnants() const { return p_remnants; }
    inline Shower_Handler *Shower() const { return p_shower; }

    inline double WeightFactor() const { return m_weightfactor; }

    inline double Sum() const { return m_sum; }

    inline int HasNLO() const  { return m_hasnlo;  }

    double GetWeight(const ATOOLS::Cluster_Amplitude &ampl,
                     const ATOOLS::nlo_type::code type,
		     const int mode) const;

    inline int SeedMode() const { return m_seedmode; }

    static size_t ExtractFlavours(PHASIC::Subprocess_Info &infos,std::string buffer);

    void CheckAssociatedContributionsSetup(
        const std::set<ATOOLS::asscontrib::type>& calculated_asscontribs) const;
  };// end of class Matrix_Element_Handler

}// end of namespace SHERPA

#endif

