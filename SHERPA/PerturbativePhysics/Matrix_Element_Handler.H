#ifndef SHERPA_PerturbativePhysics_Matrix_Element_Handler_H
#define SHERPA_PerturbativePhysics_Matrix_Element_Handler_H

#include "MODEL/Main/Model_Base.H"
#include "BEAM/Main/Beam_Spectra_Handler.H"
#include "PDF/Main/ISR_Handler.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
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
  struct nlo_type;
}

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
	m_vsfile, m_vnloqcdmode, m_vnloqcdpart, m_vnloewmode, m_vnloewpart,
	m_vnlosubv, m_vmegen, m_vrsmegen, m_vloopgen, m_vint, m_vrsint,
	m_vgpath, m_vaddname;
      MPIV_Map m_vnmaxq, m_vnminq, m_vamegicmhv,
	m_vntchan, m_vmtchan, m_voqcd, m_voew, m_vitmin, m_vrsitmin;

      std::string m_gycut;
      int m_cutcore;
      inline Processblock_Info(): m_cutcore(0) {}
    };// end of struct Processblock_Info

  private :

    PHASIC::ME_Generators  m_gens;
    PHASIC::Process_Vector m_procs;

    PHASIC::Process_Base *p_proc;

    BEAM::Beam_Spectra_Handler *p_beam;
    PDF::ISR_Handler           *p_isr;
    MODEL::Model_Base          *p_model;

    std::string m_path, m_file, m_processfile, m_selectorfile, m_respath;
    int         m_eventmode, m_seedmode, m_hasnlo, m_rsadd;

    ATOOLS::Weight_Info m_evtinfo;

    Shower_Handler   *p_shower;
    PDF::NLOMC_Base  *p_nlomc;

    double m_sum, m_ovwth, m_weightfactor;
    size_t m_globalnlomode, m_ranidx, m_fosettings;

#ifdef USING__GZIP
    ATOOLS::igzstream *p_ranin;
    ATOOLS::ogzstream *p_ranout;
#else
    std::ifstream *p_ranin;
    std::ofstream *p_ranout;
#endif

    ProcMap_Vector m_pmaps;

    void BuildDecays
    (PHASIC::Subprocess_Info &ACFS,const std::vector<std::string> &dectags);
    void BuildSingleProcessList
    (PHASIC::Process_Info &pi,Processblock_Info &pbi,
     const std::string &ini,const std::string &fin,
     const std::vector<std::string> &dectags);

    std::string MakeString(const std::vector<std::string> &in,
			   const size_t &first);

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

    std::vector<PHASIC::Process_Base*> InitializeProcess
    (const PHASIC::Process_Info &pi,PHASIC::NLOTypeStringProcessMap_Map *&pmap);

    void BuildProcesses();

    void SetRandomSeed();
    void InitNLOMC();

  public :

    Matrix_Element_Handler(const std::string &path,const std::string &file,
                           const std::string &processfile,
                           const std::string &selectorfile);

    ~Matrix_Element_Handler();

    int InitializeProcesses(MODEL::Model_Base *model,
			     BEAM::Beam_Spectra_Handler *beam,
			     PDF::ISR_Handler *isr);

    bool CalculateTotalXSecs();
    bool GenerateOneEvent();

    // inline functions
    inline PHASIC::Process_Base *Process() const { return p_proc; }

    inline PHASIC::Process_Vector AllProcesses() const { return m_procs; }
    
    inline ProcMap_Vector ProcMaps() const { return m_pmaps;}

    inline int EventGenerationMode() const { return m_eventmode; }

    inline ATOOLS::Weight_Info WeightInfo() const { return m_evtinfo; }

    inline BEAM::Beam_Spectra_Handler *GetBeam() const { return p_beam; }
    inline PDF::ISR_Handler           *GetISR() const  { return p_isr;  }

    inline void SetShowerHandler(Shower_Handler *const sh) { p_shower=sh; }

    inline Shower_Handler *Shower() const { return p_shower; }

    inline double WeightFactor() const { return m_weightfactor; }

    inline double Sum() const { return m_sum; }

    inline const std::string &Path() const { return m_path; }
    inline const std::string &File() const { return m_file; }
    
    inline int HasNLO() const  { return m_hasnlo;  }

    double GetWeight(const ATOOLS::Cluster_Amplitude &ampl,
		     const PHASIC::nlo_type::code type,
		     const int mode) const;

  };// end of class Matrix_Element_Handler

}// end of namespace SHERPA

#endif

