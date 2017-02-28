#ifndef PHASIC_Process_Process_Base_H
#define PHASIC_Process_Process_Base_H

#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Selectors/Cut_Data.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"
#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "ATOOLS/Org/CXXFLAGS.H"

namespace ATOOLS { 
  class Cluster_Leg;
  class Cluster_Amplitude; 
  class Histogram;
  class Variation_Weights;
  struct Decay_Info;
  typedef std::vector<Decay_Info* > DecayInfo_Vector;
}

namespace BEAM { class Beam_Spectra_Handler; }

namespace PDF { 
  class ISR_Handler;
  class Shower_Base;
}

namespace METOOLS { class Spin_Amplitudes; }

namespace PHASIC {

  typedef std::map<int,int> FMMap;

  class Process_Integrator;
  class Phase_Space_Handler;
  class Phase_Space_Integrator;
  class Combined_Selector;
  class ME_Generator_Base;

  class Process_Base;
  class BBar_Multi_Channel;

  typedef std::vector<PHASIC::Process_Base*> Process_Vector;

  typedef std::map<std::string,Process_Base*> StringProcess_Map;

  typedef std::map<nlo_type::code,StringProcess_Map*>
  NLOTypeStringProcessMap_Map;

#ifdef USING__Threading
  class Process_Base;

  struct MEH_TID {
    pthread_t m_id;
    Process_Base *p_proc;
    std::string m_respath;
    MEH_TID *p_prev, *p_next;
    size_t m_i, m_r;
    pthread_mutex_t m_s_mtx, *p_t_mtx, m_r_mtx;
    pthread_cond_t m_s_cnd, *p_t_cnd;
    MEH_TID(Process_Base *const proc,const size_t i):
      p_proc(proc), p_prev(this), p_next(this),
      m_i(i), m_r(0) {}
  };// end of struct MEH_TID

  typedef std::vector<MEH_TID*> MEH_TID_Vector; 
#endif

  class Process_Base {
  protected:

    Process_Base       *p_parent, *p_selected, *p_mapproc, *p_sproc, *p_proc;
    Process_Integrator *p_int;
    Combined_Selector  *p_selector;
    Cut_Data           *p_cuts;
    ME_Generator_Base  *p_gen;

    PDF::Shower_Base   *p_shower;
    ATOOLS::Variation_Weights *p_variationweights;
    bool m_variationweightsowned;
    BBar_Multi_Channel *p_mc;


    Scale_Setter_Base   *p_scale;
    KFactor_Setter_Base *p_kfactor;

    size_t      m_nin, m_nout, m_tinfo, m_mcmode, m_cmode;
    std::string m_name;

    std::vector<double> m_maxcpl, m_mincpl;

    ATOOLS::Flavour_Vector m_flavs;
    Process_Info           m_pinfo;
    ATOOLS::DecayInfo_Vector m_decins;

    ATOOLS::ME_Weight_Info m_mewgtinfo;

    double m_last, m_lastb, m_symfac, m_issymfac;
    bool   m_lookup, m_use_biweight;

    MODEL::Coupling_Map m_cpls;

    NLOTypeStringProcessMap_Map *p_apmap;

#ifdef USING__Threading
    MEH_TID *p_tid;
#endif

    static int s_usefmm;

    static void SortFlavours(Subprocess_Info &info,FMMap *const fmm=0);
    static void SortFlavours(std::vector<ATOOLS::Cluster_Leg*> &legs,
			     FMMap *const fmm=0);

  public:

    Process_Base();
    virtual ~Process_Base();

    // member functions
    Process_Base *Parent();
    Process_Base *Selected();
    size_t SelectedIndex();

    bool SetSelected(Process_Base *const proc);
    size_t SynchronizeSelectedIndex(Process_Base &);

    double LastPlus();
    double LastMinus();

    void SetSProc(Process_Base *sproc);
    void SetBBarMC(BBar_Multi_Channel *mc);

    int NaiveMapping(Process_Base *proc) const;

    std::string ShellName(std::string name="") const;

    void MPISync(Phase_Space_Integrator *const psi);
    bool MPISyncRequest();

#ifdef USING__Threading
    virtual void AddMEHThread(MEH_TID_Vector &cts,void *(*CalcFunc)(void*));
#endif

    virtual void FillProcessMap(NLOTypeStringProcessMap_Map *apmap);

    virtual void SetScale(const Scale_Setter_Arguments &args) = 0;
    virtual void SetKFactor(const KFactor_Setter_Arguments &args) = 0;

    virtual void SetFixedScale(const std::vector<double> &s);
    virtual void SetSelectorOn(const bool on);
    virtual void SetUseBIWeight(bool on);

    virtual size_t SetMCMode(const size_t mcmode);
    virtual size_t SetClusterMode(const size_t mcmode);

    virtual size_t Size() const = 0;
    virtual Process_Base *operator[](const size_t &i) = 0;

    virtual ATOOLS::Weight_Info *OneEvent(const int wmode,const int mode=0) = 0;

    virtual double Differential(const ATOOLS::Vec4D_Vector &p) = 0;
    virtual double Differential(const ATOOLS::Cluster_Amplitude &ampl,
                                int mode=0);

    virtual bool GeneratePoint();
    virtual void AddPoint(const double &value);
    virtual bool ReadIn(const std::string &pid);
    virtual void WriteOut(const std::string &pid);
    virtual void EndOptimize();
    virtual void MPICollect(std::vector<double> &sv,size_t &i);
    virtual void MPIReturn(std::vector<double> &sv,size_t &i);
    virtual void MPISync(const int mode=0);

    virtual bool IsGroup() const;

    virtual bool CalculateTotalXSec(const std::string &resultpath,
				    const bool create=false) = 0;
    virtual void SetLookUp(const bool lookup) = 0;

    virtual void Init(const Process_Info &pi,
		      BEAM::Beam_Spectra_Handler *const beamhandler,
		      PDF::ISR_Handler *const isrhandler,const int mode=0);

    virtual bool FillIntegrator(Phase_Space_Handler *const psh);
    virtual bool InitIntegrator(Phase_Space_Handler *const psh);
    virtual void UpdateIntegrator(Phase_Space_Handler *const psh);

    virtual void SetGenerator(ME_Generator_Base *const gen);
    virtual void SetShower(PDF::Shower_Base *const ps);
    virtual void SetVariationWeights(ATOOLS::Variation_Weights *const);
    void SetOwnedVariationWeights(ATOOLS::Variation_Weights *);
    ATOOLS::Variation_Weights *VariationWeights() const { return p_variationweights; }
    virtual void SetSelector(const Selector_Key &key);

    virtual void InitCuts(Cut_Data *const cuts);
    virtual void BuildCuts(Cut_Data *const cuts);

    virtual bool Trigger(const ATOOLS::Vec4D_Vector &p);

    virtual ATOOLS::NLO_subevtlist *GetSubevtList();
    virtual ATOOLS::NLO_subevtlist *GetRSSubevtList();

    virtual void SetRBMap(ATOOLS::Cluster_Amplitude *ampl);
    virtual void InitPSHandler(const double &maxerror,
			       const std::string eobs,
			       const std::string efunc);

    virtual void FillOnshellConditions();

    static void SortFlavours(ATOOLS::Cluster_Amplitude *const ampl,
			     const int mode=1);
    static void SortFlavours(Process_Info &pi,const int mode=1);

    static std::string GenerateName(const ATOOLS::Cluster_Amplitude *ampl);
    static std::string GenerateName(const ATOOLS::NLO_subevt *sub,
				    const size_t &nin);

    static std::string GenerateName(const Subprocess_Info &info);
    static std::string GenerateName(const Subprocess_Info &ii,
				    const Subprocess_Info &fi);

    // inline functions
    inline void SetParent(Process_Base *const proc) { p_parent=proc; }

    inline size_t NIn() const  { return m_nin;  }
    inline size_t NOut() const { return m_nout; }

    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_flavs; }

    inline const std::string &Name() const { return m_name; }

    inline void SetMaxOrders(const std::vector<double> &o) { m_maxcpl=o; }
    inline void SetMinOrders(const std::vector<double> &o) { m_mincpl=o; }

    inline void SetMaxOrder(const size_t &id,const double o)
    { if (m_maxcpl.size()<=id) m_maxcpl.resize(id+1,0); m_maxcpl[id]=o; }
    inline void SetMinOrder(const size_t &id,const double o)
    { if (m_mincpl.size()<=id) m_mincpl.resize(id+1,0); m_mincpl[id]=o; }

    inline const std::vector<double> &MaxOrders() const { return m_maxcpl; }
    inline const std::vector<double> &MinOrders() const { return m_mincpl; }

    inline double MaxOrder(const size_t &id) const { return m_maxcpl[id]; }
    inline double MinOrder(const size_t &id) const { return m_mincpl[id]; }

    PDF::Shower_Base *Shower() const { return p_shower; }

    inline Process_Info       &Info()       { return m_pinfo; }
    inline const Process_Info &Info() const { return m_pinfo; }

    inline ATOOLS::DecayInfo_Vector       &DecayInfos()       {return m_decins;}
    inline const ATOOLS::DecayInfo_Vector &DecayInfos() const {return m_decins;}

    inline ATOOLS::ME_Weight_Info *GetMEwgtinfo() { return &m_mewgtinfo; }
    inline void SetMEwgtinfo(const ATOOLS::ME_Weight_Info i) { m_mewgtinfo=i; }


    inline Process_Integrator *Integrator() const { return p_int; }

    virtual void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                                std::vector<std::vector<Complex> >& cols);

    inline double Last() const { return m_last; }
    inline double LastB() const { return m_lastb; }

    inline bool LookUp() const { return m_lookup; }

    inline double SymFac() const { return m_symfac; }

    inline Scale_Setter_Base *ScaleSetter(const int map=0) const
    { return map?(p_mapproc?p_mapproc->p_scale:p_scale):p_scale; }
    inline KFactor_Setter_Base *KFactorSetter(const int map=0) const
    { return map?(p_mapproc?p_mapproc->p_kfactor:p_kfactor):p_kfactor; }

    inline bool IsMapped() const { return p_mapproc; }

    inline size_t ThreadInfo() const { return m_tinfo; }

    inline size_t MCMode() const      { return m_mcmode; }
    inline size_t ClusterMode() const { return m_cmode;  }

    inline ME_Generator_Base *Generator() const { return p_gen; }

    inline Process_Base *MapProc() const  { return p_mapproc; }
    inline Process_Base *Proc() const     { return p_proc;    }

    inline Cut_Data          * Cuts() const { return p_cuts; }
    inline Combined_Selector * Selector() const 
    { return p_mapproc?p_mapproc->p_selector:p_selector; }

    template <class PType> inline PType *Get() 
    { return dynamic_cast<PType*>(this); }
    template <class PType> inline const PType *Get() const 
    { return dynamic_cast<PType*>(this); }

    inline NLOTypeStringProcessMap_Map *AllProcs() { return p_apmap; }

    inline BBar_Multi_Channel *BBarMC() const { return p_mc; }

#ifdef USING__Threading
    inline void SetTID(MEH_TID *const tid) { p_tid=tid; }
#endif

  };// end of class Process_Base

}// end of namespace ATOOLS

#endif
