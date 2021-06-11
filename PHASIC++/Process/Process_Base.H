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
  struct Decay_Info;
  typedef std::vector<Decay_Info* > DecayInfo_Vector;
  struct Weight_Info;
}

namespace BEAM { class Beam_Spectra_Handler; }

namespace PDF { 
  class ISR_Handler;
  class Shower_Base;
  class NLOMC_Base;
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

  typedef std::map<ATOOLS::nlo_type::code,StringProcess_Map*>
  NLOTypeStringProcessMap_Map;

  class Process_Base {
  protected:

    Process_Base       *p_parent, *p_selected, *p_mapproc, *p_sproc, *p_caller;
    Process_Integrator *p_int;
    Combined_Selector  *p_selector;
    Cut_Data           *p_cuts;
    ME_Generator_Base  *p_gen;

    PDF::Shower_Base   *p_shower;
    PDF::NLOMC_Base    *p_nlomc;

    BBar_Multi_Channel *p_mc;

    Scale_Setter_Base   *p_scale;
    KFactor_Setter_Base *p_kfactor;

    size_t      m_nin, m_nout, m_mcmode, m_cmode;
    std::string m_name, m_resname;

    std::vector<double> m_maxcpl, m_mincpl;

    ATOOLS::Flavour_Vector m_flavs;
    Process_Info           m_pinfo;
    ATOOLS::DecayInfo_Vector m_decins;

    ATOOLS::ME_Weight_Info m_mewgtinfo;

    ATOOLS::Weights_Map m_last, m_lastb;
    double m_symfac, m_issymfac;
    bool   m_lookup, m_use_biweight;

    bool   m_hasinternalscale;
    double m_internalscale;
    MODEL::Coupling_Map m_cpls;

    NLOTypeStringProcessMap_Map *p_apmap;

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

    /*!
     * Calculate the nominal differential
     *
     * Use the Variations_Mode argument to steer if only the nominal, or the
     * nominal plus all on-the-fly variations are calculated and returned.
     *
     * The nominal differential can be retrieved using Last, the Born part of
     * it via LastB, and the full event weights (i.e. the return value) via
     * EventWeights.  The ME weight information can be retrieved using
     * GetMEwgtinfo.
     *
     * Subclasses implementing this virtual function must make sure to set
     * update the corresponding member variables.
     */
    virtual ATOOLS::Weights_Map Differential(const ATOOLS::Vec4D_Vector&,
                                             ATOOLS::Variations_Mode) = 0;

    virtual ATOOLS::Weights_Map Differential(const ATOOLS::Cluster_Amplitude&,
                                             ATOOLS::Variations_Mode varmode=ATOOLS::Variations_Mode::nominal_only,
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

    virtual int  PerformTests();
    virtual bool InitScale();
    virtual void Init(const Process_Info &pi,
		      BEAM::Beam_Spectra_Handler *const beamhandler,
		      PDF::ISR_Handler *const isrhandler,const int mode=0);

    virtual bool FillIntegrator(Phase_Space_Handler *const psh);
    virtual bool InitIntegrator(Phase_Space_Handler *const psh);
    virtual void UpdateIntegrator(Phase_Space_Handler *const psh);

    virtual void SetGenerator(ME_Generator_Base *const gen);
    virtual void SetShower(PDF::Shower_Base *const ps);
    virtual void SetNLOMC(PDF::NLOMC_Base *const mc);
    const ATOOLS::Weights_Map& WeightsMap() const { return m_last; }
    virtual void SetSelector(const Selector_Key &key);
    virtual void SetCaller(Process_Base *const proc);

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

    static std::string BaseName(const std::string &name,
				const std::string &addname=std::string(""));

    static std::string GenerateName(const ATOOLS::Cluster_Amplitude *ampl);
    static std::string GenerateName(const ATOOLS::NLO_subevt *sub,
				    const size_t &nin);

    static std::string GenerateName(const Subprocess_Info &info);
    static std::string GenerateName(const Subprocess_Info &ii,
				    const Subprocess_Info &fi);

    // inline functions
    inline void SetParent(Process_Base *const proc) { p_parent=proc; }

    inline size_t NIn()                             const { return m_nin;  }
    inline size_t NOut()                            const { return m_nout; }
    const ATOOLS::Vec4D_Vector          &Momenta()  const;
    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_flavs; }
    inline void SetFlavours(const ATOOLS::Flavour_Vector& flavs) { m_flavs = flavs; }
    inline virtual std::vector<std::vector<int> > * Colours() const { return NULL; }

    inline const std::string &Name() const        { return m_name;    }
    inline const std::string &ResultsName() const { return m_resname; }

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
    PDF::NLOMC_Base  *NLOMC() const  { return p_nlomc;  }

    inline Process_Info       &Info()       { return m_pinfo; }
    inline const Process_Info &Info() const { return m_pinfo; }

    inline ATOOLS::DecayInfo_Vector       &DecayInfos()       {return m_decins;}
    inline const ATOOLS::DecayInfo_Vector &DecayInfos() const {return m_decins;}

    inline ATOOLS::ME_Weight_Info *GetMEwgtinfo() { return &m_mewgtinfo; }
    inline void SetMEwgtinfo(const ATOOLS::ME_Weight_Info i) { m_mewgtinfo=i; }

    inline Process_Integrator *Integrator() const { return p_int; }
    
    inline virtual const bool   HasInternalScale() const { return m_hasinternalscale; }
    inline virtual const double InternalScale()    const { return m_internalscale; }

    virtual void FillAmplitudes(std::vector<METOOLS::Spin_Amplitudes>& amps,
                                std::vector<std::vector<Complex> >& cols);

    inline double Last() const { return m_last.Nominal(); }
    inline double LastB() const { return m_lastb.Nominal(); }

    inline bool LookUp() const { return m_lookup; }

    inline double SymFac() const { return m_symfac; }
    inline double ISSymFac() const { return m_issymfac; }

    inline Scale_Setter_Base *ScaleSetter(const int map=0) const
    { return map?(p_mapproc?p_mapproc->p_scale:p_scale):p_scale; }
    inline KFactor_Setter_Base *KFactorSetter(const int map=0) const
    { return map?(p_mapproc?p_mapproc->p_kfactor:p_kfactor):p_kfactor; }

    inline bool IsMapped() const { return p_mapproc; }

    inline size_t MCMode() const      { return m_mcmode; }
    inline size_t ClusterMode() const { return m_cmode;  }

    inline ME_Generator_Base *Generator() const { return p_gen; }

    inline Process_Base *MapProc() const  { return p_mapproc; }
    inline Process_Base *Caller() const   { return p_caller;  }

    inline Cut_Data          * Cuts() const { return p_cuts; }
    inline Combined_Selector * Selector() const 
    { return p_mapproc?p_mapproc->p_selector:p_selector; }
    inline virtual const double SPrimeMin() const { return -1.; }
    inline virtual const double SPrimeMax() const { return -1.; }
    template <class PType> inline PType *Get() 
    { return dynamic_cast<PType*>(this); }
    template <class PType> inline const PType *Get() const 
    { return dynamic_cast<PType*>(this); }

    inline NLOTypeStringProcessMap_Map *AllProcs() { return p_apmap; }

    inline BBar_Multi_Channel *BBarMC() const { return p_mc; }
    virtual bool FillFinalState(const ATOOLS::Vec4D_Vector &p) { return true; }
  };// end of class Process_Base

  class Order_Flavour {
    FMMap* p_fmm;
    int Order_SVFT(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b) 
    {
      if (a.IsScalar() && !b.IsScalar()) return 1;
      if (a.IsVector() && !b.IsScalar() && 
	  !b.IsVector()) return 1;
      if (a.IsFermion() && !b.IsFermion() && 
	  !b.IsScalar() && !b.IsVector()) return 1;
      return 0;
    }
    int Order_Multi(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
    {
      if ((*p_fmm)[int(a.Kfcode())]==0 || 
	  (*p_fmm)[int(b.Kfcode())]==0) return 0;
      if ((*p_fmm)[int(a.Kfcode())]>
	  (*p_fmm)[int(b.Kfcode())]) return 1;
      return 0;
    }
    int Order_Photons(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
    {
      if (a.Strong() && a.Mass() && b.IsPhoton()) return 1;
      return 0;
    }

  public:
    int operator()(const ATOOLS::Flavour &a,const ATOOLS::Flavour &b)
    {
      if (a.Priority()>b.Priority()) return 1;
      if (a.Priority()<b.Priority()) return 0;
      if (Order_Photons(a,b)) return 1;
      if (Order_Photons(b,a)) return 0;
      if (!a.Strong()&&b.Strong()) return 1;
      if (a.Strong()&&!b.Strong()) return 0;
      if (a.Mass()>b.Mass()) return 1;
      if (a.Mass()<b.Mass()) return 0;
      if (p_fmm) {
	if (Order_Multi(a,b)) return 1;
	if (Order_Multi(b,a)) return 0;
      }
      if (Order_SVFT(a,b)) return 1;
      if (Order_SVFT(b,a)) return 0;
      if (!a.IsAnti()&&b.IsAnti()) return 1;
      if (a.IsAnti()&&!b.IsAnti()) return 0;
      return a.Kfcode()<b.Kfcode();
    }
    Order_Flavour(FMMap* fmm): p_fmm(fmm) {}
    int operator()(const Subprocess_Info &a,const Subprocess_Info &b)
    { return (*this)(a.m_fl,b.m_fl); }
    int operator()(const ATOOLS::Cluster_Leg *a,
		   const ATOOLS::Cluster_Leg *b)
    { return (*this)(a->Flav(),b->Flav()); }
  };// end of class Order_Flavour

}// end of namespace ATOOLS

#endif
