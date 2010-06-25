#ifndef PHASIC_Process_Process_Base_H
#define PHASIC_Process_Process_Base_H

#include "ATOOLS/Org/Smart_Pointer.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Selectors/Cut_Data.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Scales/KFactor_Setter_Base.H"

namespace ATOOLS { 
  class Cluster_Leg;
  class Cluster_Amplitude; 
  class Histogram;
}

namespace BEAM { class Beam_Spectra_Handler; }

namespace PDF { 
  class ISR_Handler;
  class Shower_Base;
}

namespace PHASIC {

  typedef std::map<int, int> FMMap;

  class Process_Integrator;
  class Phase_Space_Handler;
  class Combined_Selector;
  class ME_Generator_Base;
  class EEG_Multi_Channel;

  struct Weight_Info;

  class Process_Base;
  typedef std::map<std::string,Process_Base*> Process_Map;

  class ME_wgtinfo {
  public: 
    int m_nx;
    double m_w0;
    double* p_wx;
    double m_x1,m_x2,m_y1,m_y2;
    double m_renscale;
    ME_wgtinfo();
    ~ME_wgtinfo();
    ME_wgtinfo& operator*= (const double scal);
    void Flip();
    void AddMEweights(int n);
  };

  class Process_Base {
  protected:

    Process_Base       *p_parent, *p_selected, *p_mapproc;
    Process_Integrator *p_int;
    Combined_Selector  *p_selector;
    Cut_Data           *p_cuts;
    ME_Generator_Base  *p_gen;

    PDF::Shower_Base *p_shower;

    SP(Scale_Setter_Base)   p_scale;
    SP(KFactor_Setter_Base) p_kfactor;

    size_t      m_nin, m_nout, m_oqcd, m_oew;
    std::string m_name;

    ATOOLS::Flavour_Vector      m_flavs;
    Process_Info                m_pinfo;

    double m_last;
    bool   m_lookup, m_trigger;

    MODEL::Coupling_Map m_cpls;

    static void SortFlavours(Subprocess_Info &info,FMMap *const fmm=0);
    static void SortFlavours(Process_Info &pi);
    static void SortFlavours(std::vector<ATOOLS::Cluster_Leg*> &legs,
			     FMMap *const fmm=0);

  public:

    Process_Base();
    virtual ~Process_Base();

    // member functions
    Process_Base *Parent();
    Process_Base *Selected();

    virtual void SetScale(const Scale_Setter_Arguments &args) = 0;
    virtual void SetKFactor(const KFactor_Setter_Arguments &args) = 0;

    virtual void SetKFactorOn(const bool on);

    virtual size_t Size() const = 0;
    virtual Process_Base *operator[](const size_t &i) = 0;

    virtual void DeSelect() = 0;
    virtual bool SelectOne() = 0;

    virtual Weight_Info *OneEvent() = 0;
    virtual Weight_Info *WeightedEvent(const int mode=0) = 0;

    virtual double Differential(const ATOOLS::Vec4D_Vector &p) = 0;
    virtual double Differential(const ATOOLS::Cluster_Amplitude &ampl);
    virtual double Differential2() = 0;

    virtual bool GeneratePoint();
    virtual void AddPoint(const double &value);

    virtual void SetScaleSetter(const SP(Scale_Setter_Base) &scale);
    virtual void SetKFactorSetter(const SP(KFactor_Setter_Base) &kfactor);

    virtual bool IsGroup() const;

    virtual bool CalculateTotalXSec(const std::string &resultpath,
				    const bool create=false) = 0;
    virtual void SetLookUp(const bool lookup) = 0;

    virtual void Init(const Process_Info &pi,
		      BEAM::Beam_Spectra_Handler *const beamhandler,
		      PDF::ISR_Handler *const isrhandler);

    virtual bool FillIntegrator(Phase_Space_Handler *const psh);
    virtual void UpdateIntegrator(Phase_Space_Handler *const psh);

    virtual void SetGenerator(ME_Generator_Base *const gen);
    virtual void SetShower(PDF::Shower_Base *const ps);
    virtual void SetSelector(const Selector_Key &key);
    virtual void SetUpThreading();

    virtual void BuildCuts(Cut_Data *const cuts);
    virtual void UpdateCuts(const double &sp,const double &y,
			    Cut_Data *const cuts);

    virtual bool Trigger(const ATOOLS::Vec4D_Vector &p);
    virtual bool NoJetTrigger(const ATOOLS::Vec4D_Vector &p);
    virtual bool JetTrigger(const ATOOLS::Vec4D_Vector &p,
                    const ATOOLS::Flavour_Vector &fl,int n);
    virtual bool JetTrigger(const ATOOLS::Vec4D_Vector &p);

    virtual NLO_subevtlist *GetSubevtList();
    virtual ME_wgtinfo* GetMEwgtinfo() { return NULL; }

    void FillOnshellConditions();

    static void SortFlavours(ATOOLS::Cluster_Amplitude *const ampl);

    static std::string GenerateName(const ATOOLS::Cluster_Amplitude *ampl);

    static std::string GenerateName(const Subprocess_Info &info);
    static std::string GenerateName(const Subprocess_Info &ii,
				    const Subprocess_Info &fi);

    // inline functions
    inline void SetParent(Process_Base *const proc) { p_parent=proc; }

    inline void SetTrigger(const bool trigger) { m_trigger=trigger; }

    inline size_t NIn() const  { return m_nin;  }
    inline size_t NOut() const { return m_nout; }

    inline const ATOOLS::Flavour_Vector &Flavours() const { return m_flavs; }

    inline const std::string &Name() const { return m_name; }

    inline void SetOrderQCD(const size_t &oqcd) { m_oqcd=oqcd; }
    inline void SetOrderEW(const size_t &oew)   { m_oew=oew;   }

    inline size_t OrderQCD() const { return m_oqcd; }
    inline size_t OrderEW() const  { return m_oew;  }

    PDF::Shower_Base *Shower() const { return p_shower; }

    inline Process_Info       &Info()       { return m_pinfo; }
    inline const Process_Info &Info() const { return m_pinfo; }

    inline Process_Integrator *Integrator() const { return p_int; }

    inline double Last() const { return m_last; }

    inline bool LookUp() const { return m_lookup; }

    inline void SwapInOrder() 
    { std::swap<ATOOLS::Flavour>(m_flavs[0],m_flavs[1]); }

    inline const SP(Scale_Setter_Base) &ScaleSetter() const
    { return p_scale;   }
    inline const SP(KFactor_Setter_Base) &KFactorSetter() const
    { return p_kfactor; }

    inline bool Trigger() const  { return m_trigger; }
    inline bool IsMapped() const { return p_mapproc; }

    inline ME_Generator_Base *Generator() const { return p_gen;     }
    inline Process_Base      *MapProc() const   { return p_mapproc; }

    inline Cut_Data          * Cuts() const     { return p_cuts;     }
    inline Combined_Selector * Selector() const { return p_selector; }

    template <class PType> inline PType *Get() 
    { return dynamic_cast<PType*>(this); }
    template <class PType> inline const PType *Get() const 
    { return dynamic_cast<PType*>(this); }

  };// end of class Process_Base

}// end of namespace ATOOLS

#endif
