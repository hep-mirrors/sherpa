#ifndef PHASIC_Process_POWHEG_Process_H
#define PHASIC_Process_POWHEG_Process_H

#include "PHASIC++/Process/Process_Base.H"

namespace ATOOLS { class Cluster_Amplitude; }

namespace PDF {
  class Cluster_Definitions_Base;
  class NLOMC_Base;
}

namespace MODEL { class Running_AlphaS; }

namespace PHASIC {

  class ME_Generators;
  class POWHEG_Multi_Channel;
  class POWHEG_Process;

  struct ZH_Key {
    ATOOLS::NLO_subevt *p_sub;
    Process_Base *p_r;
    double m_res;
    inline ZH_Key(ATOOLS::NLO_subevt *const sub,
		  Process_Base *const r,const double &res):
      p_sub(sub), p_r(r), m_res(res) {}
  };

  typedef std::vector<ZH_Key> ZHKey_Vector;

  class POWHEG_Process : public Process_Base {
  private:

    ME_Generators &m_gens;

    Process_Base     *p_bproc, *p_viproc, *p_rproc, *p_sproc;
    PDF::NLOMC_Base  *p_powheg;

    POWHEG_Multi_Channel  *p_mc;
    MODEL::Running_AlphaS *p_as;

    double m_lastrs, m_lastb, m_lastvi, m_zhsum[2];
    int    m_smode, m_rmode;

    ATOOLS::StringIDipSet_Map m_iinfo;
    ATOOLS::StringDDipSet_Map m_dinfo;

    std::map<const size_t*,size_t> m_dmap;

    ZHKey_Vector m_zh[2];

    std::vector<double> m_rho;

    ATOOLS::Cluster_Amplitude *p_ampl;

    Process_Base *InitProcess(const Process_Info &pi,
			      const nlo_type::code nlotype,const bool real);

    Process_Base *FindProcess(const ATOOLS::Cluster_Amplitude *ampl,
			      const nlo_type::code type=nlo_type::lo,
			      const bool error=true) const;
    Process_Base *FindProcess(const ATOOLS::NLO_subevt *sub,
			      const nlo_type::code type=nlo_type::lo) const;

    ATOOLS::Cluster_Amplitude *CreateAmplitude
    (const ATOOLS::NLO_subevt *sub) const;

    double SelectZHProcess();
    double SelectBProcess();
    double SelectRProcess();
    double SelectProcess();

    double GetRho(const int mode);

  public:

    POWHEG_Process(ME_Generators &gens,NLOTypeStringProcessMap_Map *pmap);

    ~POWHEG_Process();

    void Init(const Process_Info &pi,
              BEAM::Beam_Spectra_Handler *const beam,
              PDF::ISR_Handler *const isr);
    bool InitSubtermInfo();

    bool IsGroup() const;
    size_t Size() const;
    Process_Base *operator[](const size_t &i);

    Weight_Info *OneEvent(const int wmode,const int mode=0);

    double Differential(const ATOOLS::Vec4D_Vector &p);
    double Differential2();

    bool Trigger(const ATOOLS::Vec4D_Vector &p);
    void InitCuts(Cut_Data *const cuts);
    void BuildCuts(Cut_Data *const cuts);

    bool FillIntegrator(Phase_Space_Handler *const psh);
    bool InitIntegrator(Phase_Space_Handler *const psh);
    void UpdateIntegrator(Phase_Space_Handler *const psh);

    bool CalculateTotalXSec(const std::string &resultpath,
                            const bool create=false);
    void SetLookUp(const bool lookup);
    void SetScale(const Scale_Setter_Arguments &scale);
    void SetKFactor(const KFactor_Setter_Arguments &args);
    void SetFixedScale(const std::vector<double> &s);
    void SetSelector(const Selector_Key &key);
    void SetShower(PDF::Shower_Base *const ps);

    void SetClusterDefinitions
    (PDF::Cluster_Definitions_Base *const cluster);

    void SetRBMap(ATOOLS::Cluster_Amplitude *ampl);

    double LocalKFactor(const ATOOLS::Vec4D_Vector &pb);

    ATOOLS::Cluster_Amplitude *GetAmplitude();

    inline double LastRS() const { return m_lastrs; }
    inline double LastB() const  { return m_lastb;  }
    inline double LastVI() const { return m_lastvi; }
    
    inline POWHEG_Multi_Channel *MultiChannel() { return p_mc; }

    inline const std::map<const size_t*,size_t> *DipoleMap() const
    { return &m_dmap; }

    inline size_t SMode() const { return m_smode; }

    inline void SetPOWHEG(PDF::NLOMC_Base *const powheg)
    { p_powheg=powheg; }

  };// end of class POWHEG_Process

}// end of namespace PHASIC

#endif
