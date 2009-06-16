#ifndef PHASIC_Main_Phase_Space_Handler_H
#define PHASIC_Main_Phase_Space_Handler_H

#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <list>
#ifdef USING__Threading
#include <pthread.h>
#endif

namespace ATOOLS { 
  class Blob_Data_Base; 
  class Poincare;
}

namespace BEAM {
  class Beam_Spectra_Handler;
}

namespace PDF {
  class ISR_Handler;
}

namespace PHASIC {
  class Process_Integrator;
  class Multi_Channel;

  struct psm {
    
    enum code {
      normal_call = 0,
      no_lim_isr  = 32,
      no_dice_isr = 64
    };

  };// end of struct psm

  inline const psm::code 
  operator|(const psm::code &c1,const psm::code &c2)
  { return (psm::code)((int)c1|(int)c2); }
  inline const psm::code 
  operator&(const psm::code &c1,const psm::code &c2)
  { return (psm::code)((int)c1&(int)c2); }

  struct Weight_Info {
    
    double m_weight, m_dxs, m_ntrial, m_xf1, m_xf2;
    size_t m_mode;
    
    // constructor
    Weight_Info(const size_t mode=0,const double weight=0.0,
                const double dxs=0.0,const double ntrial=0.0,
		const double xf1=0.0,const double xf2=0.0):
      m_weight(weight), m_dxs(dxs), m_ntrial(ntrial),
      m_xf1(xf1), m_xf2(xf2), m_mode(mode) {}

  };// end of struct Weight_Info

  class Phase_Space_Integrator;
  class Process_Base;

  class Phase_Space_Handler: public ATOOLS::Tag_Replacer {
  private: 

    std::string m_name;

    Process_Integrator     *p_process, *p_active;
    Phase_Space_Integrator *p_integrator;  
    Cut_Data               *p_cuts;  

    ATOOLS::Algebra_Interpreter *p_eint;
 
    static ATOOLS::Integration_Info *p_info;

    ATOOLS::Info_Key  m_isrspkey, m_isrykey, m_beamspkey, m_beamykey;
    ATOOLS::Info_Key  m_mu2key[2];

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    Multi_Channel *p_fsrchannels, *p_isrchannels, *p_beamchannels;

    ATOOLS::Flavour_Vector p_flavours;
    ATOOLS::Vec4D_Vector   p_cms, p_lab;

    ATOOLS::Poincare *p_massboost;

    int  m_nin, m_nout, m_nvec, m_inttype, m_fin_opt, m_dmode;
    bool m_initialized, m_sintegrator;

    long int m_trials, m_maxtrials, m_sumtrials, m_events;

    double m_E, m_m[2], m_m2[2], m_s, m_smin, m_psweight;
    double m_result_1, m_result_2, m_flux, m_weight, m_error; 
    
    std::vector<std::vector<double> > m_stats;

    psm::code m_cmode;

#ifdef USING__Threading
    pthread_t m_met, m_pst;
    size_t m_sig, m_uset;
    pthread_mutex_t m_sme_mtx, m_sps_mtx, m_tme_mtx, m_tps_mtx;
    pthread_cond_t m_sme_cnd, m_sps_cnd, m_tme_cnd, m_tps_cnd;

    static void *CalculateME(void *arg);
    static void *CalculatePS(void *arg);
#endif

    bool MakeIncoming(ATOOLS::Vec4D *const p);
    bool Check4Momentum(const ATOOLS::Vec4D_Vector &p);

    void CalculateME();
    void CalculatePS();

  public:

    //constructor
    Phase_Space_Handler(Process_Integrator *,PDF::ISR_Handler *,
			BEAM::Beam_Spectra_Handler *, double error=-1.);

    //destructor
    ~Phase_Space_Handler();

    // member functions
    void WriteOut(const std::string &path,const bool force=false);
    bool ReadIn(const std::string &path,const size_t exclude=0);
    bool InitIncoming();
    bool CreateIntegrators();
    bool UpdateIntegrators();
    void InitCuts();
    
    double Integrate();
    double Differential();
    double Differential(Process_Integrator *const process,
			const psm::code mode=psm::normal_call);

    Weight_Info *OneEvent(Process_Base *const proc);
    Weight_Info *WeightedEvent(Process_Base *const proc,const int mode=0);

    void AddPoint(const double xs);
    void SetUpEnhance();

    double EnhanceFactor();

    static void TestPoint(ATOOLS::Vec4D *const p,
			  const size_t &nin,const size_t &nout,
			  const ATOOLS::Flavour_Vector &flavs);

    static ATOOLS::Integration_Info *const GetInfo();

    static void DeleteInfo();

    void AddStats(const std::vector<double> &stats);

    std::string   ReplaceTags(std::string &expr) const;
    ATOOLS::Term *ReplaceTags(ATOOLS::Term *term) const;

    // inline functions
    inline const ATOOLS::Vec4D_Vector &LabPoint() const { return p_lab; }
    inline const ATOOLS::Vec4D_Vector &CMSPoint() const { return p_cms; }

    inline const ATOOLS::Flavour_Vector &Flavs() const { return p_flavours; }

    inline Cut_Data *const Cuts() const { return p_cuts; }

    inline double   Weight() const    { return m_weight;    }
    inline double   Error() const     { return m_error;     }
    inline long int MaxTrials() const { return m_maxtrials; }

    inline long int IntType() const { return m_inttype; }

    inline void SetFSRIntegrator(Multi_Channel *const fsr)
    { p_fsrchannels=fsr;  }
 
    inline Multi_Channel *const BeamIntegrator() const  
    { return p_beamchannels; }
    inline Multi_Channel *const ISRIntegrator() const   
    { return p_isrchannels;  }
    inline Multi_Channel *const FSRIntegrator() const   
    { return p_fsrchannels;  }

    inline Process_Integrator *const Process() const { return p_process; }
    inline Process_Integrator *const Active() const  { return p_active;  }

    inline double PSWeight() const { return m_psweight; }

    inline void SetError(const double error) 
    { m_error=error; }
    inline void SetMaxTrials(const long int maxtrials)
    { m_maxtrials=maxtrials; }

    inline Phase_Space_Integrator *const Integrator() const 
    { return p_integrator; }

    inline const std::vector<std::vector<double> > Stats() const
    { return m_stats; }

  };// end of class Phase_Space_Handler

  /*!
    \class Phase_Space_Handler
    \brief the main steering class for integration and event generation
  */

  /*!
    \var ATOOLS::Integration_Info *Phase_Space_Handler::p_info

    Phase_Space_Handler distributes information on the various 
    integration variables and weights through an instance of Integration_Info.
    Each Single_Channel is able to gain access to the variables via the assignment 
    of an Info_Key during its initialization.
  */

  /*!
    \var ATOOLS::Info_Key Phase_Space_Handler::m_spkey
    
    initial key to allow Phase_Space_Handler access to s'
  */

  /*!
    \var ATOOLS::Info_Key Phase_Space_Handler::m_ykey
    
    initial key to allow Phase_Space_Handler access to y
  */

}//end of namespace PHASIC

#endif
