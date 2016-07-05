#ifndef PHASIC_Main_Phase_Space_Handler_H
#define PHASIC_Main_Phase_Space_Handler_H

#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Org/Info_Key.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include "ATOOLS/Math/Histogram.H"
#include <list>

namespace ATOOLS { 
  class Blob_Data_Base; 
  class Poincare;
  class Mass_Selector;
  struct Weight_Info;
}

namespace BEAM {
  class Beam_Spectra_Handler;
}

namespace PDF {
  class ISR_Handler;
}

namespace SHERPA {
  class Variation_Weights;
}

namespace PHASIC {
  class Process_Integrator;
  class Multi_Channel;
  class Subprocess_Info;
  class Process_Info;
  class Enhance_Observable_Base;

  struct psm {
    
    enum code {
      normal_call = 0,
      no_lim_isr  = 32,
      no_gen_isr  = 64
    };

  };// end of struct psm

  inline psm::code 
  operator|(const psm::code &c1,const psm::code &c2)
  { return (psm::code)((int)c1|(int)c2); }
  inline psm::code 
  operator&(const psm::code &c1,const psm::code &c2)
  { return (psm::code)((int)c1&(int)c2); }


  class Phase_Space_Integrator;
  class Process_Base;

  class Phase_Space_Handler {
  private: 

    std::string m_name;

    Process_Integrator     *p_process, *p_active;
    Phase_Space_Integrator *p_integrator;  
    Cut_Data               *p_cuts;  

    Enhance_Observable_Base *p_enhancefunc;

    ATOOLS::Histogram *p_enhancehisto, *p_enhancehisto_current;

    static ATOOLS::Integration_Info *p_info;

    ATOOLS::Info_Key m_isrspkey, m_isrykey, m_isrxkey;
    ATOOLS::Info_Key m_beamspkey, m_beamykey;

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    Multi_Channel *p_fsrchannels, *p_isrchannels, *p_beamchannels;

    ATOOLS::Flavour_Vector p_flavours;
    ATOOLS::Vec4D_Vector   p_lab;

    ATOOLS::Poincare *p_massboost;

    int  m_nin, m_nout, m_nvec, m_fin_opt, m_dmode, m_enhancexs;
    bool m_initialized, m_sintegrator;

    long int m_maxtrials, m_killedpoints;

    double m_E, m_m[2], m_m2[2], m_s, m_smin, m_psweight, m_thkill;
    double m_result, m_error, m_abserror, m_enhance, m_osmass; 

    SHERPA::Variation_Weights * p_variationweights;

    std::vector<std::vector<double> > m_stats;

    psm::code m_cmode;

    bool m_printpspoint;

    double EnhanceFactor(Process_Base *const proc);

    bool MakeIncoming(ATOOLS::Vec4D *const p);
    bool Check4Momentum(const ATOOLS::Vec4D_Vector &p);

    void CheckSinglePoint();

    static void TestPoint(ATOOLS::Vec4D *const p,
			  ATOOLS::Vec4D_Vector cp,ATOOLS::Flavour_Vector fl,
			  const Subprocess_Info *info,size_t &n,
                          const ATOOLS::Mass_Selector* ms);

  public:

    //constructor
    Phase_Space_Handler(Process_Integrator *,double error=-1.);

    //destructor
    ~Phase_Space_Handler();

    // member functions
    void WriteOut(const std::string &path);
    bool ReadIn(const std::string &path,const size_t exclude=0);
    bool InitIncoming();
    bool CreateIntegrators();
    bool UpdateIntegrators();
    void InitCuts();
    
    double Integrate();
    double Differential();
    double Differential(Process_Integrator *const process,
			const psm::code mode=psm::normal_call);

    ATOOLS::Weight_Info *OneEvent(Process_Base *const proc,const int mode=0);

    void AddPoint(const double xs);

    void CalculateME();
    void CalculatePS();

    double Weight(ATOOLS::Vec4D_Vector &plab);

    static void TestPoint(ATOOLS::Vec4D *const p,const Process_Info *info,
                          const ATOOLS::Mass_Selector* ms,const int mode=0);
    static void TestPoint(ATOOLS::Vec4D *const p,
			  const size_t &nin,const size_t &nout,
			  const ATOOLS::Flavour_Vector &flavs,
			  const ATOOLS::Mass_Selector* ms);

    void MPISync();
    void Optimize();
    void EndOptimize();

    static ATOOLS::Integration_Info* GetInfo();

    static void DeleteInfo();

    void AddStats(const std::vector<double> &stats);

    // inline functions
    inline const ATOOLS::Flavour_Vector &Flavs() const { return p_flavours; }

    inline Cut_Data* Cuts() const { return p_cuts; }

    inline double   Error() const     { return m_error;     }
    inline long int MaxTrials() const { return m_maxtrials; }

    inline void SetFSRIntegrator(Multi_Channel *const fsr)
    { p_fsrchannels=fsr;  }
 
    inline Multi_Channel* BeamIntegrator() const  
    { return p_beamchannels; }
    inline Multi_Channel* ISRIntegrator() const   
    { return p_isrchannels;  }
    inline Multi_Channel* FSRIntegrator() const   
    { return p_fsrchannels;  }

    inline Process_Integrator* Process() const { return p_process; }
    inline Process_Integrator* Active() const  { return p_active;  }

    inline double PSWeight() const { return m_psweight; }
    inline double Enhance() const  { return m_enhance;  }

    inline void SetError(const double error) 
    { m_error=error; }
    inline void SetAbsError(const double error) 
    { m_abserror=error; }
    inline void SetMaxTrials(const long int maxtrials)
    { m_maxtrials=maxtrials; }
    void SetEnhanceObservable(const std::string &enhanceobs);
    void SetEnhanceFunction(const std::string &enhancefunc);
    void SetVariationWeights(SHERPA::Variation_Weights *vw)
    { p_variationweights=vw; }

    inline Phase_Space_Integrator* Integrator() const 
    { return p_integrator; }

    inline const std::vector<std::vector<double> > Stats() const
    { return m_stats; }

    inline void SetOSMass(const double &osmass) { m_osmass=osmass; }

    const ATOOLS::Vec4D_Vector &Momenta() { return p_lab; }

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
