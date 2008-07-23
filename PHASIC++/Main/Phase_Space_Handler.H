#ifndef Phase_Space_Handler_H
#define Phase_Space_Handler_H

#include "Cut_Data.H"
#include "Integrable_Base.H"
#include "Multi_Channel.H"
#include "CXXFLAGS.H"
#include <list>
#ifdef USING__Threading
#include <pthread.h>
#endif

namespace ATOOLS { 
  class Blob_Data_Base; 
  class Poincare;
}

namespace PHASIC {

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
    
    double        weight, procweight, xsecweight, xf1, xf2;
    unsigned long ntrial, xsecntrial;
    
    // constructor
    Weight_Info(const double _weight, 
		const double _procweight,
		const double _xsecweight,
                const double _xf1,
                const double _xf2,
		const unsigned long _ntrial,
		const unsigned long _xsecntrial):
      weight(_weight),
      procweight(_procweight),
      xsecweight(_xsecweight),
      xf1(_xf1), xf2(_xf2),
      ntrial(_ntrial),
      xsecntrial(_xsecntrial) {}

  };// end of struct Weight_Info

  class Phase_Space_Integrator;
  class Color_Integrator;
  class Helicity_Integrator;

  class Phase_Space_Handler {
  private: 

    std::string m_name;

    Integrable_Base        *p_process, *p_active;
    Phase_Space_Integrator *p_integrator;  
    ATOOLS::Cut_Data       *p_cuts;  

    static ATOOLS::Integration_Info *p_info;

    ATOOLS::Info_Key  m_isrspkey, m_isrykey, m_beamspkey, m_beamykey;
    ATOOLS::Info_Key  m_mu2key[2];

    BEAM::Beam_Spectra_Handler *p_beamhandler;
    PDF::ISR_Handler           *p_isrhandler;

    Multi_Channel *p_fsrchannels, *p_isrchannels, *p_beamchannels;
    std::list<std::string>* p_channellibnames;

    channelinfos m_beamparams, m_isrparams;

    ATOOLS::Flavour *p_flavours;
    ATOOLS::Vec4D   *p_cms, *p_lab;

    ATOOLS::Poincare *p_massboost;

    int  m_nin, m_nout, m_nvec, m_inttype, m_fin_opt;
    bool m_initialized, m_sintegrator;

    long int m_trials, m_maxtrials, m_sumtrials, m_events;

    double m_E, m_m[2], m_m2[2], m_s, m_smin, m_psweight;
    double m_result_1, m_result_2, m_flux, m_weight, m_error; 
    
    Color_Integrator    *p_colint;
    Helicity_Integrator *p_helint;

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

    bool MakeIncoming(ATOOLS::Vec4D *const p,const double mass=-1.);

    void TestIntegration(const double);
    void DropRedundantChannels();

    void Rotate(ATOOLS::Vec4D *const p,const size_t n=0);
    bool Check4Momentum(const ATOOLS::Vec4D *p);
    bool CompareCh(std::string C1,std::string C2);
    bool Compare(const ATOOLS::Vec4D *p1,const ATOOLS::Vec4D *p2);

    void MakeZChannels(const int type);

    bool MakeBeamChannels();
    bool MakeISRChannels();

    bool CreateBeamChannels();
    bool CreateISRChannels();

    bool LoadChannelLibraries(); 

    void CalculateME();
    void CalculatePS();

  public:

    //constructor
    Phase_Space_Handler(Integrable_Base *,PDF::ISR_Handler *,
			BEAM::Beam_Spectra_Handler *, double error=-1.);

    //destructor
    ~Phase_Space_Handler();

    // member functions
    Single_Channel * SetChannel(int,int,ATOOLS::Flavour *fl,
				std::string&,ATOOLS::Integration_Info *const);
    void WriteOut(const std::string &path,const bool force=false);
    bool ReadIn(const std::string &path,const size_t exclude=0);
    bool InitIncoming(const double mass=-1.);
    bool CreateIntegrators();
    bool UpdateIntegrators();
    void InitCuts();
    bool InitIntegrators();
    
    double Integrate();
    double Differential();
    double Differential(Integrable_Base *const process,
			const psm::code mode=psm::normal_call);

    ATOOLS::Blob_Data_Base *OneEvent(const double mass=-1.,const int mode=0);
    ATOOLS::Blob_Data_Base *SameEvent();

    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode=0);
    ATOOLS::Blob_Data_Base *SameWeightedEvent();

    void AddPoint(const double xs);
    void TestPoint(ATOOLS::Vec4D *const p);
    void TestPoint(ATOOLS::Vec4D *const p,int nin,int nout, ATOOLS::Flavour* flav);

    void ISRChannels(const int i,Channel_Info &ci) const; 
    void BeamChannels(const int i,Channel_Info &ci) const;

    static ATOOLS::Integration_Info *const GetInfo();

    static void DeleteInfo();

    void AddStats(const std::vector<double> &stats);

    // inline functions
    inline ATOOLS::Vec4D *const    Point() const { return p_lab;      }
    inline ATOOLS::Flavour *const  Flavs() const { return p_flavours; }
    inline ATOOLS::Cut_Data *const Cuts() const  { return p_cuts;     }

    inline double   Weight() const    { return m_weight;    }
    inline double   Error() const     { return m_error;     }
    inline long int MaxTrials() const { return m_maxtrials; }
 
    inline Multi_Channel *const BeamIntegrator() const  
    { return p_beamchannels; }
    inline Multi_Channel *const ISRIntegrator() const   
    { return p_isrchannels;  }
    inline Multi_Channel *const FSRIntegrator() const   
    { return p_fsrchannels;  }

    inline Single_Channel *const FSRIntegrator(const int i) 
    { return p_fsrchannels->Channel(i); }
    std::list<std::string>* GetChannelLibNames() 
    {return p_channellibnames;}

    inline int NumberOfBeamIntegrators() const  
    { return p_beamchannels->Number(); }
    inline int NumberOfISRIntegrators() const   
    { return p_isrchannels->Number();  }
    inline int NumberOfFSRIntegrators() const   
    { return p_fsrchannels->Number();  }

    inline channelinfos &BeamParameters() { return m_beamparams; } 
    inline channelinfos &ISRParameters()  { return m_isrparams;  } 

    inline Integrable_Base *const Process() const { return p_process; }
    inline Integrable_Base *const Active() const  { return p_active;  }

    inline double PSWeight() const { return m_psweight; }

    inline void SetError(const double error) 
    { m_error=error; }
    inline void SetMaxTrials(const long int maxtrials)
    { m_maxtrials=maxtrials; }

    inline Phase_Space_Integrator *const Integrator() const 
    { return p_integrator; }

    inline Color_Integrator    *ColorIntegrator() const    { return p_colint; }
    inline Helicity_Integrator *HelicityIntegrator() const { return p_helint; }

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
