#ifndef Phase_Space_Handler_H
#define Phase_Space_Handler_H

#include "Cut_Data.H"
#include "Integrable_Base.H"
#include "Multi_Channel.H"
#include <list>

namespace ATOOLS { class Blob_Data_Base; }

namespace PHASIC {

  struct Weight_Info {
    
    double        weight;
    unsigned long ntrial;
    
    // constructor
    Weight_Info(const double _weight, 
		const unsigned long _ntrial):
      weight(_weight), 
      ntrial(_ntrial) {}

  };// end of struct Weight_Info

  class Phase_Space_Integrator;
  class PI_Interface;

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

    Multi_Channel *p_fsrchannels, *p_zchannels, *p_kpchannels;
    Multi_Channel *p_isrchannels, *p_beamchannels;
    std::list<std::string>* p_channellibnames;

    channelinfos m_beamparams, m_isrparams, m_zparams, m_kpparams;

    ATOOLS::Flavour *p_flavours;
    ATOOLS::Vec4D   *p_cms, *p_lab;

    int  m_nin, m_nout, m_nvec, m_inttype, m_fin_opt, m_use_pi;
    bool m_initialized;

    long int m_trials, m_maxtrials, m_sumtrials, m_events;

    double m_E, m_m[2], m_m2[2], m_s, m_smin, m_psweight;
    double m_result_1, m_result_2, m_flux, m_weight, m_error; 
    
    PI_Interface* p_pi;

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
    bool MakeKMRChannels();

    bool CreateBeamChannels();
    bool CreateISRChannels();
    bool CreateKMRChannels();

    bool LoadChannelLibraries(); 
  public:

    //constructor
    Phase_Space_Handler(Integrable_Base *,PDF::ISR_Handler *,
			BEAM::Beam_Spectra_Handler *);

    //destructor
    ~Phase_Space_Handler();

    // member functions
    Single_Channel * SetChannel(int,int,ATOOLS::Flavour *fl,
				std::string&,ATOOLS::Integration_Info *const);
    void WriteOut(const std::string &path,const bool force=false);
    bool ReadIn(const std::string &path,const size_t exclude=0);
    bool InitIncoming(const double mass=-1.);
    bool CreateIntegrators();
    void InitCuts();
    
    double Integrate();
    double Differential();
    double Differential(Integrable_Base *const process,const int mode=0);

    bool OneEvent(const double mass=-1.,const int mode=0);
    bool SameEvent();

    ATOOLS::Blob_Data_Base *WeightedEvent(const int mode =0);
    ATOOLS::Blob_Data_Base *SameWeightedEvent();

    void AddPoint(const double xs);
    void TestPoint(ATOOLS::Vec4D *const p);

    void ISRChannels(const int i,Channel_Info &ci) const; 
    void BeamChannels(const int i,Channel_Info &ci) const;

    bool CreatePI();

    static ATOOLS::Integration_Info *const GetInfo();

    static void DeleteInfo();

    // inline functions
    inline ATOOLS::Vec4D *const    Point() const { return p_lab;      }
    inline ATOOLS::Flavour *const  Flavs() const { return p_flavours; }
    inline ATOOLS::Cut_Data *const Cuts() const  { return p_cuts;     }

    inline double   Weight() const    { return m_weight;    }
    inline double   Error() const     { return m_error;     }
    inline long int MaxTrials() const { return m_maxtrials; }
    inline int      UsePI() const   { return m_use_pi;  }
 
    inline Multi_Channel *const BeamIntegrator() const  
    { return p_beamchannels; }
    inline Multi_Channel *const KMRKPIntegrator() const 
    { return p_kpchannels;   }
    inline Multi_Channel *const KMRZIntegrator() const  
    { return p_zchannels;    }
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
    inline int NumberOfKMRKPIntegrators() const 
    { return p_kpchannels->Number();   }
    inline int NumberOfKMRZIntegrators() const  
    { return p_zchannels->Number();    }
    inline int NumberOfFSRIntegrators() const   
    { return p_fsrchannels->Number();  }

    inline channelinfos &BeamParameters() { return m_beamparams; } 
    inline channelinfos &ISRParameters()  { return m_isrparams;  } 
    inline channelinfos &ZParameters()    { return m_zparams;    } 
    inline channelinfos &KPParameters()   { return m_kpparams;   } 

    inline Integrable_Base *const Process() const { return p_process; }
    inline Integrable_Base *const Active() const  { return p_active;  }

    inline double PSWeight() const { return m_psweight; }

    inline void SetError(const double error) 
    { m_error=error; }
    inline void SetMaxTrials(const long int maxtrials)
    { m_maxtrials=maxtrials; }
    inline void SetUsePI(const int use_pi) 
    { m_use_pi=use_pi; }

    inline Phase_Space_Integrator *const Integrator() const 
    { return p_integrator; }

    inline PI_Interface *const PI() const { return p_pi; }

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
