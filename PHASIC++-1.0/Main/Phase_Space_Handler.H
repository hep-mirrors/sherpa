#ifndef Phase_Space_Handler_h
#define Phase_Space_Handler_h

#include "Cut_Data.H"
#include "Integrable_Base.H"
#include "Multi_Channel.H"

#include <sys/types.h>
#include <sys/stat.h>

namespace PHASIC  { class  Phase_Space_Integrator; }
namespace BEAM    { class  Beam_Spectra_Handler;   }
namespace PDF     { class  ISR_Handler;            }

namespace ATOOLS {
  class Blob_Data_Base;
}

namespace PHASIC {
  struct Weight_Info {
    double        weight;
    unsigned long ntrial;
    Weight_Info(const double _weight, const unsigned long _ntrial):
      weight(_weight), ntrial(_ntrial) {}
  };
  

  class Phase_Space_Handler {
    Integrable_Base            * proc;
    BEAM::Beam_Spectra_Handler * bh;
    PDF::ISR_Handler           * ih;
    Multi_Channel              * fsrchannels, * isrchannels, * beamchannels;
    channelinfos                 beam_params;
    channelinfos                 isr_params;
    int                          nin,nout,nvec;
    ATOOLS::Flavour            * psflavs;
    std::string                  name;
    int                          int_type;
    double                       error;
    long int                     maxtrials,sumtrials,events;
    double                       E,m1,m2,m12,m22,s,sprime,y;
    double                       sprimeB,yB,sprimeI,yI;
    double                       result1,result2,flux;
    double        m_weight; 
    unsigned long m_ntrial;
    ATOOLS::Vec4D          * p;
    Phase_Space_Integrator * psi;    
    bool  m_initialized;
    ATOOLS::Integration_Info   * p_info;
    ATOOLS::Info_Key             m_spkey, m_ykey, m_xkey;
    ATOOLS::Cut_Data           * p_cuts;

    void TestIntegration(double);
    bool MakeIncoming(ATOOLS::Vec4D *,double _mass=-1.);
    void Rotate(ATOOLS::Vec4D *);
    bool Check4Momentum(ATOOLS::Vec4D *);
    void DropRedundantChannels();
    bool Compare(ATOOLS::Vec4D *,ATOOLS::Vec4D *);
    bool MakeBeamChannels();
    bool MakeISRChannels();
    bool CreateBeamChannels();
    bool CreateISRChannels();
  public:
    Phase_Space_Handler(Integrable_Base *,PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~Phase_Space_Handler();

    void   WriteOut(std::string);
    bool   ReadIn(std::string);
    bool   CreateIntegrators();

    double Integrate();
    double Differential();
    double Differential(Integrable_Base *);
    bool   OneEvent(double _mass = -1.,int mode = 0);
    bool   SameEvent();
    ATOOLS::Blob_Data_Base *  WeightedEvent(int mode =0);
    ATOOLS::Blob_Data_Base *  SameWeightedEvent();
    void   AddPoint(double);
    void   TestPoint(ATOOLS::Vec4D *);
    bool   InitIncoming(double _mass=-1.);

    // Simple access methods.

    inline ATOOLS::Vec4D *   Point()  { return p; }
    inline ATOOLS::Flavour * Flavs()  { return psflavs; }
    inline double            Weight() { return m_weight; }
    inline unsigned long     NTrial() { return m_ntrial; }

    inline Multi_Channel  * BeamIntegrator()          { return beamchannels; }
    inline Multi_Channel  * ISRIntegrator()           { return isrchannels; }
    inline Multi_Channel  * FSRIntegrator()           { return fsrchannels; }
    inline int              NumberOfBeamIntegrators() { return beamchannels->Number(); }
    inline int              NumberOfISRIntegrators()  { return isrchannels->Number(); }
    inline int              NumberOfFSRIntegrators()  { return fsrchannels->Number(); }
    inline Single_Channel * FSRIntegrator(int i)      { return fsrchannels->Channel(i); }
    inline channelinfos   & BeamParameters()          { return beam_params; } 
    inline channelinfos   & ISRParameters()           { return isr_params; } 

    inline Integrable_Base *const Process() const { return proc; }

    inline double Error() const { return error; }

    void ISRChannels(int i,Channel_Info & ci) {
      if (i<(int)isr_params.size()) {
	ci.type       = isr_params[i].type;
	ci.parameters = isr_params[i].parameters;
	return;
      }
      else {
	ATOOLS::msg.Error()<<"Error in Phase_Space_Handler::ISRChannels("<<i<<")"<<std::endl
			      <<"  delimiter out of bounds."<<std::endl;
	abort();
      }
    }

    void BeamChannels(int i,Channel_Info & ci) {
      if (i<(int)beam_params.size()) {
	ci.type       = beam_params[i].type;
	ci.parameters = beam_params[i].parameters;
	return;
      }
      else {
	ATOOLS::msg.Error()<<"Error in Phase_Space_Handler::BeamChannels("<<i<<")"<<std::endl
			      <<"  delimiter out of bounds."<<std::endl;
	abort();
      }
    }
  };
}

#endif
