#ifndef Phase_Space_Handler_h
#define Phase_Space_Handler_h

#include "Integrable_Base.H"
#include "Multi_Channel.H"

#include <sys/types.h>
#include <sys/stat.h>

namespace PHASIC  { class  Phase_Space_Integrator; }
namespace BEAM    { class  Beam_Spectra_Handler;   }
namespace PDF     { class  ISR_Handler;            }


namespace PHASIC {
  class Phase_Space_Handler {
    Integrable_Base            * proc;
    BEAM::Beam_Spectra_Handler * bh;
    PDF::ISR_Handler           * ih;
    Multi_Channel              * fsrchannels, * isrchannels, * beamchannels;
    channelinfos                 beam_params;
    channelinfos                 isr_params;
    int                          nin,nout,nvec;
    APHYTOOLS::Flavour         * psflavs;
    std::string                  name;
    int                          int_type;
    double                       error;
    long int                     maxtrials,sumtrials,events;
    double                       E,m1,m2,m12,m22,s,sprime,y;
    double                       sprimeB,yB,sprimeI,yI;
    double                       result1,result2,flux;
    double                       m_weight; 
    AMATOOLS::Vec4D            * p;
    Phase_Space_Integrator     * psi;    
    bool  m_initialized;

    void TestIntegration(double);
    bool MakeIncoming(AMATOOLS::Vec4D *);
    void Rotate(AMATOOLS::Vec4D *);
    bool Check4Momentum(AMATOOLS::Vec4D *);
    void DropRedundantChannels();
    bool Compare(AMATOOLS::Vec4D *,AMATOOLS::Vec4D *);
    bool MakeBeamChannels();
    bool MakeISRChannels();
    bool CreateBeamChannels();
    bool CreateISRChannels();
    bool InitIncomming();
  public:
    Phase_Space_Handler(Integrable_Base *,PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~Phase_Space_Handler();

    void   WriteOut(std::string);
    bool   ReadIn(std::string);
    bool   CreateIntegrators();

    double Integrate();
    double Differential();
    double Differential(Integrable_Base *);
    bool   OneEvent(int mode =0);
    bool   SameEvent();
    double WeightedEvent(int mode =0);
    double SameWeightedEvent();
    void   AddPoint(double);
    void   TestPoint(AMATOOLS::Vec4D *);

    // Simple access methods.

    inline AMATOOLS::Vec4D *     Point()        { return p; }
    inline APHYTOOLS::Flavour *  Flavs()        { return psflavs; }
    inline double                Weight()       { return m_weight; }


    inline Multi_Channel  * BeamIntegrator()          { return beamchannels; }
    inline Multi_Channel  * ISRIntegrator()           { return isrchannels; }
    inline Multi_Channel  * FSRIntegrator()           { return fsrchannels; }
    inline int              NumberOfBeamIntegrators() { return beamchannels->Number(); }
    inline int              NumberOfISRIntegrators()  { return isrchannels->Number(); }
    inline int              NumberOfFSRIntegrators()  { return fsrchannels->Number(); }
    inline Single_Channel * FSRIntegrator(int i)      { return fsrchannels->Channel(i); }
    inline channelinfos   & BeamParameters()          { return beam_params; } 
    inline channelinfos   & ISRParameters()           { return isr_params; } 

    void ISRChannels(int i,Channel_Info & ci) {
      if (i<isr_params.size()) {
	ci.type       = isr_params[i].type;
	ci.parameters = isr_params[i].parameters;
	return;
      }
      else {
	AORGTOOLS::msg.Error()<<"Error in Phase_Space_Handler::ISRChannels("<<i<<")"<<std::endl
			      <<"  delimiter out of bounds."<<std::endl;
	abort();
      }
    }

    void BeamChannels(int i,Channel_Info & ci) {
      if (i<beam_params.size()) {
	ci.type       = beam_params[i].type;
	ci.parameters = beam_params[i].parameters;
	return;
      }
      else {
	AORGTOOLS::msg.Error()<<"Error in Phase_Space_Handler::BeamChannels("<<i<<")"<<std::endl
			      <<"  delimiter out of bounds."<<std::endl;
	abort();
      }
    }
  };
}

#endif
