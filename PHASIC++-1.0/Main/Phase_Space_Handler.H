#ifndef Phase_Space_Handler_h
#define Phase_Space_Handler_h

#include "Multi_Channel.H"

namespace PHASIC  { class  Phase_Space_Integrator; }
namespace PHASIC  { class  Phase_Space_Generator;  }
namespace PHASIC  { const int IS_XS_FLAG = -255;   }
namespace BEAM    { class  Beam_Handler;           }
namespace ISR     { class  ISR_Handler;            }
namespace AMEGIC  { class  Process_Base;           }                                       
namespace EXTRAXS { class  XS_Base;                }                                       


namespace PHASIC {
  class Phase_Space_Handler {
    AMEGIC::Process_Base * proc;
    Multi_Channel * fsrchannels;
    Multi_Channel * isrchannels;
    Multi_Channel * beamchannels;
    ISR::ISR_Handler   * ih;
    BEAM::Beam_Handler * bh;
    std::vector<Channel_Info> beam_params;
    std::vector<Channel_Info> isr_params;
    Phase_Space_Generator * psgen;
    int                  nin,nout,nvec;
    APHYTOOLS::Flavour * psflavs;
    std::string          name;
    int         int_type;
    double      error;
    long int    maxtrials,sumtrials,events;
    // Kinematics             * kin;
    double                   E,m1,m2,m12,m22,s,sprime,y;
    double                   sprimeB,yB,sprimeI,yI;
    double                   result1,result2,flux;
    AMATOOLS::vec4d        * p;
    Phase_Space_Integrator * psi;    
    bool MakeBeamChannels();
    bool MakeISRChannels();
    bool MakeFSRChannels();
    void   AddChannels(AMEGIC::Process_Base *,Multi_Channel *,
		       std::vector<Channel_Info> &,std::vector<Channel_Info> &);
    void DropRedundantChannels();
    void TestIntegration(double);
    bool Compare(AMATOOLS::vec4d *,AMATOOLS::vec4d *);
    bool MakeIncoming(AMATOOLS::vec4d *);
    void Rotate(AMATOOLS::vec4d *);


    void Init(APHYTOOLS::Flavour *);
    bool Check4Momentum(AMATOOLS::vec4d *);
  public:
    Phase_Space_Handler(AMEGIC::Process_Base *,ISR::ISR_Handler *,BEAM::Beam_Handler *);
    ~Phase_Space_Handler();
    bool   CreateChannelLibrary(std::string,std::string);
    void   CollectChannels();
    bool   CreateIntegrators();

    void   WriteOut(std::string);
    bool   ReadIn(std::string);
    double Integrate();
    double Differential();
    double Differential(AMEGIC::Process_Base *);
    double Differential(EXTRAXS::XS_Base *);
    bool   OneEvent(int mode =0);
    bool   SameEvent();
    double WeightedEvent();
    void AddPoint(double);
    void TestPoint(AMATOOLS::vec4d *);

    // Simple access methods.

    inline AMATOOLS::vec4d *     Point()        { return p; }
    inline APHYTOOLS::Flavour *  Flavs()        { return psflavs; }


    inline Multi_Channel  * BeamIntegrator()          { return beamchannels; }
    inline Multi_Channel  * ISRIntegrator()           { return isrchannels; }
    inline Multi_Channel  * FSRIntegrator()           { return fsrchannels; }
    inline int              NumberOfBeamIntegrators() { return beamchannels->Number(); }
    inline int              NumberOfISRIntegrators()  { return isrchannels->Number(); }
    inline int              NumberOfFSRIntegrators()  { return fsrchannels->Number(); }
    inline Single_Channel * FSRIntegrator(int i)      { return fsrchannels->Channel(i); }

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

// #ifndef Phase_Space_Handler_h
// #define Phase_Space_Handler_h

// #include "Multi_Channel.H"

// namespace PHASIC  { class  Phase_Space_Integrator; }
// namespace PHASIC  { class  Phase_Space_Generator;  }
// namespace PHASIC  { const int IS_XS_FLAG = -255;   }
// namespace BEAM    { class  Beam_Handler;           }
// namespace ISR     { class  ISR_Handler;            }
// namespace AMEGIC  { class  Process_Base;           }                                       
// namespace EXTRAXS { class  XS_Base;                }                                       


// namespace PHASIC {
//   /*!
//     This class manages the interaction of individual processes or process groups
//     with the integration methods. The idea behind is that these processes don't know
//     about either the generation of integration channels nor the integration procedure.
//     Therefore, the two principal classes Phase_Space_Generator and Phase_Space_Integrator
//     are invisible for processes and can be accessed only via the Phase_Space_Handler.
//     In turn this results in the existence of a specific Phase_Space_Handler for
//     each Single_Process or Process_Group. 
//     For Process_Groups, the Phase_Space_Handlers of subordinate Process_Groups or 
//     Single_Processes are erased at the end of the initialization procedure.
//   */
//   class Phase_Space_Handler {
//     /*!
//       A pointer to the process to be integrated. This will be handed over in the constructor.
//     */
//     AMEGIC::Process_Base * proc;
//     //    EXTRAXS::XS_Base     * xs;
//     /*!
//       There're three /sa Multi_Channels for the phase space integration. One is for the
//       integration over final state momenta in the c.m. frame of the hard process.
//       There the channels are produced amplitude by amplitude and added. One
//       Multi_Channel is for the integration over the initial state - basically
//       x1 and x2 or, better, s' and y, the (reduced) c.m. energy of the hard process
//       and the rapidity of the c.m. system w.r.t. the lab frame. The last one is for the
//       set-up of the beams. Similarly to the isrchannels it integrates over x1 and x2 or -
//       more conveniently, over s' and y. 
//     */
//     Multi_Channel * fsrchannels;
//     Multi_Channel * isrchannels;
//     Multi_Channel * beamchannels;
//     /*!
//       Similarly the ISR_Handler, carrying the information about the ISR strategy, is handed over 
//       from the Process_Base to the Phase_Space_Handler in its constructor. At the moment ISR 
//       strategies are available that rely on combinations of Structure_Function and No_ISR. 

//       For the Yennie-Frautschi-Suura method, however, it is not clear to me how to proceed.
//     */
//     ISR::ISR_Handler   * ih;
//     BEAM::Beam_Handler * bh;
//     /*!
//       The parameters of eventual initial state radiation channels. isr_types labels all 
//       possible /sa isrchannels (s' will be the actual c.m. energy squared), but this depends on
//       the actual ISR_Handler. For the Structure_Function approach we have:
//       0 = a simple pole according to 1/s'
//       1 = a breit-wigner wih paramters given by the according isr_masses and isr_widths
//       3 = a typical leptonic 1/(s-s') pole
//     */
//     std::vector<Channel_Info> beam_params;
//     std::vector<Channel_Info> isr_params;
//     /*!
//       The Phase_Space_Generator object is instantiated and owned by the Phase_Space_Handler.
//       It is responsible for the creation of channel libraries, their usage and the set-up
//       of a process specific integrator.
//     */   
//     Phase_Space_Generator * psgen;
//     /*! 
//       Process-specifics : Number of legs, flavours, name. All these parameters are 
//       filled with the knowledge given by the process, so to some extend the same 
//       information exists twice.
//     */
//     int                  nin,nout,nvec;
//     APHYTOOLS::Flavour * psflavs;
//     std::string          name;
//     /*! 
//       Phase space specifics : 
//       - integration method, potential jetfinder, on-the-flight analysis switch
//       - integration error, and jetmeasure
//     */
//     int         int_type;
//     double      error;
//     long int    maxtrials,sumtrials,events;
//     /*!
//       The actual integration :
//       - Parameters of the incoming beams
//       - Parameters of the initial state
//       - Four-momenta to be generated
//       - Phase_Space_Integrator for the actual integration
//     */
//     // Kinematics             * kin;
//     double                   E,m1,m2,m12,m22,s,sprime,y;
//     double                   sprimeB,yB,sprimeI,yI;
//     double                   result1,result2,flux;
//     AMATOOLS::vec4d        * p;
//     Phase_Space_Integrator * psi;    
//     /*!
//       Method used to analyse the individual (Final-State-) channels and to construct 
//       ISR channels accordingly. The idea is to fill the classification vectors for the
//       different ISR channels, as exemplified by the Structure_Function case: 
//       By default a simple 1/s' pole with isr_type = 0 and mass and width = 0 is added. 
//       If the process has a leptonic initial state, isr_type = 3 for the typical 1/(s-s') 
//       pole is added, again with mass and width equal to zero. Then all active final 
//       state channels are iterated over and type, mass and width are filled with the 
//       method ISRInfo (/sa Channel_Base).
//       Essentially this method fills the parameters via the method ISRType for each channel.
//       Comparing with all other elements in the vectors so far, the corresponding
//       characteristics, i.e. type, mass, and idth, are added to the vectors, if they
//       give rise to a new channel. At the moment this is only true for type = 1,
//       i.e. Breit-Wigner rtesonances in the s-channel.
//       This method is called by CreateIntegrator of the Phase_Space_Handler
//       and calling the method CreateISRChannels of the Phase_Space_Generator.
//     */
//     bool MakeBeamChannels();
//     bool MakeISRChannels();
//     bool MakeFSRChannels();
//     /*!
//       The following method is used when setting up the integrator of a Process_Group.
//       It recursively adds all channels (via a new instantiation) of all processes that 
//       have no partner. With this condition the number of channels is at least a bit reduced 
//       to those that are potentially different. A similar thing happens for the initial state part. 
//       Iterating over the same processes, types, masses and widths are added to corresponding vectors.
//     */
//     void   AddChannels(AMEGIC::Process_Base *,Multi_Channel *,
// 		       std::vector<Channel_Info> &,std::vector<Channel_Info> &);
//     /*!
//       Method to drop channels from the multichannel. The method is to
//       generate on phase space ponit and to check it. 
//       Then dropping might happen for a number of reasons :
//       - a channel produces a zero weight.
//       - if two channels produce the same momenta with the same set of
//         random numbers /sa Compare
//       - for the number of resonating propagators : this is commented out
//     */	
//     void DropRedundantChannels();
//     /*!
//       A cheap and simple test integration of channels.
//     */
//     void TestIntegration(double);
//     /*!
//       Compares whether two sets of momenta are identical (up to permutations).
//       If so a "true" is returned.
//     */
//     bool Compare(AMATOOLS::vec4d *,AMATOOLS::vec4d *);
//     /*!
//       Initializes the initial state part of the process:
//       - defines default (no ISR) initial state vectors.
//     */
//     bool MakeIncoming(AMATOOLS::vec4d *);
//     /*!
//       This method might be used for the usage of the integration with ISR being switched on.
//       Let us assume, we have, say, p pbar ->2 jets to be integrated. Obviuosly, there are
//       subprocesses like u dbar -> u dbar, and obviously the paton density of the two quarks
//       is different for proton and anti-proton. A nice way to have the u-quark stemming from
//       both the p and the pbar is to use the same value for the matgrix element but two different
//       sets of parton densities (from left and right). This basically boils down to rotating
//       the full by 180 degrees to ensure some non-symmetric cuts. In fact that's why we may need
//       this method in the future.
//     */
//     void Rotate(AMATOOLS::vec4d *);


//     void Init(APHYTOOLS::Flavour *);
//     bool Check4Momentum(AMATOOLS::vec4d *);
//   public:
//     /*!
//       Explicit constructor for the PS_Handler of a Single_Process or
//       a process group. In the latter case, all outgoing flavours are
//       Flavour(kf::none) :
//       - Reads phase space data : ycut, error, jetfinder, integration method.
//       - Initializes a set of four-vectors.
//       - Creates a specific, process-dependent selector 
//       - Initializes the initial state, incoming vectors, etc..
//       - makes a safety copy of the flavours.
//       /sa CreateSelector(), InitInitialState()
//     */
//     Phase_Space_Handler(AMEGIC::Process_Base *,ISR::ISR_Handler *,BEAM::Beam_Handler *);
//     Phase_Space_Handler(EXTRAXS::XS_Base *,ISR::ISR_Handler *,BEAM::Beam_Handler *);
//     ~Phase_Space_Handler();
//     /*!
//       Creates the channels of the multichannel with help of the
//       Phase_Space_Generator. The latter is instantiated in this method.
//       Then its method Construct is called for all channels. This method checks whether
//       a corresponding file in rthe directory "ptype" with characteristic "pID" is
//       available. If not, it is created and Construct returns a "true".
//       If no new channels were produced, CreateChannelLibrary returns a "true", otherwise
//       it returns a "false".
//     */
//     bool   CreateChannelLibrary(std::string,std::string);
//     /*!
//       The following method is used when setting up the integrator of a Process_Group.
//       It adds all channels of all processes that have no partner via the method AddChannels that
//       is fully recursive.
//     */
//     void   CollectChannels();
//     /*!
//       Initializes the specific integrator to be used in this integration
//       of a Single_Process. This can be a single channel like Rambo or Sarge, 
//       or a combination of both, or any combination of the explicit multichannel 
//       with or without them. If there are two incoming particles, and if
//       isr is to be taken into account, according ISR channels are to be created.
//     */
//     bool   CreateIntegrators();

//     void   WriteOut(std::string);
//     bool   ReadIn(std::string);
//     /*!
//       Here the actual integration takes place, essentially only
//       the Phase_Space_Integrator is called after setting the ISR-range.
//     */
//     double Integrate();
//     /*!
//       A wrapper for the other Differential method when used for actual integrations.
//     */
//     double Differential();
//     /*!
//       The actual value of the cross section is calculated here :
//       - First, the incoming beam is set up (via an integration over sprimeB and yB,
//         performed by the beamchannels.
//       - Following this the initial state characteristics (sprimeI and yI) are determined with 
//         the isr integrator, following that suitable initial state vectors are constructed.
//       	Here, InitInitialState of the ISR_Handler is used, the use is quite limited to
// 	processes with particles that have identical masses.
//       - Then, the fsr integrator determines the outgoing momenta.
//       - Before selection cuts etc. are applied, the full system is boosted into the
//         lab system. It should be stressed that boosts might be mandatory in the selector
// 	as well. As an example consider the jet-finder, where for the case of e+e- 
// 	collisions and the Durham algorithm the relevant frame is always the c.m.
// 	frame of the hard process or the case of DIS, where the frame in question should
// 	be the photon's Breit frame.
//       - If the selector is passed, i.e. if all momenta are inside the valid region,
//         a scale is determined which by default is the hard c.m. energy squared.
// 	Eventually this is replaced, in any case, initial PDF's are precalculated if needed.
//       - Then the vectors are boosted back into the hard c.m. frame and the value of the
//         matrix element squared is evaluated. 
//       - Finally, the weights of the three integrators are calculated. Their product, multiplied
//         with the value of the ME squared - obtained via Differential - is the result for 
// 	parton1(2) coming from beam1(2). 
//       - Essentially, after flipping the momenta, a similar procedure (calculation of PDF's,
//         folding in the ME's) is repeated via Differential2 and yields the total result  
// 	that is returned to the Phase_Space_Integrator. 
//     */
//     double Differential(AMEGIC::Process_Base *);
//     double Differential(EXTRAXS::XS_Base *);
//     /*! 
//       Generates a one-weighted event for a Single_Process or a (non-atomar)
//       Process_Group. The alogrithm is as follows:
//       - Select a process within the group (according to the relative xsec).
//       - For up to maxtrials trials generate momentuim configurations, including ISR.
//       - Check the configuration whether or not it passes the cuts, if not, try again.
//       - Generate the differential cross section at this phase space point and the
//         corresponding phase space weight, again including ISR-weights.
//       - Check if the product of weight and differential cross section is larger
//         than the maximum accumulated before during the integration phase.
//         If this is the case, set the new maximum accordingly and stop the
//         event generation with a success.
//       - If the combined value is smaller than the maximum, accept or reject
//         w.r.t. this maximum.
//     */
//     bool   OneEvent(int mode =0);
//     bool   SameEvent();
//     /*!
//       Generates weighted events for a Single_Process or a (non-atomar)
//       Process_Group. The alogrithm is as the one for one-weighted events
//       with the only difference that no rejection according to the maximum
//       is performed. Therefore, no loop of trial is needed.
//     */
//     double WeightedEvent();
//     /*!
//       Here, the differential cross sections are added up sub process-wise.
//      */
//     void AddPoint(double);
//     /*!
//       With help of Rambo, one set of momenta is created that is used
//       to check the amplitude for string and gauge invariance.
//       /sa Single_Process::Tests
//     */
//     void TestPoint(AMATOOLS::vec4d *);

//     // Simple access methods.

//     inline AMATOOLS::vec4d *     Point()        { return p; }
//     inline APHYTOOLS::Flavour *  Flavs()        { return psflavs; }


//     inline Multi_Channel  * BeamIntegrator()          { return beamchannels; }
//     inline Multi_Channel  * ISRIntegrator()           { return isrchannels; }
//     inline Multi_Channel  * FSRIntegrator()           { return fsrchannels; }
//     inline int              NumberOfBeamIntegrators() { return beamchannels->Number(); }
//     inline int              NumberOfISRIntegrators()  { return isrchannels->Number(); }
//     inline int              NumberOfFSRIntegrators()  { return fsrchannels->Number(); }
//     inline Single_Channel * FSRIntegrator(int i)      { return fsrchannels->Channel(i); }

//     void ISRChannels(int i,Channel_Info & ci) {
//       if (i<isr_params.size()) {
// 	ci.type       = isr_params[i].type;
// 	ci.parameters = isr_params[i].parameters;
// 	return;
//       }
//       else {
// 	AORGTOOLS::msg.Error()<<"Error in Phase_Space_Handler::ISRChannels("<<i<<")"<<std::endl
// 			      <<"  delimiter out of bounds."<<std::endl;
// 	abort();
//       }
//     }

//     void BeamChannels(int i,Channel_Info & ci) {
//       if (i<beam_params.size()) {
// 	ci.type       = beam_params[i].type;
// 	ci.parameters = beam_params[i].parameters;
// 	return;
//       }
//       else {
// 	AORGTOOLS::msg.Error()<<"Error in Phase_Space_Handler::BeamChannels("<<i<<")"<<std::endl
// 			      <<"  delimiter out of bounds."<<std::endl;
// 	abort();
//       }
//     }
//   };
// }
// #endif
