#ifndef Phase_Space_Handler_h
#define Phase_Space_Handler_h

#include "Multi_Channel.H"

namespace PHASIC  { class  Phase_Space_Integrator; }
namespace PHASIC  { const int IS_XS_FLAG = -255;   }
namespace BEAM    { class  Beam_Spectra_Handler;   }
namespace ISR     { class  ISR_Handler;            }
namespace AMEGIC  { class  Process_Base;           }                                       
namespace EXTRAXS { class  XS_Base;                }                                       


namespace PHASIC {
  class Phase_Space_Handler {
    AMEGIC::Process_Base       * proc;
    Multi_Channel              * fsrchannels;
    Multi_Channel              * isrchannels;
    Multi_Channel              * beamchannels;
    BEAM::Beam_Spectra_Handler * bh;
    ISR::ISR_Handler           * ih;
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
    AMATOOLS::Vec4D            * p;
    Phase_Space_Integrator     * psi;    

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

  public:
    Phase_Space_Handler(AMEGIC::Process_Base *,ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *);
    ~Phase_Space_Handler();

    void   WriteOut(std::string);
    bool   ReadIn(std::string);
    bool   CreateIntegrators();

    double Integrate();
    double Differential();
    double Differential(AMEGIC::Process_Base *);
    double Differential(EXTRAXS::XS_Base *);
    bool   OneEvent(int mode =0);
    bool   SameEvent();
    double WeightedEvent();
    void   AddPoint(double);
    void   TestPoint(AMATOOLS::Vec4D *);

    // Simple access methods.

    inline AMATOOLS::Vec4D *     Point()        { return p; }
    inline APHYTOOLS::Flavour *  Flavs()        { return psflavs; }


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

  /*! 
    \file 
    \brief contains the class PHASIC::Phase_Space_Handler
  */

  /*!
    \class Phase_Space_Handler
    \brief This class manages the interaction of individual processes or process groups
    with the integration methods.

    This class manages the interaction of individual processes or process groups
    with the integration methods. The idea behind is that these processes don't know
    about either the generation of integration channels nor the integration procedure.
    Therefore, the two principal classes Phase_Space_Generator and Phase_Space_Integrator
    are invisible for processes and can be accessed only via the Phase_Space_Handler.
    In turn this results in the existence of a specific Phase_Space_Handler for
    each Single_Process or Process_Group. 

    For Process_Groups, the Phase_Space_Handlers of subordinate Process_Groups or 
    Single_Processes are erased at the end of the initialization procedure.
  */

  /*!
    \var AMEGIC::Process_Base * proc;
    \brief   A pointer to the process to be integrated. This will be handed over in the constructor.
  */
  /*!
     \var Multi_Channel * fsrchannels;
     \var Multi_Channel * isrchannels;
     \var Multi_Channel * beamchannels;

     There're three /sa Multi_Channels for the phase space integration. One is for the
     integration over final state momenta in the c.m. frame of the hard process.
     There the channels are produced amplitude by amplitude and added. One
     Multi_Channel is for the integration over the initial state - basically
     x1 and x2 or, better, s' and y, the (reduced) c.m. energy of the hard process
     and the rapidity of the c.m. system w.r.t. the lab frame. The last one is for the
     set-up of the beams. Similarly to the isrchannels it integrates over x1 and x2 or -
     more conveniently, over s' and y. 
  */

  /*!
     \var ISR::ISR_Handler   * ih;
     \var BEAM::Beam_Handler * bh;

    Similarly the ISR_Handler, carrying the information about the ISR strategy, is handed over 
    from the Process_Base to the Phase_Space_Handler in its constructor. At the moment ISR 
    strategies are available that rely on combinations of Structure_Function and No_ISR. 

    For the Yennie-Frautschi-Suura method, however, it is not clear to me how to proceed.
  */
  /*!
    \var std::vector<Channel_Info> beam_params;
    \var std::vector<Channel_Info> isr_params;
    The parameters of eventual initial state radiation channels. isr_types labels all 
    possible /sa isrchannels (s' will be the actual c.m. energy squared), but this depends on
    the actual ISR_Handler. For the Structure_Function approach we have:
     0 = a simple pole according to 1/s'
     1 = a breit-wigner wih paramters given by the according isr_masses and isr_widths
     3 = a typical leptonic 1/(s-s') pole
  */
  /*!
    \var Phase_Space_Generator * psgen;
      The Phase_Space_Generator object is instantiated and owned by the Phase_Space_Handler.
      It is responsible for the creation of channel libraries, their usage and the set-up
      of a process specific integrator.
  */   
     
  /*! 
    \var int                  nin,nout,nvec;
    \var APHYTOOLS::Flavour * psflavs;
    \var std::string          name;
      Process-specifics : Number of legs, flavours, name. All these parameters are 
      filled with the knowledge given by the process, so to some extend the same 
       information exists twice.
  */
  /*! 
    \var int         int_type;
    \var double      error;
    \var long int    maxtrials,sumtrials,events;

      Phase space specifics : 
       - integration method, potential jetfinder, on-the-flight analysis switch
       - integration error, and jetmeasure
   */
   /*!
     \var double                   E,m1,m2,m12,m22,s,sprime,y;
     \var double                   sprimeB,yB,sprimeI,yI;
     \var double                   result1,result2,flux;
     \var AMATOOLS::Vec4D        * p;
     \var Phase_Space_Integrator * psi;    

     The actual integration :
       - Parameters of the incoming beams
       - Parameters of the initial state
       - Four-momenta to be generated
       - Phase_Space_Integrator for the actual integration
   */

  /*!
    \fn bool Phase_Space_Handler::MakeBeamChannels();
    \fn bool Phase_Space_Handler::MakeISRChannels();
    \fn bool Phase_Space_Handler::MakeFSRChannels();
    
    Method used to analyse the individual (Final-State-) channels and to construct 
    ISR channels accordingly. The idea is to fill the classification vectors for the
    different ISR channels, as exemplified by the Structure_Function case: 
    By default a simple 1/s' pole with isr_type = 0 and mass and width = 0 is added. 
    If the process has a leptonic initial state, isr_type = 3 for the typical 1/(s-s') 
    pole is added, again with mass and width equal to zero. Then all active final 
    state channels are iterated over and type, mass and width are filled with the 
    method ISRInfo (/sa Channel_Base).
    Essentially this method fills the parameters via the method ISRType for each channel.
    Comparing with all other elements in the vectors so far, the corresponding
    characteristics, i.e. type, mass, and idth, are added to the vectors, if they
    give rise to a new channel. At the moment this is only true for type = 1,
    i.e. Breit-Wigner rtesonances in the s-channel.
    This method is called by CreateIntegrator of the Phase_Space_Handler
    and calling the method CreateISRChannels of the Phase_Space_Generator.
  */
  /*!
    \fn  void   Phase_Space_Handler::AddChannels(AMEGIC::Process_Base *,Multi_Channel *,
 		       std::vector<Channel_Info> &,std::vector<Channel_Info> &);
      
    The following method is used when setting up the integrator of a Process_Group.
    It recursively adds all channels (via a new instantiation) of all processes that 
    have no partner. With this condition the number of channels is at least a bit reduced 
    to those that are potentially different. A similar thing happens for the initial state part. 
    Iterating over the same processes, types, masses and widths are added to corresponding vectors.
  */
  /*!
    \fn void Phase_Space_Handler::DropRedundantChannels();
    Method to drop channels from the multichannel. The method is to
    generate on phase space ponit and to check it. 
    Then dropping might happen for a number of reasons :
       - a channel produces a zero weight.
       - if two channels produce the same momenta with the same set of
         random numbers /sa Compare
       - for the number of resonating propagators : this is commented out
  */	

  /*!
    \fn   void Phase_Space_Handler::TestIntegration(double);
       A cheap and simple test integration of channels.
  */
  /*!
    \fn bool Phase_Space_Handler::Compare(AMATOOLS::Vec4D *,AMATOOLS::Vec4D *);

    Compares whether two sets of momenta are identical (up to permutations).
    If so a "true" is returned.
  */
  /*!
    \fn bool Phase_Space_Handler::MakeIncoming(AMATOOLS::Vec4D *); 

    Initializes the initial state part of the process:
    - defines default (no ISR) initial state vectors.
  */
  /*!
    \fn      void Phase_Space_Handler::Rotate(AMATOOLS::Vec4D *);

    This method might be used for the usage of the integration with ISR being switched on.
    Let us assume, we have, say, p pbar ->2 jets to be integrated. Obviuosly, there are
    subprocesses like u dbar -> u dbar, and obviously the paton density of the two quarks
    is different for proton and anti-proton. A nice way to have the u-quark stemming from
    both the p and the pbar is to use the same value for the matgrix element but two different
    sets of parton densities (from left and right). This basically boils down to rotating
    the full by 180 degrees to ensure some non-symmetric cuts. In fact that's why we may need
    this method in the future.
  */
  /*
     void Phase_Space_Handler::Init(APHYTOOLS::Flavour *);
     bool Phase_Space_Handler::Check4Momentum(AMATOOLS::Vec4D *);
  */
  /*!
    \fn Phase_Space_Handler::Phase_Space_Handler(AMEGIC::Process_Base *,ISR::ISR_Handler *,BEAM::Beam_Handler *);
    \fn Phase_Space_Handler::Phase_Space_Handler(EXTRAXS::XS_Base *,ISR::ISR_Handler *,BEAM::Beam_Handler *);
    \fn Phase_Space_Handler::~Phase_Space_Handler();

    Explicit constructor for the PS_Handler of a Single_Process or
    a process group. In the latter case, all outgoing flavours are
    Flavour(kf::none) :
       - Reads phase space data : ycut, error, jetfinder, integration method.
       - Initializes a set of four-vectors.
       - Creates a specific, process-dependent selector 
       - Initializes the initial state, incoming vectors, etc..
       - makes a safety copy of the flavours.
       /sa CreateSelector(), InitInitialState()
  */
  /*!
    \fn bool   Phase_Space_Handler::CreateChannelLibrary(std::string,std::string);

    Creates the channels of the multichannel with help of the
    Phase_Space_Generator. The latter is instantiated in this method.
    Then its method Construct is called for all channels. This method checks whether
    a corresponding file in rthe directory "ptype" with characteristic "pID" is
    available. If not, it is created and Construct returns a "true".
    If no new channels were produced, CreateChannelLibrary returns a "true", otherwise
    it returns a "false".
  */
  /*!
    \fn void   Phase_Space_Handler::CollectChannels();

    The following method is used when setting up the integrator of a Process_Group.
    It adds all channels of all processes that have no partner via the method AddChannels that
    is fully recursive.
  */

  /*!
    \fn bool   Phase_Space_Handler::CreateIntegrators();

    Initializes the specific integrator to be used in this integration
    of a Single_Process. This can be a single channel like Rambo or Sarge, 
    or a combination of both, or any combination of the explicit multichannel 
    with or without them. If there are two incoming particles, and if
    isr is to be taken into account, according ISR channels are to be created.
  */

//     void   Phase_Space_Handler::WriteOut(std::string);
//     bool   Phase_Space_Handler::ReadIn(std::string);

  /*!
    \fn double Phase_Space_Handler::Integrate();
    
    Here the actual integration takes place, essentially only
    the Phase_Space_Integrator is called after setting the ISR-range.
  */

  /*!
    \fn double Phase_Space_Handler::Differential();
    A wrapper for the other Differential method when used for actual integrations.
  */

  /*!
    \fn double Phase_Space_Handler::Differential(AMEGIC::Process_Base *);
    The actual value of the cross section is calculated here :
      - First, the incoming beam is set up (via an integration over sprimeB and yB,
        performed by the beamchannels.
      - Following this the initial state characteristics (sprimeI and yI) are determined with 
        the isr integrator, following that suitable initial state vectors are constructed.
       	Here, InitInitialState of the ISR_Handler is used, the use is quite limited to
 	processes with particles that have identical masses.
      - Then, the fsr integrator determines the outgoing momenta.
      - Before selection cuts etc. are applied, the full system is boosted into the
        lab system. It should be stressed that boosts might be mandatory in the selector
 	as well. As an example consider the jet-finder, where for the case of e+e- 
 	collisions and the Durham algorithm the relevant frame is always the c.m.
 	frame of the hard process or the case of DIS, where the frame in question should
 	be the photon's Breit frame.
      - If the selector is passed, i.e. if all momenta are inside the valid region,
        a scale is determined which by default is the hard c.m. energy squared.
 	Eventually this is replaced, in any case, initial PDF's are precalculated if needed.
      - Then the vectors are boosted back into the hard c.m. frame and the value of the
        matrix element squared is evaluated. 
      - Finally, the weights of the three integrators are calculated. Their product, multiplied
        with the value of the ME squared - obtained via Differential - is the result for 
	parton1(2) coming from beam1(2). 
      - Essentially, after flipping the momenta, a similar procedure (calculation of PDF's,
        folding in the ME's) is repeated via Differential2 and yields the total result  
 	that is returned to the Phase_Space_Integrator. 
  */
  /*! 
    \fn bool   Phase_Space_Handler::OneEvent(int mode =0);
    
    Generates a one-weighted event for a Single_Process or a (non-atomar)
    Process_Group. The alogrithm is as follows:
      - Select a process within the group (according to the relative xsec).
      - For up to maxtrials trials generate momentuim configurations, including ISR.
      - Check the configuration whether or not it passes the cuts, if not, try again.
      - Generate the differential cross section at this phase space point and the
        corresponding phase space weight, again including ISR-weights.
      - Check if the product of weight and differential cross section is larger
        than the maximum accumulated before during the integration phase.
        If this is the case, set the new maximum accordingly and stop the
        event generation with a success.
      - If the combined value is smaller than the maximum, accept or reject
        w.r.t. this maximum.
  */

//     bool   SameEvent();

  /*!
    \fn double Phase_Space_Handler::WeightedEvent();
       
    Generates weighted events for a Single_Process or a (non-atomar)
    Process_Group. The alogrithm is as the one for one-weighted events
    with the only difference that no rejection according to the maximum
    is performed. Therefore, no loop of trial is needed.
  */

  /*!
    \fn void Phase_Space_Handler::AddPoint(double);
    
    Here, the differential cross sections are added up sub process-wise.
  */

  /*!
    \fn void Phase_Space_Handler::TestPoint(AMATOOLS::Vec4D *);
    
    With help of Rambo, one set of momenta is created that is used
    to check the amplitude for string and gauge invariance.
    /sa Single_Process::Tests
  */


     // Simple access methods.

//     inline AMATOOLS::Vec4D *     Point()        { return p; }
//     inline APHYTOOLS::Flavour *  Flavs()        { return psflavs; }
//     inline Multi_Channel  * BeamIntegrator()          { return beamchannels; }
//     inline Multi_Channel  * ISRIntegrator()           { return isrchannels; }
//     inline Multi_Channel  * FSRIntegrator()           { return fsrchannels; }
//     inline int              NumberOfBeamIntegrators() { return beamchannels->Number(); }
//     inline int              NumberOfISRIntegrators()  { return isrchannels->Number(); }
//     inline int              NumberOfFSRIntegrators()  { return fsrchannels->Number(); }
//     inline Single_Channel * FSRIntegrator(int i)      { return fsrchannels->Channel(i); }

//     void ISRChannels(int i,Channel_Info & ci) 
//     void BeamChannels(int i,Channel_Info & ci) 

} // namespace PHASIC 

#endif
