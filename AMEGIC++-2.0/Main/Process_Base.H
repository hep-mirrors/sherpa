#ifndef Process_Base_H
#define Process_Base_H

#include "Polarisation.H"

#include "Primitive_Analysis.H"
#include "Primitive_Observable_Base.H"

#include "Phase_Space_Handler.H"
#include "Multi_Channel.H"
#include "Integrable_Base.H"

#include "Phase_Space_Generator.H"

#include "Beam_Spectra_Handler.H"
#include "ISR_Handler.H"
#include "Selector.H"
#include "Flavour.H"
#include "Vector.H"

#include "Interaction_Model_Base.H"
#include "Topology.H"

#include "Message.H"
#include <string>

namespace AMEGIC {
  class Process_Base;
  class Single_Process;

  class Process_Base : public PHASIC::Integrable_Base {
  protected:
    std::string                      m_name, m_resdir;
    int                              m_gen_str;

    int                              m_nvec,   m_nin,    m_nout,    m_nex;
    APHYTOOLS::Flavour             * p_fl,   * p_flin, * p_flout, * p_ex_fl;
    AMEGIC::Pol_Info               * p_pl,   * p_plin, * p_plout;
    int                            * p_b;
    AMATOOLS::Vec4D                * p_moms;

    bool                             m_atoms, m_analyse, m_tables;
    long int                         m_n;
    double                           m_totalxs, m_totalerr, m_totalsum, m_totalsumsqr, m_max, m_rfactor;
    double                           m_last, m_lastdxs, m_lastlumi;
    int                              m_kfactorscheme,m_scalescheme;
    int                              m_orderQCD, m_orderEW,m_nstrong,m_neweak;
    double                           m_scale,m_scalefactor;
    double                           m_isrthreshold;


    PHASIC::Phase_Space_Handler    * p_ps; 
    Phase_Space_Generator          * p_psgen;

    BEAM::Beam_Spectra_Handler     * p_beam;
    ISR::ISR_Handler               * p_isr;
    APHYTOOLS::Cut_Data            * p_cuts;
    APHYTOOLS::Selector_Base       * p_sel;
    APHYTOOLS::Primitive_Analysis  * p_analysis;

    AMEGIC::Process_Base           * p_selected;
    
    /*------------------------------------------------------------------------------

      Naming

      ------------------------------------------------------------------------------*/
  public:
    bool           CheckExternalFlavours(int,APHYTOOLS::Flavour *,
					 int,APHYTOOLS::Flavour *);
  protected:
    std::string *  GenerateNames(int,APHYTOOLS::Flavour *,AMEGIC::Pol_Info *,
				 int,APHYTOOLS::Flavour *,Pol_Info *,
				 std::string &,std::string &,std::string &);
    void           Reshuffle(int,APHYTOOLS::Flavour *,AMEGIC::Pol_Info *);
    bool           IsFile(std::string);
    /*------------------------------------------------------------------------------

      Process initialization

      ------------------------------------------------------------------------------*/
  public:
    Process_Base::Process_Base();
    Process_Base(int,int,APHYTOOLS::Flavour *,
		 ISR::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		 int,int,int,int,int,double,Pol_Info * _pl = 0,
		 int _nex = 0,APHYTOOLS::Flavour * _ex_fl = 0);
    virtual ~Process_Base();

    virtual int    InitAmplitude(Interaction_Model_Base *,AMEGIC::Topology *,
				 AMATOOLS::Vec4D *&,std::vector<double> &,
				 std::vector<AMEGIC::Single_Process *> &)            =0;
    virtual void   InitAnalysis(std::vector<APHYTOOLS::Primitive_Observable_Base *>) =0;
    virtual bool   SetUpIntegrator()                                                 =0; 
    virtual void   AddChannels(Process_Base *,PHASIC::Multi_Channel *,
    			       PHASIC::channelinfos &,PHASIC::channelinfos &);
    virtual void   UpdateCuts(double,double);
    virtual void   InitCuts();
    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    virtual void             SetName(std::string);
    virtual void             SetResDir(std::string);
    virtual void             SetAtoms(bool);
    virtual void             SetTables(bool);

    virtual Process_Base  *  Partner()                                             =0;
    virtual Process_Base  *  Selected()                                            =0;
    virtual void             SelectOne()                                           =0;
    virtual void             DeSelect()                                            =0;
    virtual void             Add(Process_Base *)                                   =0;
    virtual bool             Find(std::string,Process_Base *&)                     =0;
    virtual void             WriteOutXSecs(std::ofstream &)                        =0;
    virtual int              Size()                                                =0;
    virtual Process_Base  *  operator[] (int idx)                                  =0;

    virtual void             SetBeam(BEAM::Beam_Spectra_Handler *);
    virtual void             SetISR(ISR::ISR_Handler *);
    virtual void             SetCuts(APHYTOOLS::Cut_Data *);
    virtual void             SetSelector(APHYTOOLS::Selector_Base *);
    virtual void             SetMomenta(AMATOOLS::Vec4D *);

    virtual void             SetNStrong(int);
    virtual void             SetNEWeak(int);

    virtual void             Empty()                                               =0;
    virtual void             SetTotalXS(int)                                       =0;
    virtual void             SetScale(double)                                      =0;
    virtual void             SetISRThreshold(double)                               =0;
    virtual void             SetTotal(double);
    virtual void             SetMax(double);


    /*------------------------------------------------------------------------------

      Calculating total cross sections and single event generation

      ------------------------------------------------------------------------------*/
    virtual bool   CalculateTotalXSec(std::string _resdir=std::string(""))         =0;
    virtual bool   LookUpXSec(double,bool,std::string)                             =0;
    virtual void   RescaleXSec(double);
    virtual bool   PrepareXSecTables()                                             =0;
    virtual void   AddPoint(const double value)                                    =0;
    virtual double Differential(AMATOOLS::Vec4D *)                                 =0;
    virtual double Differential2()                                                 =0;
    virtual double DSigma(AMATOOLS::Vec4D *,bool)                                  =0;
    virtual double DSigma2()                                                       =0;
    virtual double operator()(AMATOOLS::Vec4D *) { return 0.; }
    virtual double Scale(AMATOOLS::Vec4D *);
    virtual double Scale();
    virtual bool   OneEvent()                                                      =0;
    virtual bool   SameEvent()                                                     =0;
    virtual double WeightedEvent()                                                 =0;

    /*------------------------------------------------------------------------------

      Access methods

      ------------------------------------------------------------------------------*/
    virtual int                             Nin();
    virtual int                             Nout();
    virtual int                             Nvec();
    virtual APHYTOOLS::Flavour            * Flavs();
    virtual AMATOOLS::Vec4D               * Momenta();

    virtual int                             NStrong();
    virtual int                             NEWeak();

    virtual std::string                     Name();
    virtual std::string                     ResDir();
    virtual std::string                     LibName();
    virtual bool                            Atoms();
    virtual bool                            Tables();

    virtual int                             NumberOfDiagrams();
    virtual AMEGIC::Point                 * Diagram(int i);
    virtual bool                            IsFreeOfFourVertex(AMEGIC::Point * _p);

    virtual BEAM::Beam_Spectra_Handler    * Beam();
    virtual ISR::ISR_Handler              * ISR();
    virtual APHYTOOLS::Cut_Data           * Cuts();
    virtual APHYTOOLS::Selector_Base      * Selector();
    virtual APHYTOOLS::Primitive_Analysis * Analysis();
    virtual Phase_Space_Generator         * PSGenerator();
    virtual void                            ISRChannels(int,PHASIC::Channel_Info & );
    virtual void                            BeamChannels(int,PHASIC::Channel_Info & );
    virtual int                             NumberOfBeamIntegrators();
    virtual int                             NumberOfISRIntegrators();
    virtual int                             NumberOfFSRIntegrators();
    virtual PHASIC::Multi_Channel         * FSRIntegrator();
    virtual PHASIC::Single_Channel        * FSRIntegrator(int i);
    virtual int                             ISRNumber();
    virtual void                            ISRInfo(int,int &,double &,double &);
    virtual int                             BeamNumber();
    virtual void                            BeamInfo(int,int &,double &,double &);

    virtual double                          Total();
    virtual double                          TotalError();
    virtual double                          Max();
    virtual double                          Last();
    virtual double                          LastXS();
    virtual double                          LastLumi();
    virtual double                          KFactor(double);
    virtual double                          ISRThreshold();
    

    /*------------------------------------------------------------------------------

      Helpers

      ------------------------------------------------------------------------------*/

    virtual void               PrintDifferential() = 0;

  };
}

#endif







/*! 
  The mother class of all processes : Single_Process, Process_Group, All_Processes.
  
  Every class, which is to be integrated by the Phase_Space_Integrator is derived 
  from this class. 
*/
/*!
  Characteristics for processes: name, number, flavours.
  
  Treatment of the flavours for groups :
  Basically, the incoming flavours should be identical or similar enough
  (like for instance partons in a proton) to allow integration "in one go",
  otherwise they should be set to "none".
  The number of outgoing flavours must be identical, they are to be set to
  "none" in most cases.
*/
/*!
  The set of four momenta, used or instance in InitAmplitude and InitLibrary.
  There it is important to have the very same four vectors.
*/
/*!
  This flag tells whether it is a group of groups or single processes
  to be integrated separately (atoms = 1) or wether the groups or
  single processes can be integrated "in one go" (atoms = 0).
*/
/*!
  This flag steers whether tables/histograms are filled "on the flight"
  of integration (i.e. whether distributions of differential cross sections are evaluated) 
  or whether just one number (i.e. the total cross section) is calculated.
  Depending on the method that starts the integration tables = 0 (CalculateTotalXSec)
  or tables = 1 (CalculateDifferentialXSecs).
*/
/*! 
  Cross sections and their errors. 
  There's a double strategy : 
  First of all the relative contribution of individual channels is optimized within 
  the multichannel such that the overall variance is minimized. To do this, individual 
  results, variance for each /sa Single_Channel are to be stored within each Single_Channel. 
  The total result and its error is accumulated in the /sa Phase_Space_Integrator and in 
  the full /sa Multi_Channel as well. It should be noted that this /sa Multi_Channel
  can also be a /sa Single_Channel only, both are derived from the /sa Channel_Base. 
  It should be stressed here that, if ISR is on, there are two Multi_Channels,
  one for th final state momenta and one for the initial state piece of the
  integration.
  
  On the other hand, results have to be distributed over the processes that are
  integrated. This is especially true for /sa Process_Groups. There the total result
  coincides with the one of the integration stored in /sa Phase_Space_Integrator
  and its specific /sa Channel_Base. However, this total result has to be distributed
  over the /sa Single_Process(es), therefore they have to accumulate their own
  result, variance, number of points etc..
*/

    /*!
      The results for the partial cross section at the last phase space point.
      This is important for instance for the integration, if one /sa Single_Process
      can be mapped onto another one. Then the last amplitude squared lastdxs of this 
      partner will be used and folded with the PDF (lastlumi) that depends on the particular 
      initial state. We have   last = lastlumi * lastdxs.
    */
    /*!
      The scale that enters (for instance) the PDF's. It is clearly process-dependent,
      and I think we'll have to work on this.
    */
    /*!
      These are actual cuts that can be imposed on the final state particles.
      At the moment, energies, invariant masses of pairs of momenta and their relative
      angles are enabled.
    */
    /*!
      This is the selector or a set thereof for actual cuts on the phase space.
      In principle selectors can be added, such that Single_Selectors (like the
      Standard_Selectors : ranges for invariant masses, energies, angles or like
      the Jet_Finder) can be combined in a quite arbitrary fashion.
    */ 
    /*!
      The ISR_Base for the actual process. It is set according to the information of
      the file ISR.dat, at the moment located in the subdirectory Amegic/Testrun.
      At the moment two ISR strategies are implemented : /sa No_ISR for fixed energy 
      evaluations and /sa Structure_Function for the usage of parton densities.
    */
    /*!
      The specific /sa Phase_Space_Handler related to the process. For its instantiation
      from a Process_Base it needs the very same Process_Base ("this" in the instantiation)
      and the related /sa ISR_Base. 
      The /sa Process_Base is needed to generate the fsr channel library in case it is
      a /sa Single_Process and of course this link is needed for the actual integration.
      The /sa ISR_Base is needed to produce the isr channel library and, again, for the
      actual integration, since boosts into lab and c.m. frames, the flux of the incoming 
      particles etc. are handled by the actual /sa ISR_Base. 
    */
    /*!
      This is for groups only. The selection mechanism is via
      SelectOne() and DeSelect() respectively. The Process_Base that is selected
      can be accessed via Selected().
    */
    /*!
      This is for Process_Groups only. The partner points to a Single_Process that had,
      for the same set of four momenta, exactly the same ME (or the same integrand) when 
      doing the gauge and string tests. Algorithmically it boils down to the following:
      Before calling InitAmplitude for the Single_Processes, for each group with
      nout outgoing particles nout momenta are generated to be used in the
      gauge and string tests - and later in InitLibrary. The results of the test
      are stored and then checked for potential mapping onto each other. Then, only
      libraries are initialized for those Single_Processes, who have partner=this,
      the default setting.
      This pointer is widely used for two reasons:
      First, it saves time when doing total cross sections for /sa Process_Groups,
      second it minimizes the number of libaries that have to be read in.
    */
    /*!
      The constructor initializes the process in question. In any case, all
      parameters relevant for integration are set to their starting value, basically 0.
      However, it can be envisionend that we will store cross section information in the
      near future and read it in.

      In Single_Process a little bit more is done via the constructor: 
      Names are generated, as process-identifier, for the libraries and for the
      subdirectory where the libaries will be stored. The process-identifier will
      consist of the flavours, that's unambiguous. For the subdirectory
      the name will have the form "P{nin}_{nout}", and for the library files we will
      have {nin}_{nout}_{diagrams}_{helicities}_{identifier=number}. The constructor
      also calculates the norm and symmetry factors due to spins, flavours etc..
    */
    /*!
      InitAmplitude initializes the amplitudes of all (sub-)processes in the group, or,
      respectively, of the Single_Process. 
      In the latter case, the momenta are generated via the Phase_Space_Handler, if they
      were not already generated before. Then, the /sa Helicity-object, the /sa
      Basic_Sfuncs and the /sa String_Handler are instantiated, before the /sa
      Amplitude_Handler is instantiated handing over also the Topology. In case the
      Amplitude_Handler did not generate any Diagrams, a "false" is returned. Otherwise
      the Basic_Sfuncs and the String_Handler are initialized and gauge and string tests
      are performed. If they failed a "false" is returned, otherwise a "true" will be 
      returned after the follwoing step. The test provides a result that is compared with 
      other results stored in the <double>-vector. In case any identical result is found, 
      the partner will be set accordingly. Otherwise both result and a pointer to the 
      Single_Process are added to the respective extending vectors. Finally

      In the case of Process_Groups, Init_Amplitude essentially iterates over all
      constituents. The only tricky point is related to whether the group is atomic or not,
      i.e. whether the constituents are integrated together or not. In the latter case,
      before each call to a constituent the momenta are reset, in order to keep the
      results independent of each other.
    */
    /*!
      Adding processes to the group of processes in Process_Group.
      This is somewhat constrained for Process_Groups. There, the first element, i.e.
      the first Process_Base to be added defines the number of incoming and outgoing
      legs, and, therefore, momenta for the integration. If a Process_Base to be added
      does not coincide in these numbers it is rejected. This handling is different
      for /sa All_Processes.
    */
    /*!
      This is used only by All_Processes. See explanations there.
    */
    /*!
      Setting up the integrators, either directly by creating or reading from the 
      libraries or as a sum of integrators of individual processes.

      SetUpIntegrators creates the integration channel(s) with help
      of the Phase_Space_Handler. The individual channels of the processes
      are added up and then they are checked for double occurrences to be 
      deleted from the common multi-channel.
    */
    /*!    
      Initializes the process dependent Cut_Data. Here, they are just
      initialized, i.e. the corresponding vectors are allocated. T
      They will have to be filled by the selector.

      At the moment, the cuts contain energies, angular regions with respect 
      to the beam, invariant masses, and mutual angles.
      /sa Selector, Jet_Finder, Cut_Data
    */
    /*!
      This method steers the calculation of total cross sections.
      Then the Phase_Space_Handler is employed to integrate the process via
      Integrate() and the result is stored internally via SetTotalXS():
    */
    /*!
      Adding a point, i.e. the contribution of one phase space point to the totalsum
      yielding the total cross section at the end. The version here is the one
      for Single_Processes. For Process_Groups the value will be split in due proportion
      and then added via calling AddPoint of the constituents.
    */
    /*!
      From the accumulated sum, this method generates total cross sections and
      their errors.
    */
    /*!
      Prints out the contribution at one specific point, the differential contribution
      from the matrix element and the parton luminosity.
    */
    /*!
      Differential yields a differential cross section in dependence on the
      four vectors, the phase space point.
      Differential calls DSigma with lookup set to false that will calculate the differential 
      cross section at this point including the ISR weight. In DSigma there's the chance
      to refer to identical differential cross sections, that may or may not
      have identical partonic luminosities, like for instance dd->dd vs. uu->uu 
      in proton-proton collisions. In pure QCD the MEs of these subprocesses yield 
      the same result and are just folded with different PDFs, however, there's no 
      need to calculate the ME twice. In any case, the Phase_Space_Handler calls
      Differential, and, depending on whether this is a non-atomar Process_Group
      or not, DSigma is called with the lookup flag equals 1 (groups) or 0
      (single_processes).
    */
    /*!
      The differential cross section including the isr-weight. In case the boolean 
      "lookup" is true and the reference is not the identical process (partner != this), 
      DSigma will use the differential ME of this very reference (via Lastdxs). 
      Then it will fold it with the corresponding parton luminosity via isr->Weight(flin).
      It should be noted that the parton luminosities at this specific point have been
      calculated beforehand, namely in the method "Differential" of the Phase_Space_Handler,
      and this is also the place, where the scale an the energy-fractions of the partons enter.
      At this stage it is only their weight that is extracted, since it depends on
      specific partons.
    */
    /*!
      This method returns the value of the amplitude squared.
    */
    /*!
      This method partially kills ingredients that are not needed any longer.
    */
    /*!
      Scale sets the scale for the isr-weight and eventually, the couplings.
      In the default version, the scale is the c.m. energy squared. For Process_Groups,
      the scale is propagated to the lower levels in the hierarchy via SetScale.      
    */
    /*!
      These methods yield single weighted or unweighted events.
      For Single_Processes the strategy is as follows:
      The Phase_Space_Handler is connected with the suitable integrator set up in the 
      calculation of the total cross section for the process. Then repeatedly, up to 
      _maxtrials_ times, the PS_Handler calls the integrator to produce
      a point in phase space. This might involve initial state radiation as well.
      If the point passes all cuts handled by the selector, a phase space weight
      and the ME are evaluated. The point is taken or rejected according to a 
      hit-or-miss rejection with the corresponding maximal value of the subprocess.

      For Process_Groups this is a little bit mroe tricky and depends on whether the
      group is atomar or not. In any case, first one subprocess is selected.
      If the group is atomar, the according single event of the selected atom is
      returned. If the group is non-atomar, the Phase_Space_Handler again is called,
      calling the integrator etc. . The only difference is that in contrast to above
      the ME and the maximum of the selected atom is used.

      I'll have to work on this again !!!!




    virtual int     NumberOfDiagrams() { 
      AORGTOOLS::msg.Error()<<"Virtual Method : Process_Base::NumberOfDiagrams()"<<std::endl; 
    }
    virtual Point * Diagram(int i) { 
      AORGTOOLS::msg.Error()<<"Virtual Method : Process_Base::Diagram()"<<std::endl; 
    }
    virtual bool    IsFreeOfFourVertex(Point * _p) { 
      AORGTOOLS::msg.Error()<<"Virtual Method : Process_Base::IsFreeOfFourVertex()"<<std::endl; 
    }
    */
