#ifndef Single_Process_H
#define Single_Process_H

#include "Process_Base.H"

#include "Polarisation.H"
#include "Helicity.H"
#include "Amplitude_Handler.H"
#include "Standard_Selector.H"
#include "Spin_Correlation_Tensor.H"
#include <vector>

namespace AMEGIC {
  class Single_Process : public Process_Base {
  private:
    std::string             m_ptypename,m_libname,m_pslibname;
    Polarisation            m_pol;
    double                  m_Norm;
    double                  m_save_max;
    bool                    m_newlib;
    double                  m_iresult;
    int                     m_libnumb;
    double                  m_sfactor;

    Helicity              * p_hel;
    Basic_Sfuncs          * p_BS;
    Amplitude_Handler     * p_ampl;
    String_Handler        * p_shand;

    Single_Process        * p_partner;
    /*------------------------------------------------------------------------------

      Constructors

      ------------------------------------------------------------------------------*/
  public:

    Single_Process(int,int,ATOOLS::Flavour *,
		   PDF::ISR_Handler * =NULL,BEAM::Beam_Spectra_Handler * =NULL,
		   ATOOLS::Selector_Data * =NULL,
		   int=1,int=99,int=99,int=0,int=0,
		   double=-1.,Pol_Info* _pl=0,int _nex=0,ATOOLS::Flavour * _ex_fl=0,
		   int usepi=0,double ycut=-1.,double error=-1.,std::string e_func="1");
    Single_Process(Process_Info*,int,int,ATOOLS::Flavour *,
		   PDF::ISR_Handler *,BEAM::Beam_Spectra_Handler *,
		   ATOOLS::Selector_Data *,
		   int,int,int,int,int,
		   double,Pol_Info* _pl,int _nex,ATOOLS::Flavour * _ex_fl,
		   int usepi,double ycut,double error,std::string e_func="1");
    ~Single_Process();
    void Add(Process_Base *) {
      ATOOLS::msg.Error()<<"Error in Single_Process. Method Add called."<<std::endl;
      abort();
    }

    /*------------------------------------------------------------------------------

      Generic stuff for initialization of Single_Processes

      ------------------------------------------------------------------------------*/
  private:
    void           PolarizationNorm();
    double         SymmetryFactors();
    double         SBSymmetryFactor(ATOOLS::Flavour* flout,size_t nout);
    void           PrepareTerminate();
    bool           CheckMapping(const Process_Base * proc);

    /*------------------------------------------------------------------------------

      Initializing libraries, amplitudes, etc.

      ------------------------------------------------------------------------------*/
  protected :
    int                 Tests();
    int                 CheckLibraries();
    int                 CheckStrings(Single_Process*);
    void                WriteLibrary();
    std::string         CreateLibName();
    void                CreateMappingFile();
    bool                CreateChannelLibrary();
  public:
    bool                NewLibs() {return m_newlib;}
    bool                FoundMappingFile(std::string &,std::string &);
    int                 InitAmplitude(Interaction_Model_Base *,Topology *,ATOOLS::Vec4D *&,
				      std::vector<Process_Base *> &,std::vector<Process_Base *> &,
				      int &,int &,int &);
    int                 InitAmplitude(Interaction_Model_Base *,Topology *);
    void                InitDecay(Topology* top); 
    bool                SetUpIntegrator();
    bool                IsFreeOfFourVertex(Point * _p) { return 1-p_ampl->ExistFourVertex(_p); }
    String_Handler    * GetStringHandler()             { return p_shand;}
    Amplitude_Handler * GetAmplitudeHandler()          { return p_ampl;}
    Helicity *          GetHelicity()                  { return p_hel; }    
    double              Result()                       { return m_iresult; } 

    /*------------------------------------------------------------------------------

      Process management

      ------------------------------------------------------------------------------*/
  public:
    bool             SelectOne();
    bool             SelectOneFromList();
    void             DeSelect();
    bool             ReSelect(int);
    size_t           Size();
    bool             Find(std::string,Process_Base *&);
    void             WriteOutXSecs(std::ofstream & _to);
    Process_Base   * operator[] (int idx);
    std::string      LibName()                          { return m_libname;     }
    std::string      PSLibName()                        { return m_pslibname;   }
    void             SetTotal(int flag, int depth=0);
    void             FixISRThreshold();
    void             Empty();
    void             SetPartner(Single_Process * _p)    { p_partner = _p;       }
    Process_Base   * Partner()                          { return p_partner;     }
    void             Minimize();

    /*------------------------------------------------------------------------------

      Calculating total cross sections

      ------------------------------------------------------------------------------*/
  public:
    bool           CalculateTotalXSec(std::string _resdir=std::string(""));
    bool           LookUpXSec(double,bool,std::string);
    bool           PrepareXSecTables();
    void           AddPoint(const double value); 
    double         Differential(const ATOOLS::Vec4D *);
    double         Differential2();
    double         DSigma(const ATOOLS::Vec4D *,bool);
    double         DSigma2();
    double         operator()(const ATOOLS::Vec4D *);
    //! Return the spin correlations for the current process.
    virtual ATOOLS::Spin_Correlation_Tensor* GetSpinCorrelations();
    ATOOLS::Blob_Data_Base * OneEvent(double = -1.);
    ATOOLS::Blob_Data_Base * SameEvent();
    ATOOLS::Blob_Data_Base * WeightedEvent(const int mode=0);
    ATOOLS::Blob_Data_Base * WeightedEventNS(const int mode=0);
    ATOOLS::Blob_Data_Base * SameWeightedEvent();
    int            NumberOfDiagrams();
    Point        * Diagram(int i);
    void           ResetMax(int);
    void           OptimizeResult();

    /*------------------------------------------------------------------------------

      Helpers
      
      ------------------------------------------------------------------------------*/
  public:
    void           PrintDifferential();
    /*------------------------------------------------------------------------------

      Stuff to sample over helicities
      
      ------------------------------------------------------------------------------*/
  private:
    bool                    m_helsample,m_inithelsample;
    std::vector<int>        m_helnumbers;
    std::vector<double>     m_helalphas, m_helalphasaves;
    std::vector<double>     m_helresults, m_helresults2;
    std::vector<long int>   m_helthrows;
    long int                m_throws;
    int                     m_helnumber;
    double                  m_helresult, m_helresult2;

    void           InitializeHelicityWeights();
    const int      SelectedHelicity();
    void           AddToHelicity(const double,const int);
    void           OptimizeHelicityWeights();
    double         operator()(const ATOOLS::Vec4D *,const int);
  public:
    void           SetHelicitySampling(bool _helsample) { m_helsample = _helsample; }
  };
}



#endif
/*! 
  \class Single_Process
  \brief Incorporates the handling of single processes
 



  \fn int Single_Process::InitAmplitude(Topology *,Vec4D *&,
                                        vector<double> &,vector<Process_Base *> &)
  \brief Constructs the amplitudes and deals with the libraries

  The return values signify :
  <ul>
  <li> -2  Error in amplitude creation. Tests went wrong.
  <li> -1  No amplitudes created. Process should be deleted.
  <li>  0  New library created, has to be compiled and linked. This is a
           potential outcome of the call to InitLibrary.
  <li>  1  No new library created. For this, there are two potential reasons :
           <ol> 
	   <li> It was possible to map the actual amplitude onto another one,
	        without giving any hint whether this has produced a new library
		file. In this case the libname needed for the creation of
		the mapping file through CreateMappingFile is set to the
		partners library name.
           <li> The library named in the .map file was already compiled and could
	        be found, or the process could be mapped onto another library
                that already existed. Both outcomes are possible via InitLibrary.
	   </ol>
  </ul>	

  The tasks in this method in their sequence are :
  <ul>
  <li> Momenta to perform the tests are either handed over or newly created.
  <li> Helicity-, Basic_Sfuncs- and the String_Handler object are instantiated.
       They are needed to then initialize the Amplitude_Handler.
  <li> The construction of the Amplitude_Handler already creates the diagrams.
       It is checked whether any diagrams were calculated. If not, the process
       should be deleted from the list.
       /todo I'll have to do & check this carefully.
  <li> The Polarisation-, Basic_Sfuncs- and the String_Handler are initialised
       such that the Tests of the amplitude can be performed. 
  <li> The Tests produce a result for the amplitude value at the phase space
       point in question, i.e. for the set of momenta. This result is compared
       with a list of results handed over into that method stemming from other
       Singel_Processes to find eventual mapping options. 
       If such an option exists, the partner is set and the .map file is
       filled accordingly. The method is left with return value 1.
       If there is no mapping option, the result is added to the vector of results
       and the pointer of this Single_Process is added to the list of corresponding
       pointers to allow for checks of subsequently initialised processes. Then 
       InitLibrary is called, checking for libraries yielding the same result
       and eventually writing out the library file.
  </ul>



  \fn int Tests(double&)
  The return values signify :
  <ul>
  <li> 2 A .map file has been found. The result of the strings stored
         in this file coincide with the one just calculated. Then the String_Handler is
	 updated and one can skip any further test.
  <li> 1 All tests have been performed and work out. No suitable .map file yielding
         the same result has been found. 
  <li> 0 Either the string test or the gauge test misfired.
  </ul>
  Tests of the amplitude : Gauge and string tests. In doing so, a test result for
  the given set of momenta is obtained that will be used in the libary initialisation,
  InitLibrary. 


  \var double Norm
  The norm of the individual process. It consists of :
  <ul>
  <li> a symmetry piece for the particles calculated via /sa SymmetryFactors(),
  <li> a spin and colour piece determined via SpinAverage(nin,flin) in Polarisation
  <li> a piece 3/(2 m^2) for each massive vector boson of mass m stemming from the
       treatment of the polarisation vectors.
  </ul>

  \var string ptypename 
  It is used to determine into which directory libraries are stored, its form 
  is P{nin}_{nout}.

  \var string libname 
  The library name is the basic name (up to an identiying number at the end) of all
  library files, it is {nin}_{nout}_{number of helicities}_{number of diagrams}

  \var Polarisation * pol
  The Polarisation class deals with all global aspects of external polarisations, in principle
  \todo I would like to have some treatment of polarized cross sections there ... 

  \var Helicity * hel
  The Helicity class is used for the actual sum over helicities. 

  \var Amplitude_Handler * ampl
  It actually manages the generation of all single amplitudes. This process is started and 
  finished when the constructor is called. The Amplitude_Handler also manages the evaluation 
  of the matrix element squared through the two methods differential - one for use with 
  strings and one without. For the former methods are provided that proagate amplitude-related 
  information like for instance couplings into the strings, i.e. into the String_Handler.

  \var Basic_Sfuncs * BS
  They form the core of the calculation. From them, all complicated structures are deduced. 
  The Basic_Sfuncs correspond to the bracket symbols [ij] and <ij> notorious in the helicity 
  formalism.

  \var String_Handler * shand
  It is responsible for creating the strings to be stored in the library files., The strings
  are also used for quicker evaluation.

*/
/*!
  Initialize basic classes, mainly the Phase_Space_Handler. It will be needed
  to generate a test point and cross check with existing libaries.
  All other stuff will be set to zero.
*/
/*!
  Calculates the norm of the process due to a norm stemming from the polarizations
  and from the symmetry factors of the process due to the flavours.
*/
/*!
  The dummy version of creating a void selector with no effect whatsoever.
*/
/*!
  A little helper, just checking whether a specific file exists at all.
  This method is also available elsewhere, we might move it to the ATOOLS.
*/ 
/*!
  When being called, the constructor copies the flavours into flin and flout, thus
  filling in these vectors. Then, names are generated, for the process, the subdirectory
  and the naming root of the library. Finally, some global factors are initialized via
  PolarizationNorm and the cuts and the Phase_Space_Handler are initialized.
*/
/*!
  InitLibrary is called from All_Processes with the result of the Tests evaluated before.
  If InitLibrary is  called, a cross-reference was not possible (or atoms was set, but this
  would be quite stupid in most cases). The first task within InitLibrary is to extend
  the library name with the number of graphs and the helicities.
  
  Comment on the description :
  It would be nice to add the number of Zfunctions and couplings as well, just to be more 
  selective later on. Actually this was a good idea, because this is what killed me in this method.
  Maybe this is what you, Ralf, could do. Another idea I had is that sometimes we have
  exactly the same amplitudes, but different couplings. An example is e+e- -> ddbar vs.
  e+e- -> uubar. It would be real cool to check whether different couplings alone would
  do the trick. This would reduce the number of libraries to be linked even further.
  Another item on my private shopping list would be to set the couplings in dependence
  on kinematics, an inclusion of running alpha_s would be extremely nice, especially for
  hadronic initial states !

  However, the next step is to instantiate a second String_Handler in this method for
  testing purposes with the String_Generator as an argument. Then testnames are
  constructed, i.e. names for libraries with the same parameters but different identifier,
  an incrementing number : {nin}_{nout}_{nhelicities}_{ngraphs}_{number}
  Then the dummy - String_Handler looks for libraries and tries to evaluate their amplitudes
  at the same set of vectors that were actually stored before hand. If results are identical,
  the originaL String_Handler will load in the corresponding library and this method
  exits with return value "true". If no match was found, a library name is generated by
  taking the parameters and incrementing the identifier number until no file was found.
  This might look a little bit awkward, but it is meant for the situation when we talk
  about many similar processes. Then it is possible that other library files were written out
  in the same run but are different to the Single_process under consideration.
  However, having generated the name the library is written out and "false" is returned,
  signifying that some more compilation has to be performed. 
*/
/*!
  This is to set up the integrator. First Multi_Channels for the initial and the final state 
  momenta are initialized and the corresponding final state channel library is tested via 
  CreateChannelLibrary of the Phase_Space_Handler. If the channels are already there and linked 
  CreateIntegrator of Phase_Space_Handler will stick them together, test them and eventually erase 
  some channels. In fact, CreateIntegrator is also responsible for the initial state part.
  In addition to the steps outlined above for the final state part, CreateIntegrator will
  iterate over the final state channels and it will fill the three vectors isr_types, isr_masses, 
  and isr_widths via suitable methods within the Phase_Space_Handler. That way initial state 
  channels are added "hard wired" into another multi channel integrator for the
  initial state part of the integral.
*/
/*!
  CalculateTotalXSec does three things: It sets up the selector that has to be specified
  (for nice examples check some predefined process_groups like, e.g. Drell-Yan or 
  Four_Fermions), then it employs the Phase_Space_Integrator via the Phase_Space_Handler
  to integrate the amplitude over the phase space, and finally it sets the total cross
  section.
*/
/*!
  Differential essentially calls DSigma with lookup = "false"
*/
/*!
  DSigma is the "true" evaluation of the cross section. Depending on the lookup and
  on the existence of a "partner"-process, the result of either the partner-process or
  the operator are evaluated and, eventually, folded with the pdf of the initial state.
*/
/*!
  This method returns the value of the amplitude squared.
*/
/*!
  Here I'll have to work & test !!!
*/
/*!
  InitLibraries initializes the amplitude libraries. It checks all other libaries
  eventually existent in the subdirectory whether they produce the same value for the given
  set of vectors. In fact, before that partner processes are established for 
  Process_Groups to make sure that one look-up is sufficient for all Single_Processes
  whose partner is "this". In case such a library is found the according libary is
  loaded into the String_Handler otherwise a new library is saved. The saving proceeds
  in two steps :
  - First, the output file is set for the String_Handler via SetOutput(filename),
  - then, the full library is created via ampl->Output(filename). This includes
  as well colour factors etc., therefore the output cannot proceed via the
  String_Handler alone. 
  In case a libary or a partner was existent, "true" is returned, otherwise, i.e. if a 
  new libary was saved, a "false" is returned.
  In this method there is one of the remaining problems : I have to initialize
  a new String_Handler for the comparison with a copy of the String_Generator.
*/

