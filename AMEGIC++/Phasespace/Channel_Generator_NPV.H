#ifndef AMEGIC_Phasespace_Channel_Generator_NPV_H
#define AMEGIC_Phasespace_Channel_Generator_NPV_H

#include "AMEGIC++/Phasespace/Channel_Generator_Base.H"


namespace AMEGIC {

typedef std::map<std::string,std::string> Decls;
  
  typedef std::vector<std::string> String_List;

  

  class Channel_Generator_NPV : public Channel_Generator_Base {
    Decls declarations;
    
    int         extrachannelflag,newchannel;
    int         tcount,acount,m_aid;
    String_List m_idc;
    std::string m_idstr,m_mapstr;
    std::vector<Point *> m_topos,m_pclist;
    void        Step0(int,Point *,int&,std::ofstream&);
    bool        QCDAntenna(int,Point*,int&,std::ofstream&,int);

    void        GenerateMassChain(int,Point *,Point *,int&,std::ofstream&);
    void        GenerateMassFwd(int,Point*,int&,std::ofstream&);
    void        GenerateDecayChain(int flag,Point* p,int& rannum,std::ofstream& sf,
				   std::vector<std::string>, std::vector<std::string>);
    std::string LinkedMasses(Point *);
    void        SetProps(Point *,Point **,Point **,int&);
    void        CalcSmin(int,const char*,std::string,std::ofstream&,Point*);
    void        CalcTSmin(int,String_List&,std::ofstream&);
    void        IdentifyProps(Point *);
    void        BackLinks(Point *,Point * &);
    void        InitT(Point *);
    bool        Massive(ATOOLS::Flavour fl) { return (fl.Mass()!=0.) ? 1:0; }
    int         AntennaS(Point* p);
    void        GenerateTopos(Point*);
    Point*      CopyPoints(Point*);
    Point*      TransformTS(Point*);
    Point*      GetMirrorTopo(Point*);
    int         MarkNP(Point*);
    void        MRPScan();
    std::string Order(std::string);
    std::string IString(int);
    void        AddToVariables(int,const std::string&,const std::string&,const int&,std::ofstream&); 
    bool        CheckVariables(int flag,const std::string& lhs,const int& type);
    void        ClearDeclarations() { declarations.clear(); } 
    double      PMassSum(Point*,std::vector<int>*);
    std::string GetFlMass(Point* p);
  public:
    Channel_Generator_NPV(int,int,Point *,int);
    ~Channel_Generator_NPV();
    int         MakeChannel(int&,int,std::string&,std::string&);
    Point     * GetPointlist()             { return plist; }
    //void        SetName(std::string _name) { name = _name; }
    std::string CreateChannelID(int);
    int         NumberOfChannels() { return m_topos.size(); }
  };
      //! The usual identifiers, number of incoming and outgoing legs.
    /*!
      This is a point list, i.e. a copy of the point list of the amplitude for
      which this channel is constructed. The copy is done with help of the class
      topology.
    */
    /*!
      The number of t-channel propagators helps in selecting the correct "principal"
      topology. This number is deterimined via the method IdentifyProps.
      Actually this might help also for efficient ISR channels in non-s channel
      cases later on.
    */
    /*!
      Step0 basically plays the role of distributing the different basic topologies for the final 
      state integral. Depending on tcount, the number of t-channel propagators in the amplitude,
      StepS (0 t's) or StepNT is called, for both the "Momenta"- and the "Weight"-mode. For StepS 
      some initialization work is done in Step0, mainly bacause StepS is called by the other 
      channels as well and because the invariant mass of the propagator in the s-channel has to 
      be known then.
    */

    /*!
      StepS initializes - if need be - the two invariant masses of the decay products by suitable 
      propagator terms via GenerateMasses. Then it decides on the decay type, either
      anisotropic or isotropic, depending on the flavour constellation.
      Finally it calls StepS for its decay products.
    */
    /*!
      StepNT is the start of a recursive method to fill multiple T-channels. Assume you have a 
      _ _ _    "comb-like structure", StepNT treats it by deciding which point is which type, 
      | | | |  s- or t-channel via SetProps. Via GenerateMasses the s-channel props receive 
      | | | |  their masses. However, StepNT then defines four vectors of legs, two incoming 
               and two outgoing each. Starting from the edges of the comb, either the left- 
      or the right-most s-channel prop and all the other s-channel props are treated as two 
      outgoing partner for a single t-channel. The single s-channel will be treated as belonging 
      to an incoming leg to continue to distribute the bulk of all the other s-channel 
      props in the next recursion step.
      This recursion, shuffling one s-channel leg into one incoming leg per step, is doen by 
      SingleTStep.
    */
    /*!
      SingleTStep initializes phase space boundaries and calls then a propagator structure 
      peaking at smin for the bulk of s-channel particles. Having two outgoing masses, the 
      TChannel method - Momenta or Weight - is called. Then a new s-channel prop to be singled
      out is selected.
    */
    /*!
      For a number of points, GenerateMasses defines the minimal invariant masses. It then works 
      consecutively and tries to find a sequence of propagators to be set such that the more 
      resonating a propagator is the earlier its mass will be selected to maximise the chance
      of it contributing significantly.
    */
    /*!
      Linked masses proceeds iterativley from a given point and creates a string that consists of 
      all numbers of outgoing particles that are off-springs of this point.
    */
    /*!
      SetProps selects t- and s-channel propagators in comb-like structures.
    */
    /*!
      CalcSmin determines the minimal invariant mass for a single propagator, defined by the sum 
      of the minimal invariant masses of all pairs of offsprings.
    */
    /*!
      CalcTSmin calculates the minimal invariant mass for a bulk of propagators 
      that has been already equipped with invariant masses.
    */
    /*!
      IdentifyProps checks for t-channel propagators.
      Starting by setting all t-flags in the points to zero via InitT, it first initializes the 
      previous-links in the points via BackLinks until it finds the second incoming endpoint 
      (with b=-1). Starting from there it walks over the previos links setting all
      t-flags on the way and incrementing the t-channel prop-counter step by step.
    */
    /*!
      A simple check whether a specific flavour is massive or not.
    */

    /*
      void    SetDouble(std::string,double);
      double  GetDouble(std::string);
      void    SetVector(std::string,Vec4D);
      Vec4D   GetVector(std::string);
      void    Clean();
      void    PrintSettings();
      ostream Translate(ATOOLS::Flavour);
    */
    /*!
      Copies the pointlist to make sure that we can work on it, identifies the types of 
      intermediate lines - whether they are s- or t-channel.
    */
    /*!
      This is the method seen from outside to initalize the channel in the library. It basically 
      produces the C-file of the channel and the methods of the file are initialized. MakeChannel 
      is called with a counters for the amplitude. 
      The individual methods, basically Weight and Momenta are filled by calling the method Step0 
      with different flags to switch in their respective "Momenta"- or "Weight"-mode. A counter 
      max_numb is passed through all these methods to count and collect the resonant propagators and 
      their flavours. Havong done that, with flag = 2, Step0 is passed again to fill in the
      ISRtype method. There, if tcount = 0, the isrtype will be set to 1, a Breit-Wigner pole
      with corresponding mass and width (that may be zero, such cases will be filtered in
      the MakeISRChannels). For t-channel type structures I'd like to have a similar structure,
      this is to be done. Finally MakeChannel returns the number of random numbers needed.
    */
}  
#endif

