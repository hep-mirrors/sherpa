#ifndef Multi_Channel_h
#define Multi_Channel_h

#include "Info_Key.H"

namespace PHASIC {

  class Single_Channel;
  class PI_Interface;

  struct Channel_Info {
    int                  type;
    std::vector <double> parameters;

    bool operator==(const Channel_Info & ci) {
      if (type != ci.type)                           return 0;
      if (parameters.size() != ci.parameters.size()) return 0;
      for (unsigned int i=0;i<parameters.size();i++) {
	if (parameters[i] != ci.parameters[i])       return 0;
      }
      return 1;
    }
    bool operator!=(const Channel_Info & ci) { return (!(operator==(ci))); }
  };

  typedef std::vector<PHASIC::Channel_Info> channelinfos;

  class Multi_Channel {
  protected :
    int                           nin, nout;
    ATOOLS::Flavour             * fl;
    std::string                   name;
    long int                      n_points,n_contrib;
    double                        m_weight;
    double                        rans[2];
    double                        s1xmin, * s1, * s2;
    std::vector<Single_Channel *> channels, m_best;
    bool m_readin, m_fixalpha;
    int m_lastdice,m_optcnt;
    double m_pol;
  public:
    Multi_Channel(std::string);
    ~Multi_Channel();
    void             Add(Single_Channel *);  
    void             DropAllChannels();    
    void             DropChannel(int);     
    Single_Channel * Channel(int); 
    void    Reset();
    void    ResetOpt();
    void    Optimize(double);
    void    MPIOptimize(double);

    /*!
      Endoptimize replaces, after an ordinary optimization step, the set of a-priori weights
      with the best set so far which has been stored during the optimization procedure.
      Channels, whose a-priori weight drops under a critical value are discarded
      in the follwoing by setting their a-priori weights to 0.
    */
    std::string SelectChannel();

    void    EndOptimize(double);
    void    AddPoint(double);
    void    GeneratePoint(int,ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *);
    void    GenerateWeight(int,ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    void    GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    void    GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,PI_Interface*);
    void    GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);     

    void    GeneratePoint(int,ATOOLS::Vec4D *,double *);
    void    GenerateWeight(int,ATOOLS::Vec4D *);
    void    GeneratePoint(ATOOLS::Vec4D *);
    void    GenerateWeight(ATOOLS::Vec4D *);     

    void    GeneratePoint(int n,double & sprime,double & y,int mode);
    void    GenerateWeight(int n,double sprime,double y,int mode);
    void    GeneratePoint(double & sprime,double & y,int mode);
    void    GenerateWeight(double sprime,double y,int mode);

    void    GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,int mode);
    void    GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,int mode,
			  PI_Interface *pi);
    void    GenerateWeight(int mode);

    void    ISRInfo(int,int &,double &,double &);

    inline  Single_Channel *operator[](size_t i) { return channels[i]; }

    void    Print();
    void    SetRange(double * sprimerange,double * yrange); 
    void    GetRange(); 
    void    WriteOut(std::string);
    bool    ReadIn(std::string);

    inline void FixAlpha() { m_fixalpha=true; }

    inline int                  Number()    { return channels.size(); }
    inline int                  N()         { return n_points; }
    inline int                  ValidN()    { return n_contrib; }
    inline int                  Nin()       { return nin; }
    inline int                  Nout()      { return nout; }
    inline ATOOLS::Flavour     *Flavs()     { return fl; }
    inline double               Weight()    { return m_weight; }
    inline std::string          Name()      { return name; }
    std::string                 ChID(int);
    Single_Channel*             LastChannel() { if (m_lastdice>=0) return Channel(m_lastdice); return 0;}
    void                        NoDice()    { m_lastdice=-1; }
    inline const std::vector<Single_Channel*> &Best() const { return m_best; }

  };


  /*!
    This is the Multi_Channel, AMEGICs preferred integration
    tool. The idea is to have a vector of Single_Channels
    that are repeatedly called according to their a priori 
    weights (alpha in the channels). These weights are optimized
    during integration such that the overall variance is minimized.
  */
  /*!
    Some basic ingredients for all channels :
    numbers of legs and their flavours (not needed for isr)
    and the name of the channel. This name might be given explicitly 
    (for instance for ISR_Channels) or it may be constructed from amplitudes 
    in the Phase_Space/Channel_Generator.
  */
  /*!
    Stuff for the immediate integration 
  */
  /*!
    Managing the multichannel, adding and dropping channels as well as
    access to individual channels.
  */
  /*!
    Reset does a complete reset of the multi-channel. All a-priori weights are reset
    to uniform probabilities, the individual channels are reset, and the actual weights
    of all channels are set to 0. Furthermore, the result and the result2 within
    the multi-channel are set to 0, the minimal spread so far is set to some 
    arbitrary high value. If not existing so far, spread-vectors are initialized.      
  */
  /*!
    Here the counters within an iterations are reset to 0. 
    I doubt whether we need this.
  */
  /*!
    Optimzing the a-priori weights of the multi-channel.
    The idea is to compare the individual variance with their average and
    reshuffle the a-priori weights such that all variances would equal the
    average. What I find strange in this method at the moment is the occurence
    of s2 ... . I'll have to check this ... .
  */
  /*!
    Endoptimize replaces, after an ordinary optimization step, the set of a-priori weights
    with the best set so far which has been stored during the optimization procedure.
    Channels, whose a-priori weight drops under a critical value are discarded
    in the follwoing by setting their a-priori weights to 0.
  */
  /*!
    A value is added to the multi-channel, and, in due proportion, to the
    individual channels. This amounts to incrementing the results and the results squared
    of the channels, for the multi-channel, the global result(2) is incremented,
    for the individual channels within, the res(1,2,3) are incremented that are valid
    only during one iteration.
  */
  /*!
    The variance of a statistical distribution (values accumulated in result and their
    squares in result2) over a number of trials. This is equivalent to one Standard deviation
    and thus a measure for the quality of an integral-estimate.
  */
  /*!
    Methods for the generation of phase space points
  */
  /*!
    This is used for ISR Channels to set the y and s' range.
  */
  /*!
    Simple access methods
  */
}


#endif
