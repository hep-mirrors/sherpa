#ifndef PHASIC_Channels_Single_Channel_H
#define PHASIC_Channels_Single_Channel_H

#include "ATOOLS/Org/Message.H"
#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Info_Key.H"

namespace PHASIC {

  class Single_Channel {
  protected:

    std::string      name;

    long int n_points, n_contrib;
    double   weight, result, result2;
    double   res1, res2, res3;
    double   mres1, mres2, mres3;
    double   alpha, alpha_save;

    int     nin, nout;
    double *ms;
    double  sprimerange[3], yrange[2], p_ylimits[2], m_weight;
    double  m_Q, sprime, y, y0, ymin, ymax, Eprime, E1, E2;

    int      rannum;
    double  *rans;

  public:

    // constructor
    Single_Channel();
    Single_Channel(Single_Channel *);
    Single_Channel(int,int,const ATOOLS::Flavour *);
    Single_Channel(int,int,const ATOOLS::Flavour *,ATOOLS::Integration_Info * const) {};

    // destructor
    virtual ~Single_Channel();

    virtual void AddPoint(double);
    virtual void Reset(double value);
    void ResetOpt();

    virtual void Reset() {}

    virtual void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    virtual void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *);
    virtual void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *);

    virtual void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,
			       const double *rans,const int mode=0);
    virtual void GenerateWeight(const int mode=0);

    virtual void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);
    virtual void CalculateLimits();

    virtual void SetRange(double * _sprimerange,double * _yrange);
    virtual void GetRange(); 
    virtual void ISRInfo(int &,double &,double &); 
    virtual void ISRInfo(std::vector<int> &ts,
			 std::vector<double> &ms,std::vector<double> &ws) const;
    virtual int  ChNumber(); 
    virtual void SetChNumber(int);

    virtual std::string      Name();     
    virtual std::string      ChID();
    virtual size_t           Dimension() const;

    virtual size_t NChannels() const; 

    virtual int N();
    virtual int ValidN();
    virtual int Nin();    
    virtual int Nout();  

    virtual double Res1();  
    virtual double Res2();  
    virtual double Res3();  
    virtual double Weight();
    virtual double Alpha(); 
    virtual double AlphaSave();
    virtual double Result();   

    virtual void IncrementN();               
    virtual void SetRes1(double r);         
    virtual void SetRes2(double r);         
    virtual void SetRes3(double r);         
    virtual void SetN(long int n);          
    virtual void SetName(std::string name); 
    virtual void SetWeight(double weight);  
    virtual void SetAlpha(double alpha);    
    virtual void SetAlphaSave(double alpha);

    virtual int  OType();

    virtual void Optimize();
    virtual void EndOptimize();
    virtual bool OptimizationFinished()  { return false; }
    virtual void WriteOut(std::string);
    virtual void ReadIn(std::string);

    virtual void MPISync();

    void CopyMPIValues();

    inline void AddMPIVars(const double &r1,
			   const double &r2,
			   const double &r3)
    { mres1+=r1; mres2+=r2; mres3+=r3; }
    inline void SetMPIVars(const double &r1,
			   const double &r2,
			   const double &r3)
    { mres1=r1; mres2=r2; mres3=r3; }

    inline double MRes1() { return mres1; }
    inline double MRes2() { return mres2; }
    inline double MRes3() { return mres3; }

  };// end of class Single_Channel

  /*! 
    \enum chr::type
    \brief declares the mode for setting phase space integration channel ranges
  */
  
  /*!
    This class is the mother class for all channels that
    are constructed depending on the process.
  */
  /*!
    Some basic ingredients for all channels :
    numbers of legs, their flavours and masses squared (not needed for isr)
    and the name of the channel. This name might be given explicitly 
    (for instance for ISR_Channels) or it may be constructed from amplitudes 
    in the Phase_Space/Channel_Generator.
  */
  /*!
    Stuff for the immediate integration 
  */
  /*!
    Specifics for isr and beamchannels
  */
  /*!
    Explicit constructor. 
    It initializes outgoing masses and the vector of random numbers :
    For them the following construction is applied :
    - two rans for the first pair of outgoing particles (orientation)
    - three rans for each new particle, one for the intermediate propagator,
      two for the orientation of the splitting
  */
  /*!
    Adding a value to the Single_Channel. It basically increments n,
    eventuall n_contrib, and the various result, result2, and res{i} in
    appropriate manners.
  */
  /*!
    Resetting the channels completely for a new integration.
  */
  /*!
    Resetting for the optimization involves a bit more. However, the accumulated
    results are kept. So this involves only what was accumulated between two iteration
    steps for the optimzation procedure. See the Multi_Channel for more information.
  */
  /*!
    Point and weight generation for fsr channels.
  */
  /*!
    Point and weight generation for isr channels.
  */

}

#endif




