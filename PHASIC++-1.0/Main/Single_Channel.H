#ifndef Single_Channel_h
#define Single_Channel_h

#include "Message.H"
#include "Cut_Data.H"
#include "Flavour.H"
#include "Vector.H"
#include "Info_Key.H"

namespace PHASIC {

  /*!
    This class is the mother class for all channels that
    are constructed depending on the process.
  */
  class Single_Channel {
  protected:
    int                  nin, nout;
    ATOOLS::Flavour * fl;
    double             * ms;
    std::string          name;
    long int             n_points,n_contrib;
    double               weight,result,result2;
    double               res1,res2,res3;
    double               alpha,alpha_save;
    int                  rannum;
    double             * rans;
    double               sprimerange[3], yrange[2];
    double               sprime,y,y0,ymin,ymax,Eprime,E1,E2;
  public:
    Single_Channel():ms(0),rans(0) { };
    Single_Channel(Single_Channel *);
    Single_Channel(int,int,ATOOLS::Flavour *);
    virtual ~Single_Channel();
    void AddPoint(double);
    void Reset(double value);
    void ResetOpt();
    virtual void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const double *rans,const int mode=0);
    virtual void GeneratePoint(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey,const int mode=0);
    virtual void GenerateWeight(const int mode=0);
    
    virtual void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);
    virtual void GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    virtual void GeneratePoint(ATOOLS::Vec4D *,ATOOLS::Cut_Data *,double *); 
    virtual void GenerateWeight(ATOOLS::Vec4D *,ATOOLS::Cut_Data *);
    virtual void GeneratePoint(ATOOLS::Vec4D *);
    virtual void GeneratePoint(ATOOLS::Vec4D *,double *);
    virtual void GenerateWeight(ATOOLS::Vec4D *);
    virtual void GeneratePoint(double &,double &,int,double *);
    virtual void GeneratePoint(double &,double &,int);
    virtual void GenerateWeight(double,double,int);
    /*!
      Some more service methods.
    */
    virtual void SetRange(double * _sprimerange,double * _yrange);

    virtual void GetRange();
    virtual void ISRInfo(int &,double &,double &);
    virtual int CountResonances(ATOOLS::Flavour*&);
    virtual int ChNumber();
    virtual void SetChNumber(int);
    virtual std::string ChID();

    virtual int                  N()         { return n_points; }
    virtual int                  ValidN()    { return n_contrib; }
    virtual int                  Nin()       { return nin; }
    virtual int                  Nout()      { return nout; }
    virtual ATOOLS::Flavour    * Flavs()     { return fl; }
    virtual double               Res1()      { return res1; }
    virtual double               Res2()      { return res2; }
    virtual double               Res3()      { return res3; }
    virtual double               Weight()    { return weight; }
    virtual double               Alpha()     { return alpha; }
    virtual double               AlphaSave() { return alpha_save; }
    virtual double               Result()    { return result; }
    virtual std::string          Name()      { return name; }

    virtual void IncrementN()                { n_points++; }
    virtual void SetRes1(double _r)          { res1       = _r; }
    virtual void SetRes2(double _r)          { res2       = _r; }
    virtual void SetRes3(double _r)          { res3       = _r; }
    virtual void SetN(long int _n)           { n_points   = _n; }
    virtual void SetValidN(long int _n)      { n_contrib  = _n; }
    virtual void SetName(std::string _name)  { name       = _name; }
    virtual void SetWeight(double _weight)   { weight     = _weight; }
    virtual void SetAlpha(double _alpha)     { alpha      = _alpha; }
    virtual void SetAlphaSave(double _alpha) { alpha_save = _alpha; }
  };
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




