#ifndef Single_Channel_h
#define Single_Channel_h

#include "Message.H"
#include "Cut_Data.H"
#include "Flavour.H"
#include "Vector.H"

namespace PHASIC {

  /*!
    This class is the mother class for all channels that
    are constructed depending on the process.
  */
  class Single_Channel {
  protected:
    int                  nin, nout;
    APHYTOOLS::Flavour * fl;
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
    Single_Channel(int,int,APHYTOOLS::Flavour *);
    virtual ~Single_Channel() {
      if (ms)   { delete[] ms; }
      if (rans) { delete[] rans; }
    }
    void AddPoint(double);
    void Reset(double value);
    void ResetOpt();
    virtual void GeneratePoint(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *);
    virtual void GeneratePoint(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *,double *) {
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GeneratePoint(Vec4D,Cut_Data,double)."<<std::endl; 
    }
    virtual void GenerateWeight(AMATOOLS::Vec4D *,APHYTOOLS::Cut_Data *) {
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GenerateWeight(Vec4D,Cut_Data)"<<std::endl; 
    }

    virtual void GeneratePoint(AMATOOLS::Vec4D *);
    virtual void GeneratePoint(AMATOOLS::Vec4D *,double *) { 
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GeneratePoint(Vec4D,double)."<<std::endl; 
    }
    virtual void GenerateWeight(AMATOOLS::Vec4D *) {
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GenerateWeight(Vec4D)"<<std::endl; 
    }
    virtual void GeneratePoint(double &,double &,int,double *) {
      AORGTOOLS::msg.Error()<<"Virtual Method : "
			    <<"Single_Channel::GeneratePoint(double,double,double)"<<std::endl; 
    }
    virtual void GeneratePoint(double &,double &,int) {
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GeneratePoint(double,double)"<<std::endl; 
    }
    virtual void GenerateWeight(double,double,int) {
      AORGTOOLS::msg.Error()<<"Virtual Method : ";
      AORGTOOLS::msg.Error()<<"Single_Channel::GenerateWeight(double,double)"<<std::endl; 
    }
    /*!
      Some more service methods.
    */
    virtual void SetRange(double * _sprimerange,double * _yrange) {
      for (int i=0;i<2;i++) {
	sprimerange[i] = _sprimerange[i];
	yrange[i]      = _yrange[i];
      }
      sprimerange[2] =  _sprimerange[2];
    };

    virtual void GetRange() {
      AORGTOOLS::msg.Out()<<"  sprime : "<<sprimerange[0]<<" "<<sprimerange[1]<<" / "<<sprimerange[2]<<" / "
			  <<"  y : "<<yrange[0]<<" ... "<<yrange[1]<<std::endl;
    }

    virtual void ISRInfo(int &,double &,double &) {
      AORGTOOLS::msg.Error()<<"Virtual Method : Single_Channel::ISRInfo()"<<std::endl;
    }
    virtual int CountResonances(APHYTOOLS::Flavour*&) { 
      AORGTOOLS::msg.Error()<<"Virtual Method : Single_Channel::CountResonances()"<<std::endl; 
      return 0;
    }
    virtual int ChNumber() {
      AORGTOOLS::msg.Error()<<"Virtual Method : Single_Channel::ChNumber()"<<std::endl;
      return 0;
    }
    virtual void SetChNumber(int) {
      AORGTOOLS::msg.Error()<<"Virtual Method : Single_Channel::SetChNumber()"<<std::endl;
    }

    virtual int                  N()         { return n_points; }
    virtual int                  ValidN()    { return n_contrib; }
    virtual int                  Nin()       { return nin; }
    virtual int                  Nout()      { return nout; }
    virtual APHYTOOLS::Flavour * Flavs()     { return fl; }
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




