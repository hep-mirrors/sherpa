#ifndef PHASIC_Channels_Single_Channel_H
#define PHASIC_Channels_Single_Channel_H

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "PHASIC++/Selectors/Cut_Data.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Info_Key.H"

namespace PHASIC {

  class Single_Channel {
  protected:
    size_t   m_nin, m_nout;
    int m_incisr, m_status;
    double * p_ms;
    size_t   m_rannum;
    double * p_rans;
    double   m_res1, m_res2, m_mres1, m_mres2;
    double   m_alpha, m_alpha_save, m_weight;
    std::string m_name;
  public:

    // constructor
    Single_Channel();
    Single_Channel(Single_Channel *);
    Single_Channel(size_t,size_t,const ATOOLS::Flavour *);
    Single_Channel(size_t,size_t,const ATOOLS::Flavour *,ATOOLS::Integration_Info * const) {};

    // destructor
    virtual ~Single_Channel();

    virtual void AddPoint(double);
    virtual void Reset(double value);
    void ResetOpt();

    virtual void Reset() {}

    virtual void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *,double *);
    virtual void GeneratePoint(ATOOLS::Vec4D *,Cut_Data *);
    virtual void GenerateWeight(ATOOLS::Vec4D *,Cut_Data *,bool=true);

    virtual void GeneratePoint(const double * rns);
    virtual void GenerateWeight(const int & mode=0);

    virtual void CalculateLimits(ATOOLS::Info_Key &spkey,ATOOLS::Info_Key &ykey);
    virtual void CalculateLimits();

    virtual inline const size_t &      Nin() const       { return m_nin; }
    virtual inline const size_t &      Nout() const      { return m_nout; }
    virtual inline const size_t &      Dimension() const { return m_rannum; }
    virtual inline const size_t        NChannels() const { return 1; }
    virtual inline const std::string & Name() const      { return m_name; }
    virtual inline const double &      Res1() const      { return m_res1; }
    virtual inline const double &      Res2() const      { return m_res2; }
    virtual inline const double &      MRes1() const     { return m_mres1; }
    virtual inline const double &      MRes2() const     { return m_mres2; }
    virtual inline const double &      Weight() const    { return m_weight; }
    virtual inline const double &      Alpha() const     { return m_alpha; }
    virtual inline const double &      AlphaSave() const { return m_alpha_save; }
    virtual inline const int           OType() const     { return 0; }
    const std::string ChID() const;


    virtual inline void SetRes1(double _r)          { m_res1       = _r; }
    virtual inline void SetRes2(double _r)          { m_res2       = _r; }
    inline void AddMPIVars(const double &r1,const double &r2) { m_mres1+=r1; m_mres2+=r2; }
    inline void SetMPIVars(const double &r1,const double &r2) { m_mres1=r1;  m_mres2=r2;  }
    virtual inline void SetName(std::string _name)  { m_name       = _name; }
    virtual inline void SetWeight(double _weight)   { m_weight     = _weight; }
    virtual inline void SetAlpha(double _alpha)     { m_alpha      = _alpha; }
    virtual inline void SetAlphaSave(double _alpha) { m_alpha_save = _alpha; }

    virtual inline void Optimize()              {}
    virtual inline void EndOptimize()           {}
    virtual inline void WriteOut(std::string)   {}
    virtual inline void ReadIn(std::string)     {}
    virtual inline bool OptimizationFinished()  { return false; }

    virtual void ISRInfo(int & type,double & mass,double & width) {
      type = 0; mass = width = 0.0;
    }
    virtual void ISRInfo(std::vector<int> &ts,
			 std::vector<double> &ms,std::vector<double> &ws) const {}
    virtual int  ChNumber(); 
    virtual void SetChNumber(int);



    virtual void MPISync();

    inline virtual void CopyMPIValues() {
      m_res1 += m_mres1;
      m_res2 += m_mres2;
      m_mres1 = m_mres2 = 0.0;
      m_status = 1;
    }
  };

    inline double MRes1() { return mres1; }
    inline double MRes2() { return mres2; }

    inline double SPrime() { return sprime; }
    inline double Y()      { return y;      }

    inline int IncludesISR() { return m_incisr; }
    inline int Status() { return m_status; }

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




