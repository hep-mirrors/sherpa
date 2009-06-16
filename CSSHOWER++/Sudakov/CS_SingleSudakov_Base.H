#ifndef CSSHOWER_Sudakov_CS_SingleSudakov_Base_H
#define CSSHOWER_Sudakov_CS_SingleSudakov_Base_H

#include "CSSHOWER++/Showers/Splitting_Function_Base.H"

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/MathTools.H"

namespace CSSHOWER {

  class CS_SingleSudakov_Base : public ATOOLS::Function_Base {
  protected:
    
    const static double s_Nc;
    const static double s_CF;
    const static double s_CA;
    const static double s_TR;
    const static double s_Tf;
    
    cstp::code m_type;
    ATOOLS::Flavour    m_fl;
    int        m_beam;
    double     m_mk2,m_Q2;
    double     m_kt2min, m_kt2max;
        
  public:
    CS_SingleSudakov_Base() : m_type(cstp::none), m_fl(ATOOLS::Flavour(kf_none)) {}
    double  Lambda(const double,const double,const double);

    virtual ~CS_SingleSudakov_Base() {}

    virtual double operator()(double kt2prod,double kt2dec) = 0;
    virtual double operator()(double kt2dec) = 0;
    virtual double Delta(double kt2prod,double kt2dec) = 0;
    
    virtual double KT2min()          const { return m_kt2min; }
    virtual double KT2max()          const { return m_kt2max; }
    virtual void   SetBeam(const int    beam)    { m_beam   = beam; }
    virtual void   SetSpecMass(const double mk2) { m_mk2    = mk2; }
    virtual void   SetQ2(const double Q2)        { m_Q2     = Q2; }
    virtual void SetKT2min(double kt2min)        { m_kt2min = kt2min; }
    virtual void SetKT2max(double kt2max)        { m_kt2max = kt2max; }

    virtual cstp::code GetType()               const { return m_type; } 
    virtual ATOOLS::Flavour GetFlavour()       const { return m_fl; } 
    virtual int        Beam()                  const { return m_beam; } 
  };
}

#endif
