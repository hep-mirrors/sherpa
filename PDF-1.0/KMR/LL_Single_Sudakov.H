#ifndef LL_Single_Sudakov_H
#define LL_Single_Sudakov_H

#include "NLL_Branching_Probability_Base.H" 
#include "NLL_Sudakov_Base.H"
#include "Info_Key.H"

namespace SHERPA {

  class LL_Single_Sudakov : public NLL_Sudakov_Base  {
  private:

    double m_qmax, m_qmin;

    Sudakov::code m_calcmode;
    Sudakov::code m_cutmode;

    NLL_Branching_Probability_Base *p_bp;
    ATOOLS::Gauss_Integrator        m_gauss;

    std::string      m_outpath;
    ATOOLS::Info_Key m_lastsud;

  public:

    // constructor
    LL_Single_Sudakov(NLL_Branching_Probability_Base *,int);

    // destructor
    ~LL_Single_Sudakov();

    // member functions 
    bool Initialize(double qmin=0.,double qmax=0.);

    double operator()(double Q,double q);
    double operator()(double Q);

    double Log(double Q,double q);

    void AssignKey(ATOOLS::Integration_Info *const info);

    // inline functions 
    inline double Qmax() { return m_qmax; }
    inline double Qmin() { return m_qmin; }
    
    inline void SetQmax(double q) { m_qmax = q; p_bp->SetQmax(q); }
    inline void SetQmin(double q) { m_qmin = q; p_bp->SetQmin(q); }
    
    inline void   SetOutPath(const std::string _m_outpath) { m_outpath=_m_outpath;       }
    inline double IntGamma(double q,double Q)              { return p_bp->IntGamma(q,Q); }

  };

}// end of namespace SHERPA

#endif
