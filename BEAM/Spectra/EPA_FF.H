#ifndef BEAM_Spectra_EPA_FF_H
#define BEAM_Spectra_EPA_FF_H

#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Lookup_Tables.H"

namespace ATOOLS { class Settings; }

namespace BEAM {
  enum class EPA_ff_type {
    point      = 0,
    Gauss      = 1,
    dipole     = 2,
    hcs        = 11,
    smooth_hcs = 12,
    WoodSaxon  = 13
  };
  
  class EPA_FF_Base {
  protected:
    ATOOLS::Flavour m_beam;
    double          m_mass, m_mass2, m_R;
    double          m_q2min, m_q2max, m_pt2min, m_pt2max;
    double          m_xmin, m_xmax, m_bmin, m_bmax;
    ATOOLS::OneDim_Table * p_Nred_x;
    ATOOLS::TwoDim_Table * p_N_xb, * p_Inv_xb;

    size_t m_approx, m_analytic;
    virtual void Fill_Nxb_Table(ATOOLS::axis & xaxis,ATOOLS::axis & baxis);
    virtual void Fill_Invxb_Table();
    virtual void Fill_Nredx_Table();
    virtual void TestInvxb();
  public:
    EPA_FF_Base(const ATOOLS::Flavour & beam);
    ~EPA_FF_Base() {
      if (p_Nred_x) { delete p_Nred_x; p_Nred_x = nullptr; }
      if (p_N_xb)   { delete p_N_xb;   p_N_xb   = nullptr; }
      if (p_Inv_xb) { delete p_Inv_xb; p_Inv_xb = nullptr; }
    }
    virtual void  FillTables(const size_t & nxbins=0,const size_t & nbbins=0) = 0;
    virtual const double operator()(const double & x,const double & Q2=-1.)   = 0;
    virtual const double N(const double & x)                                  = 0;
    virtual const double SelectB(const double & x=0.);
    virtual const double N(const double & x,const double & b) { return -1.; }
    virtual const double ReducedN(const double & x)           { return -1.; }
    virtual const double FF(const double & Q2)                { return 1.; }

    virtual void SetRadius(const double & R)      { m_R      = R; }
    virtual void SetQ2Min(const double & q2min)   { m_q2min  = q2min; } 
    virtual void SetQ2Max(const double & q2max)   { m_q2max  = q2max; }
    virtual void SetQ2Range(const double & q2min,const double & q2max) {
      m_q2min = q2min; m_q2max = q2max;
    }
    virtual void SetPT2Max(const double & pt2max) { m_pt2max = pt2max; }
    virtual void SetPT2Range(const double & pt2min,const double & pt2max) {
      m_pt2min = pt2min; m_pt2max = pt2max;
    }
    virtual void SetSwitch(const std::string & tag,const size_t & sw) {
      if (tag=="approximation") m_approx   = sw;
      if (tag=="analytic")      m_analytic = sw;
    }
    virtual       void   SetParam(const std::string & tag,const double & val) {}

    
    virtual const double            Radius()    const { return m_R; }
    virtual const ATOOLS::Flavour & Flav()      const { return m_beam; } 
    virtual const double            Mass()      const { return m_beam.Mass(); } 

    virtual ATOOLS::TwoDim_Table  * GetN_xb()   const { return p_N_xb; }
    virtual ATOOLS::OneDim_Table  * GetNred_x() const { return p_Nred_x; }
    
    inline virtual const double Q2min(const double & x) const {
      return ATOOLS::Max(m_q2min,ATOOLS::sqr(m_mass*x)/(1-x));
    }
    // Use maximal angle for the scattered electron in Q^2_max,
    // compare with hep-ph/9610406 and hep-ph/9310350 for the 
    // inclusion of term ~(x m)^2
    inline virtual const double Q2max(const double & x) const {
      return ( (m_pt2max>0.) ?
	       ATOOLS::Min(m_q2max,ATOOLS::sqr(m_mass*x)/(1-x)+m_pt2max*(1.-x)) :
	       m_q2max );
    }
  };

  class N_xb_int : public ATOOLS::Function_Base {
  private:
    EPA_FF_Base * p_ff;
    double m_b, m_x;
  public:
    N_xb_int(EPA_FF_Base * ff) : p_ff(ff) {}
    inline void SetXB(const double & x,const double & b) { m_x = x; m_b = b; }
    double operator()(double y);
  };

  class Nred_x_int : public ATOOLS::Function_Base {
  private:
    EPA_FF_Base * p_ff;    
    double m_x;
  public:
    Nred_x_int(EPA_FF_Base * ff) : p_ff(ff) {}
    inline void SetX(const double & x) { m_x = x; }
    double operator()(double b);
  };

  
  class EPA_Point : public EPA_FF_Base {
  private:
    void Fill_Nxb_Table(ATOOLS::axis & xaxis,ATOOLS::axis & baxis);
  public:
    EPA_Point(const ATOOLS::Flavour & beam);
    void   FillTables(const size_t & nxbins=0,const size_t & nbbins=0);
    const  double operator()(const double & x,const double & Q2=-1.);
    const  double N(const double & x);
    const  double SelectB(const double & x);
    const  double N(const double & x,const double & b);
    const  double ReducedN(const double & x);
  };

  class EPA_Dipole : public EPA_FF_Base {
  private:
    double m_Lambda2, m_mu2, m_aDip, m_bDip, m_cDip;
    void SetABC();
    const double phi(const double & y,const double & arg);
  public:
    EPA_Dipole(const ATOOLS::Flavour & beam);
    void FillTables(const size_t & m_nxbins=0,const size_t & m_nbbins=0);
    const  double operator()(const double & x,const double & Q2=-1.);
    const  double N(const double & x);
    const  double N(const double & x,const double & b);
    const  double ReducedN(const double & x);
    inline void   SetParam(const std::string & tag,const double & value) {
      if (tag=="Lambda2") m_Lambda2 = value;
      if (tag=="Mu2")     m_mu2     = value;
      SetABC();
    }
  };

  class EPA_Gauss : public EPA_FF_Base {
  private:
    size_t m_approx;
    double m_Q02, m_mu2;
  public:
    EPA_Gauss(const ATOOLS::Flavour & beam);
    void FillTables(const size_t & nxbins=0,const size_t & nbbins=0);
    const double operator()(const double & x,const double & Q2);
    const double N(const double & x);
    const double N(const double & x,const double & b);
    const double ReducedN(const double & x);
    const double FF(const double & Q2) { return exp(-Q2/m_Q02); }
    inline void   SetParam(const std::string & tag,const double & value) {
      if (tag=="Q02") m_Q02 = value;
      if (tag=="Mu2") m_mu2 = value;
    }
    inline void  SetSwitch(const std::string & tag,size_t sw) {
      if (tag=="approximation") m_approx = sw;
    }
  };

  class EPA_WoodSaxon : public EPA_FF_Base {
  private:
    class Rho_argument : public ATOOLS::Function_Base {
    private:
      double m_R, m_d;
    public:
      Rho_argument(const double & R,const double & d) : m_R(R), m_d(d) {}
      double operator()(double r) { return 4.*M_PI*r*r/(1.+exp((r-m_R)/m_d)); }
    };
    class WS_potential : public ATOOLS::Function_Base {
    private:
      double m_R, m_d, m_q;
    public:
      WS_potential(const double & R,const double & d) : m_R(R), m_d(d), m_q(1.) {}
      inline void SetQ(const double & q) { m_q = q; }
      double operator()(double r) {
	if (m_q*r<1.e-6) return 4.*M_PI*r*r/(1.+exp((r-m_R)/m_d));
	return 4.*M_PI*sin(m_q*r)/m_q*r/(1.+exp((r-m_R)/m_d));
      }
    };
    class N_argument : public ATOOLS::Function_Base {
    private:
      EPA_WoodSaxon * p_ws;
      double m_x;
    public:
      N_argument(EPA_WoodSaxon * ws) : p_ws(ws), m_x(0.) {}
      inline void SetX(const double & x) { m_x = x; }
      double operator()(double q2)       { return (*p_ws)(m_x,q2)/q2; } 
    };

    double m_d, m_rho0;
    ATOOLS::OneDim_Table * p_FF_Q2, * p_N;
    
    double CalculateDensity();
    void   InitFFTable(const double & q2min,const double & q2max);
    void   InitNTable(const double & xmin,const double & xmax);
  public:
    EPA_WoodSaxon(const ATOOLS::Flavour & beam);
    ~EPA_WoodSaxon() {
      if (p_N)     { delete p_N;      p_N     = nullptr; }
      if (p_FF_Q2) { delete p_FF_Q2;  p_FF_Q2 = nullptr; }
    }
    void FillTables(const size_t & nxbins=0,const size_t & nbbins=0);
    const double operator()(const double & x,const double & Q2);
    const double N(const double & x);
    const double N(const double & x,const double & b);
    const double ReducedN(const double & x);
    const double FF(const double & Q2) { return (*p_FF_Q2)(Q2); }
    inline void  SetParam(const std::string & tag,const double & value) {
      if (tag=="WSd") m_d = value;
    }
  };
}// end of namespace BEAM

#endif
