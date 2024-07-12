#ifndef AMISIC_Tools_Interaction_Probability_H
#define AMISIC_Tools_Interaction_Probability_H

#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Perturbative/MI_Processes.H"
#include "AMISIC++/Tools/Lookup_Tables.H"
#include "ATOOLS/Math/Histogram.H"

namespace ATOOLS {
  class Gauss_Integrator;
}

namespace AMISIC {
  class Interaction_Probability {
    Matter_Overlap * p_mo;
    MI_Processes   * p_procs;
    axis           * p_sbins, * p_bbins;
    OneDim_Table   * p_k, * p_norm;
    TwoDim_Table   * p_diffxsec;
    double           m_xs_test, m_xs_hard, m_bmax, m_smin;
    bool             m_dynamic;

    void FixKandSmin();
    void InitializeTables();
    void InitializeTable(const size_t & sbin);
    bool CheckTables();
    void OutputTables();
    
    bool m_test, m_ana;
    std::map<std::string,ATOOLS::Histogram *> m_histos;
    void InitAnalysis();
    void Analyse();
    void FinishAnalysis();
  public:
    Interaction_Probability();
    ~Interaction_Probability();
    
    void   Initialize(Matter_Overlap * mo,MI_Processes * processes,
		      axis * sbins);
    double SelectB(const double & s);
    
    inline const double & Bmax() const { return m_bmax; }
    inline const double & Smin() const { return m_smin; }
    //////////////////////////////////////////////////////////////////////////
    // Differential hard cross section per unit area in impact parameter:
    // (*p_diffxsec)(s,b) = d\sigma_{hard}(s)/d^2b =
    // \int dx_1 dx_2 f(x_1)f(x_2) O(b;x_1,x_2) \hat\sigma(\hat{s})
    // which will become the full hard cross section when integrated over b.
    //////////////////////////////////////////////////////////////////////////
    inline double DiffXSec(const double & s,const double & b) {
      return (*p_diffxsec)(s,b);
    }
    inline double fb(const double & s,const double & b) {
      if ((*this)(s,b)<=1.e-12) return -1.;
      return (*p_diffxsec)(s,b)/p_procs->GetXSecs()->XSratio(s);
    }
    //////////////////////////////////////////////////////////////////////////
    // Probability for at least one interaction, Eq. (SZ, 24)
    //////////////////////////////////////////////////////////////////////////
    inline double operator()(const double & s,const double & b) {
      return 1.-exp(-(*p_diffxsec)(s,b));
    }
    //////////////////////////////////////////////////////////////////////////
    // The s-dependent radius multiplier
    //////////////////////////////////////////////////////////////////////////
    inline double K(const double & s) const { return (*p_k)(s); }
    //////////////////////////////////////////////////////////////////////////
    // The s-dependent norm for the diffxsec
    // (*p_norm)(s) = \int d^2b \int dx_1 dx_2 f(x_1)f(x_2) {
    //                          O(b;x_1,x_2) \hat\sigma(\hat{s})
    //                          exp[1-O(b;x_1,x_2)\hat\sigma(\hat{s})] }
    //////////////////////////////////////////////////////////////////////////
    inline double Norm(const double & s) const { return (*p_norm)(s); }
  };

  //////////////////////////////////////////////////////////////////////////
  // Integrand for a product/convolution of matter overlap and partonic
  // cross section to yield the total hard cross section:
  // int dx_1 dx_2 dpt^2 2 pi b O(b; x1,x2) dsigma^hat(s^hat,pt^2)/dpt^2
  // where the product/convolution is obtained from a pre-produced
  // look-up table.
  // This is being used for testing purposes.
  //////////////////////////////////////////////////////////////////////////
  class XS_Integrand : public ATOOLS::Function_Base {
    TwoDim_Table * p_lookup;
    double m_s;
  public:
    XS_Integrand(TwoDim_Table * lookup,const double & s) :
      p_lookup(lookup), m_s(s) {}
    ~XS_Integrand() {};
    double operator()(double b) { return 2.*M_PI*b*(*p_lookup)(m_s,b); }
    inline void SetS(const double & s) { m_s = s; }
  };

  //////////////////////////////////////////////////////////////////////////
  // Integrand for the total interaction probability constructed as one
  // minus the exponential of the negative product/convolution of matter
  // overlap and partonic cross section:
  // 2 pi b {1 - exp[-int dx_1 dx_2 dpt^2 O(b;x1,x2) dsigma^hat/dpt^2] }
  // where the product/convolution is obtained from a pre-produced
  // look-up table in the case of dynamic overlaps or directly from 
  // the matter overlap and the partonic cross section for static
  // matter overlaps.
  // This is being used to fix the (potentially s-dependent) prefactors
  // k for the matter radii of the colliding hadronic systems.
  //////////////////////////////////////////////////////////////////////////
  class PInt_Dyn_Integrand : public ATOOLS::Function_Base {
    TwoDim_Table * p_lookup;
    double m_s;
  public:
    PInt_Dyn_Integrand(TwoDim_Table * lookup) : p_lookup(lookup), m_s(0.) {}
    ~PInt_Dyn_Integrand() {};
    double operator()(double b) {
      return 2.*M_PI*b*(1. - exp(-(*p_lookup)(m_s,b)));
    }
    inline void SetS(const double & s) { m_s = s; }
  };


  //////////////////////////////////////////////////////////////////////////
  // Integrand for the normalisation of the differential hard cross
  // per unit impact parameter, i.e.
  // 2 pi b \int dx_1 dx_2 f(x_1)f(x_2) {
  //                          O(b;x_1,x_2) \hat\sigma(\hat{s})
  //                          exp[1-O(b;x_1,x_2)\hat\sigma(\hat{s})] }.
  // This is being used to create an analogue to the impact parameter
  // enhancement in the MPIs.
  //////////////////////////////////////////////////////////////////////////
  class PInt_Norm_Integrand : public ATOOLS::Function_Base {
    TwoDim_Table * p_lookup;
    double m_s;
  public:
    PInt_Norm_Integrand(TwoDim_Table * lookup) : p_lookup(lookup), m_s(0.) {}
    ~PInt_Norm_Integrand() {};
    double operator()(double b) {
      return 2.*M_PI*b*(*p_lookup)(m_s,b)*(1. - exp(-(*p_lookup)(m_s,b)));
    }
    inline void SetS(const double & s) { m_s = s; }
  };
}

#endif
