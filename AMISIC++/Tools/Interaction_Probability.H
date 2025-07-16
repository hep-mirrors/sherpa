#ifndef AMISIC_Tools_Interaction_Probability_H
#define AMISIC_Tools_Interaction_Probability_H

#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/Lookup_Tables.H"

namespace AMISIC {
  class MI_Processes;
  
  class Interaction_Probability {
    Matter_Overlap   m_mo;
    OneDim_Table   * p_k, * p_integral, * p_expO, * p_fc;
    bool   m_test;
    
    void   FixK(MI_Processes * processes);
    double Integral(const double & k,const int & diff=0);
    double NewtonRaphson(const double & xsratio);
    void   FixOExp();
    void   OutputTables(MI_Processes * processes);
  public:
    Interaction_Probability();
    ~Interaction_Probability();
    
    void   Initialize(REMNANTS::Remnant_Handler * remnant_handler,
		      MI_Processes * processes);
    
    inline Matter_Overlap * GetOverlap()                        { return &m_mo; }
    /////////////////////////////////////////////////////////////////////////////////
    // Probability for at least one interaction, Eq. (24)
    /////////////////////////////////////////////////////////////////////////////////
    inline double operator()(const double & s,const double & b) { return 1.-exp(-(*p_k)(s) * m_mo(b)); }
    /////////////////////////////////////////////////////////////////////////////////
    // b-dependent modification factor, Eq. (28)
    /////////////////////////////////////////////////////////////////////////////////
    inline double fb(const double & s,const double & b)         { return m_mo(b)/(*p_expO)(s); }
    /////////////////////////////////////////////////////////////////////////////////
    // b-dependent modification factor, Eq. (28)
    /////////////////////////////////////////////////////////////////////////////////
    inline double expO(const double & s)                        { return (*p_expO)(s); }
    /////////////////////////////////////////////////////////////////////////////////
    // Correction factor for demaning at least one interaction, Eq. (31)
    /////////////////////////////////////////////////////////////////////////////////
    inline double fc(const double & s)                          { return (*p_fc)(s); }
  };


  /////////////////////////////////////////////////////////////////////////////////
  // Interaction probability, Eq. (24) as operator, needs to be integrated for
  // Eqs.(26, 32)
  /////////////////////////////////////////////////////////////////////////////////
  class P_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    P_Integrand(Matter_Overlap * mo, const double & k) : p_mo(mo), m_k(k) {}
    ~P_Integrand() {};
    double operator()(double b);
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Integrand O(b) exp[-kO(b)], used by the Newton-Raphson method to fix k
  /////////////////////////////////////////////////////////////////////////////////
  class OtimesExp_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    OtimesExp_Integrand(Matter_Overlap * mo, const double & k) : p_mo(mo), m_k(k) {}
    ~OtimesExp_Integrand() {};
    double operator()(double b);
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Integrand of numerator in Eqs. (29, 31) as operator, given by O(b) Pint(b)
  /////////////////////////////////////////////////////////////////////////////////
  class OtimesP_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    OtimesP_Integrand(Matter_Overlap * mo, const double & k) : p_mo(mo), m_k(k) {}
    ~OtimesP_Integrand() {};
    double operator()(double b);
  };
}

#endif
