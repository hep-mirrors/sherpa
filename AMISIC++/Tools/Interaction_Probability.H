#ifndef AMISIC_Tools_Interaction_Probability_H
#define AMISIC_Tools_Interaction_Probability_H

#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/Lookup_Tables.H"

namespace AMISIC {
  class MI_Processes;
  
  class Interaction_Probability {
    Matter_Overlap * p_mo;
    OneDim_Table   * p_k, * p_integral, * p_expO, * p_fc;
    bool             m_test;
    
    void   InitializeDynamic(MI_Processes * processes,axis * sbins);
    void   InitializeStatic(MI_Processes * processes,axis * sbins);
    void   FixK4Static(MI_Processes * processes,axis * sbins);
    void   FixOExp(axis * sbins);
    double Integral(const double & k,const int & diff=0);
    double NewtonRaphson(const double & xsratio);
    
    void   OutputTables(MI_Processes * processes,axis * sbins);
  public:
    Interaction_Probability();
    ~Interaction_Probability();
    
    void   Initialize(Matter_Overlap * mo,MI_Processes * processes,
		      axis * sbins);
    
    //////////////////////////////////////////////////////////////////////////
    // Probability for at least one interaction, Eq. (SZ, 24)
    //////////////////////////////////////////////////////////////////////////
    inline double operator()(const double & s,const double & b) {
      return 1.-exp(-(*p_k)(s) * (*p_mo)(b));
    }
    //////////////////////////////////////////////////////////////////////////
    // b-dependent modification factor, Eq. (SZ, 28)
    //////////////////////////////////////////////////////////////////////////
    inline double fb(const double & s,const double & b) {
      return (*p_mo)(b)/(*p_expO)(s);
    }
    //////////////////////////////////////////////////////////////////////////
    // b-dependent modification factor, Eq. (SZ, 28)
    //////////////////////////////////////////////////////////////////////////
    inline double expO(const double & s) { return (*p_expO)(s); }
    //////////////////////////////////////////////////////////////////////////
    // Correction factor for demaning at least one interaction, Eq. (SZ, 31)
    //////////////////////////////////////////////////////////////////////////
    inline double fc(const double & s) { return (*p_fc)(s); }
  };


  //////////////////////////////////////////////////////////////////////////
  // Interaction probability, Eq. (SZ, 24) as operator, needs to be 
  // integrated for Eqs.(SZ, 26 & SZ, 32)
  //////////////////////////////////////////////////////////////////////////
  class P_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    P_Integrand(Matter_Overlap * mo, const double & k) : p_mo(mo), m_k(k) {}
    ~P_Integrand() {};
    double operator()(double b);
  };

  //////////////////////////////////////////////////////////////////////////
  // Integrand for Interaction_Probablity from Eq. (CS, 13,15), needs to be 
  // integrated and used to fix k in Eqs. (CS, 14-15).
  //////////////////////////////////////////////////////////////////////////
  class nbar_DynIntegrand : public ATOOLS::Function_Base {
    MI_Processes * p_procs;
    OneDim_Table   m_table;  
  public:
    nbar_DynIntegrand(MI_Processes * procs,axis * bbins);
    ~nbar_DynIntegrand() {};
    double Initialize();
    double operator()(double b);
  };

  //////////////////////////////////////////////////////////////////////////
  // Integrand O(b) exp[-kO(b)], used by the Newton-Raphson method to fix k
  //////////////////////////////////////////////////////////////////////////
  class OtimesExp_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    OtimesExp_Integrand(Matter_Overlap * mo, const double & k) :
      p_mo(mo), m_k(k) {}
    ~OtimesExp_Integrand() {};
    double operator()(double b);
  };

  //////////////////////////////////////////////////////////////////////////
  // Integrand of numerator in Eqs. (SZ, 29 & SZ,  31) as operator,
  // given by O(b) Pint(b)
  //////////////////////////////////////////////////////////////////////////
  class OtimesP_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap * p_mo;
    double           m_k;
  public:
    OtimesP_Integrand(Matter_Overlap * mo, const double & k) :
      p_mo(mo), m_k(k) {}
    ~OtimesP_Integrand() {};
    double operator()(double b);
  };
}

#endif
