#ifndef AMISIC_Tools_Interaction_Probability_H
#define AMISIC_Tools_Interaction_Probability_H

#include "AMISIC++/Tools/Matter_Overlap.H"
#include "AMISIC++/Tools/Lookup_Tables.H"

namespace AMISIC {
  class MI_Processes;
  
  class Interaction_Probability {
    Matter_Overlap * p_mo;
    OneDim_Table   * p_prefKs, * p_oexpvalues, * p_integrals;
    double m_prefK, m_integral, m_integralB, m_bnorm, m_oexpvalue, m_fc;
    bool   m_output;
    
    bool FixPrefactor(MI_Processes * processes);
    void CalculateIntegral();
    void CalculateBNorm();
    void CalculateOExpValue(MI_Processes * processes);
  public:
    Interaction_Probability();
    ~Interaction_Probability();
    
    void   Initialize(MI_Processes * processes);
    void   Update(double xsecratio);
    double operator()(const double & b);
    
    inline Matter_Overlap * GetOverlap()               { return p_mo; }
    inline double Prefactor(const double & s) const    { return (*p_prefKs)(s); }
    inline double Integral(const double & s) const     { return (*p_integrals)(s); }
    inline double OverlapExpectation(const double & s) { return (*p_oexpvalues)(s); }
    inline double Bnorm(const double & s) const        { return m_bnorm/(*p_integrals)(s); }
  };


  class IP_Integrand : public ATOOLS::Function_Base {
    Interaction_Probability * p_ip;
    int m_expo;
  public:
    IP_Integrand(Interaction_Probability * ip,const int & expo=1) :
      p_ip(ip), m_expo(expo) {}
    ~IP_Integrand() {};
    double operator()(double b);
  };

  class O_ExpV_Integrand : public ATOOLS::Function_Base {
    Matter_Overlap          * p_mo;
    Interaction_Probability * p_ip;
  public:
    O_ExpV_Integrand(Interaction_Probability * ip) :
      p_ip(ip), p_mo(ip->GetOverlap()) {}
    ~O_ExpV_Integrand() {};
    double operator()(double b);
  };
}

#endif
