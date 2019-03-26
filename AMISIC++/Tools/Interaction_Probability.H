#ifndef AMISIC_Tools_Interaction_Probability_H
#define AMISIC_Tools_Interaction_Probability_H

#include "AMISIC++/Tools/Matter_Overlap.H"

namespace AMISIC {
  class Interaction_Probability {
    Matter_Overlap * p_mo;
    double m_prefK, m_bmax, m_integral, m_oexpvalue, m_fc;
    
    bool FixPrefactor(const double & xsecratio);
    void CalculateIntegral();
    void CalculateOExpValue();
  public:
    Interaction_Probability();
    ~Interaction_Probability() { delete p_mo; }

    void   Initialize(const double & xsecratio);
    double operator()(const double & b);

    inline void SetPrefactor(const double & prefK)   { m_prefK = prefK; }
    inline const double & Prefactor() const          { return m_prefK; }
    inline const double & Integral() const           { return m_integral; }
    inline const double & OverlapExpectation() const { return m_oexpvalue; }
    inline Matter_Overlap * GetOverlap()             { return p_mo; }
  };


  class IP_Integrand : public ATOOLS::Function_Base {
    Interaction_Probability * p_ip;
  public:
    IP_Integrand(Interaction_Probability * ip) : p_ip(ip) {}
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
