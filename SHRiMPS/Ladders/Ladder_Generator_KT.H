#ifndef SHRIMPS_Ladders_Ladder_Generator_KT_H
#define SHRIMPS_Ladders_Ladder_Generator_KT_H

#include "SHRiMPS/Ladders/Ladder_Generator_Base.H"
#include "SHRiMPS/Ladders/Ladder.H"
#include "SHRiMPS/Cross_Sections/Sigma_Partonic.H"
#include "MODEL/Main/Strong_Coupling.H"
#include "ATOOLS/Math/Histogram.H"

namespace SHRIMPS {
  struct rap_order {
    bool operator()(const double & y1, const double & y2) const {
      return (ATOOLS::dabs(y1)>ATOOLS::dabs(y2));
    }
  };
  
  class Ladder_Generator_KT : public Ladder_Generator_Base {
  private:
    double          m_sigmahat;
    ATOOLS::Vec4D   m_qini[2], m_sumQ, m_sumK, m_lastk;
    ATOOLS::Flavour m_flavs[2];
    std::map<double, ATOOLS::Vec4D, rap_order> m_rapidities;
    bool          MakeTrialInitialPartons();
    bool          MakeTrialLadder();
    void          MakeTrialRapidities(const size_t & N);
    void          FillRapidities();
    bool          SelectKTs();
    ATOOLS::Vec4D SelectKT(const double & y,const double & seff);
    void          ConstructPropagators();
    void          SelectPropagatorColours();
    void          CalculateWeight();
  public:
    Ladder_Generator_KT();
    ~Ladder_Generator_KT();
    
    Ladder * operator()(const ATOOLS::Vec4D & pos);    
  };
}
#endif
