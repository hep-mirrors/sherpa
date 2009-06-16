#ifndef HADRONS_Current_Library_VA_0_PP_H
#define HADRONS_Current_Library_VA_0_PP_H

#include "HADRONS++/Current_Library/Current_Base.H"
#include "HADRONS++/PS_Library/ResonanceFlavour.H"
#include "HADRONS++/Main/Tools.H"

namespace HADRONS {
  class VA_0_PP : public Current_Base {
    //     int    m_nutau, m_pion_ch, m_pion0;    because m_pion0 = 0, m_pion_ch = 1
    bool   m_pionmode;              // either 0->2pi or 0->2kaon
    int    m_ff;                    // ff model
    double m_global;                // GF * V_CKM * Clebsch Gordon
    double m_fpi;                   // pion decay constant
    double m_frho, m_grpp;          // rho-pipi coupling
    ResonanceFlavour m_R, m_RR, m_RRR;  // resonances
    double m_beta, m_gamma;         // rel. strength
    double m_gammaR, m_gammaRR, m_gammaRRR;
    double m_m, m_m2;               // mass of outgoing particle
    double m_m2_pi, m_m2_K;         // mass of pion, kaon
    Complex A( double x, double y );
    Complex FormFactor( double s );
  public:
    VA_0_PP(ATOOLS::Flavour* flavs, int n, int* indices, std::string name) :
      Current_Base(flavs, n, indices, name) {};
    void SetModelParameters( struct GeneralModel _md );
    void Calc(const ATOOLS::Vec4D * moms);
  };
}
#endif
