#ifndef HADRON_RESCATTERING_XSecs_NNfits_H
#define HADRON_RESCATTERING_XSecs_NNfits_H

#include "HADRON_RESCATTERING/XSecs/HPR1R2.H"

namespace HADRON_RESCATTERING {
  class NNfits {
    HPR1R2  m_hpr1r2;
    double  m_mp, m_mp2, m_mn, m_mn2;
    double  m_bp, m_bn, m_alphaP, m_s0;
  public:
    NNfits();
    ~NNfits() {}
    // pp collisions
    double pptot(const double & s);
    double ppel(const double & s);
    double ppB(const double & s);
    // pn collisions
    double pntot(const double & s);
    double pnel(const double & s);
    double pnB(const double & s);
    // ppbar collisions
    double ppbartot(const double & s);
    double ppbarel(const double & s);
    double ppbarB(const double & s);
    // B-Bbar annihilations
    double ppbarAnnihil(const double & s);
    double pnbarAnnihil(const double & s);
    // B-Bbar charge exchange
    double ppbarCEX(const double & s);
    double ppbarCEXB(const double & s);
  };
}

#endif

