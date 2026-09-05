// The EEX (exclusive exponentiation) beta expansion for a single dipole.
//
// Split out of Dipole.C, which was 968 lines. This block - EEX, Beta1/2/3,
// VirtualEEX, the three Hard overloads and the two xi helpers - is a
// self-contained unit: it reads the dipole's eikonal momenta, gamma and photon
// list and returns a number. It touches no frame/boost state and no photon
// bookkeeping, which is what the rest of Dipole.C is. Same class, second
// translation unit, so the physics path is unchanged - only where it is
// written down. Eikonal() and CalculateGamma() stay in Dipole.C.
//
// The expansion. beta_n is the n-photon term with all its lower-order
// subtractions removed, so that summing beta_1..beta_n over the generated
// photons and dividing by the crude eikonal weight gives the correction to the
// Born. Each beta_n therefore evaluates beta_(n-1), beta_(n-2), ... as part of
// building itself, and it evaluates them at LOWER orders - a subtraction term
// must not carry the virtual dressing of the term it is subtracted from.
//
// That "order to evaluate at" is a parameter here. It used to be the member
// m_betaorder, which the callees wrote: Beta2 decremented it and restored it,
// Beta3 set it to -1 and then -2 and never restored it, and EEX had to
// reassign it after every pair and every triple to repair the damage before
// the next iteration. Two consequences, both now gone: the reassignments in
// EEX's loops were load-bearing rather than defensive, and any caller reaching
// Beta1() without going through EEX() first - Define_Dipoles::EEXRealVirtual
// and NLO_Base's collinear-ratio histogram both do - read whatever order the
// last unrelated call had left behind.

#include "YFS/Main/Dipole.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace ATOOLS;
using namespace YFS;

// Defined in Dipole.C, written there by CalculateGamma(). NOT per-dipole: deli
// is the initial-state 0.5*gamma and delf the final-state one, and Beta1/Beta2
// below read BOTH whatever this dipole's own type is - an ISR dipole's beta_1
// carries the FSR dressing and vice versa. See the comment at their
// definition before considering making them members; doing so zeroes the cross
// term.
extern double delf;
extern double deli;

namespace {
  // Orders used for the subtraction terms inside Beta3. Any value below 2
  // disables the virtual dressing in Hard() and the delta factors in Beta1/
  // Beta2, which is the point; the specific negatives are the ones the
  // member-mutating version left behind, kept so the arithmetic is unchanged.
  constexpr int SUBTRACTION_ORDER_BETA2 = -1;
  constexpr int SUBTRACTION_ORDER_BETA1 = -2;
}

double Dipole::EEX(const int betaorder){
  double real=0;
  if(m_dipolePhotonsEEX.size()==0) return real;
  CalculateGamma();
  // NB the single-photon term runs over m_dipolePhotons while the pair and
  // triple terms run over m_dipolePhotonsEEX. For an ISR dipole the two hold
  // the same vectors (Boost() assigns one to the other). For an FSR dipole
  // they do NOT: Boost() pushes each photon into m_dipolePhotonsEEX and then
  // boosts it in m_dipolePhotons, so beta_1 sees boosted photons and
  // beta_2/beta_3 see the pre-boost ones. Preserved as found - correcting it
  // would move numbers - but it is not obviously intended.
  if(betaorder >= 1) {
    for(auto k: m_dipolePhotons){
      real += Beta1(k, betaorder)/Eikonal(k);
    }
  }
  if(betaorder >= 2 ) {
    for (size_t j = 1; j < m_dipolePhotonsEEX.size(); j++) {
      for (size_t i = 0; i < j; i++) {
        const Vec4D &k1 = m_dipolePhotonsEEX[j];
        const Vec4D &k2 = m_dipolePhotonsEEX[i];
        real += Beta2(k1,k2,betaorder)/Eikonal(k1)/Eikonal(k2);
      }
    }
  }
  if(betaorder >= 3){
    for (size_t j = 1; j < m_dipolePhotonsEEX.size(); j++) {
      for (size_t i = 0; i < j; i++) {
        for (size_t k = 0; k < i; k++) {
          const Vec4D &k1 = m_dipolePhotonsEEX[j];
          const Vec4D &k2 = m_dipolePhotonsEEX[i];
          const Vec4D &k3 = m_dipolePhotonsEEX[k];
          real += Beta3(k1,k2,k3,betaorder)
                  /Eikonal(k1)/Eikonal(k2)/Eikonal(k3);
        }
      }
    }
  }
  if(IsNan(real)){
    msg_Error()<<"YFS EEX is NaN at order "<<betaorder<<std::endl;
  }
  return real;
}

double Dipole::Beta1(const Vec4D &k, int order){
  const double hard = Hard(k, order);
  const double eik  = Eikonal(k);
  if(Type()==dipoletype::initial) {
    // beta11
    if(order==2) return hard*(1+delf)-eik*(1+deli)*(1+delf);
    // beta12
    if(order==3) return hard-eik*(1+delf+0.5*delf*delf)
                             *(1+deli+0.5*deli*deli);
    return (hard-eik)*(1+delf);
  }
  if(Type()==dipoletype::final) {
    if(order==2) return hard*(1+deli)-eik*(1.+deli)*(1+delf);
    if(order==3) return hard*(1+deli+0.5*deli*deli)
                        -eik*(1+delf+0.5*delf*delf)
                            *(1+deli+0.5*deli*deli);
    return hard-eik;
  }
  return hard-eik;
}

double Dipole::Beta2(const Vec4D &k1, const Vec4D &k2, int order){
  const double eik1 = Eikonal(k1);
  const double eik2 = Eikonal(k2);
  const double delta = (order==3) ? (1+delf)*(1+deli) : 1.;
  // The two single-photon subtractions go in one order lower: they remove
  // what beta_1 already accounted for, without its dressing.
  //
  // Summed as `hard += (-a-b-c)` rather than `Hard()-a-b-c` to keep the
  // floating-point association the member-mutating version had: the two
  // group the additions differently and would disagree in the last bits.
  double hard = Hard(k1,k2);
  hard += -eik1*Beta1(k2, order-1)
          -eik2*Beta1(k1, order-1)
          -eik1*eik2;
  return hard*delta;
}

double Dipole::Beta3(const Vec4D &k1, const Vec4D &k2, const Vec4D &k3, int order){
  if(Type()!=dipoletype::initial) return 0;
  const double eik1 = Eikonal(k1);
  const double eik2 = Eikonal(k2);
  const double eik3 = Eikonal(k3);
  // Two separate `+=` groups, as before: see the note in Beta2 on keeping the
  // floating-point association.
  double hard = Hard(k1,k2,k3);
  hard += -eik1*Beta2(k3,k2,SUBTRACTION_ORDER_BETA2)
          -eik2*Beta2(k3,k1,SUBTRACTION_ORDER_BETA2)
          -eik3*Beta2(k1,k2,SUBTRACTION_ORDER_BETA2)
          -eik1*eik2*eik3;
  hard += -eik2*eik3*Beta1(k1,SUBTRACTION_ORDER_BETA1)
          -eik1*eik3*Beta1(k2,SUBTRACTION_ORDER_BETA1)
          -eik1*eik2*Beta1(k3,SUBTRACTION_ORDER_BETA1);
  return hard;
}

double Dipole::VirtualEEX(const int betaorder){
  // For ISR+FSR virtuals are taken in for ISRxFSR not ISR+FSR
  double virt{0};
  if(betaorder==1)      virt = 0.5*m_gamma;
  else if(betaorder==2) virt = 0.5*m_gamma + 0.125*m_gamma*m_gamma;
  else if(betaorder==3) virt = 0.5*m_gamma + 0.125*m_gamma*m_gamma
                               + pow(m_gamma,3)/48;
  return virt;
}

double Dipole::Hard(const Vec4D &k, int order){
  const double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];
  const double a = k*m_eikmomentum[0]/p1p2;
  const double b = k*m_eikmomentum[1]/p1p2;
  const double ap = a/(1.+a+b);
  const double bp = b/(1.+a+b);
  double delta = 0;
  if (Type() == dipoletype::initial) {
    const double z = (1-a)*(1-b);
    if(order>=2 && !RealOnly()){
      delta += 0.5*m_gamma
              +m_alpi*(log(a)*log(1-b)+log(b)*log(1-a)
                      +DiLog(a) + DiLog(b)
                      -0.5*sqr(log(1-a))-0.5*sqr(log(1-b))
                      +1.5*log(1-a)+1.5*log(1-b)
                      +0.5*a*(1-a)/(1+sqr(1-a))
                      +0.5*b*(1-b)/(1+sqr(1-b)));
    }
    if(order>=3){
      delta += 0.125*sqr(m_gamma)*(1-log(z))
             +sqr(m_gamma)/24 *sqr(log(z));
    }
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-b))*(1+delta);
  }
  if (Type() == dipoletype::final) {
    const double z = (1-ap)*(1-bp);
    if(order>=2){
      delta += 0.5*m_gamma+0.25*m_gamma*log(z);
    }
    return 0.5*Eikonal(k)*(sqr(1-ap)+sqr(1-bp))*(1+delta);
  }
  if (Type() == dipoletype::ifi) {
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-bp));
  }
  return 0;
}

double Dipole::Hard(const Vec4D &k1, const Vec4D &k2){
  const double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];

  const double a1 = k1*m_eikmomentum[0]/p1p2;
  const double a2 = k2*m_eikmomentum[0]/p1p2;

  const double b1 = k1*m_eikmomentum[1]/p1p2;
  const double b2 = k2*m_eikmomentum[1]/p1p2;

  const double eta1 = a1/(1+a1+b1);
  const double eta2 = a2/(1+a2+b2);

  const double zeta1 = b1/(1+a1+b1);
  const double zeta2 = b2/(1+a2+b2);

  const double etap1 = eta1/(1+eta2);
  const double etap2 = eta2/(1+eta1);

  const double zetap1 = zeta1/(1+zeta2);
  const double zetap2 = zeta2/(1+zeta1);

  const double ap1 = a1/(1.-a2);
  const double bp1 = b1/(1.-b2);

  const double ap2 = a2/(1.-a1);
  const double bp2 = b2/(1.-b1);

  const double v1 = a1+b1;
  const double v2 = a2+b2;
  double hard;
  if (Type() == dipoletype::initial) {
    if(v1 > v2) hard = xi(a1,ap2,bp2) + xi(b1,ap2,bp2);
    else        hard = xi(a2,ap1,bp1) + xi(b2,ap1,bp1);
    return Eikonal(k1)*Eikonal(k2)*hard;
  }
  if (Type() == dipoletype::final) {
    if(v1 > v2) hard = xi(eta1,etap2,zetap2) + xi(zeta1,etap2,zetap2);
    else        hard = xi(eta2,etap1,zetap1) + xi(zeta2,etap1,zetap1);
    return Eikonal(k1)*Eikonal(k2)*hard;
  }
  // No two-photon ifi term. The third branch here used to test
  // dipoletype::initial a second time, so it was unreachable, and it fell
  // through to the return below without using what it computed.
  return 0;
}

double Dipole::Hard(const Vec4D &k1, const Vec4D &k2, const Vec4D &k3){
  const double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];

  const double a1 = k1*m_eikmomentum[0]/p1p2;
  const double a2 = k2*m_eikmomentum[0]/p1p2;
  const double a3 = k3*m_eikmomentum[0]/p1p2;

  const double b1 = k1*m_eikmomentum[1]/p1p2;
  const double b2 = k2*m_eikmomentum[1]/p1p2;
  const double b3 = k3*m_eikmomentum[1]/p1p2;

  const double eta1 = a1/(1+a1+b1);
  const double eta2 = a2/(1+a2+b2);
  const double eta3 = a3/(1+a3+b3);

  const double zeta1 = b1/(1+a1+b1);
  const double zeta2 = b2/(1+a2+b2);
  const double zeta3 = b3/(1+a3+b3);

  const double etap1 = eta1/(1+eta2);
  const double zetap1 = zeta1/(1+zeta2);

  const double etap2 = eta2/(1+eta1);
  const double zetap2 = zeta2/(1+zeta1);

  const double etap3  = eta3/(1+eta1+eta3);
  const double zetap3 = zeta3/(1+zeta1+zeta2);

  const double ap1 = a1/(1.-a2);
  const double bp1 = b1/(1.-b2);

  const double ap2 = a2/(1.-a1);
  const double bp2 = b2/(1.-b1);

  const double ap3 = a3/(1-a1-a2);
  const double bp3 = b3/(1-b1-b2);

  const double v1 = a1+b1;
  const double v2 = a2+b2;
  double hard;
  if (Type() == dipoletype::initial) {
    if(v1 > v2) hard = xi(a1,ap2,bp2,ap3,bp3) + xi(b1,ap2,bp2,ap3,bp3);
    else        hard = xi(a2,ap1,bp1,ap3,bp3) + xi(b2,ap1,bp1,ap3,bp3);
    return Eikonal(k1)*Eikonal(k2)*Eikonal(k3)*hard;
  }
  if (Type() == dipoletype::final) {
    if(v1 > v2) hard = xi(eta1,etap2,zetap2,etap3,zetap3)
                       + xi(zeta1,etap2,zetap2,etap3,zetap3);
    else        hard = xi(eta2,etap1,zetap1,etap3,zetap3)
                       + xi(zeta2,etap1,zetap1,etap3,zetap3);
    return Eikonal(k1)*Eikonal(k2)*Eikonal(k3)*hard;
  }
  return 0;
}

double Dipole::xi(const double &alp, const double &beta, const double &gamma){
  return 0.25*sqr(1.-alp)*(sqr(1.-beta)+sqr(1.-gamma));
}

double Dipole::xi(const double &alp, const double &a1, const double &b1, const double &a2, const double &b2){
  return 0.125*sqr(1.-alp)*(sqr(1.-a1)+sqr(1.-b1))*(sqr(1.-a2)+sqr(1.-b2));
}
