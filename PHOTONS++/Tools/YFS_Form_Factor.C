#include "PHOTONS++/Tools/YFS_Form_Factor.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "PHOTONS++/Main/Photons.H"

#define LOG_2 0.69314718055994530942

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

namespace PHOTONS {

  class IG1: public ATOOLS::Function_Base {
  private:

    YFS_Form_Factor *p_ff;

  public:

    inline IG1(YFS_Form_Factor *const ff): p_ff(ff) {}

    double operator()(double x)
    {
      Vec4D px1 = 0.5*((p_ff->P1()+p_ff->P2())-x*(p_ff->P1()-p_ff->P2()));
      Vec4D px2 = 0.5*((p_ff->P1()+p_ff->P2())+x*(p_ff->P1()-p_ff->P2()));
      return p_ff->G(-x)/(px1*px1)+p_ff->G(x)/(px2*px2);
    }
    double operator()() { return m_defval; }

  };// end of class IG1

  class IG2: public ATOOLS::Function_Base {
  private:

    YFS_Form_Factor *p_ff;

  public:

    inline IG2(YFS_Form_Factor *const ff): p_ff(ff) {}

    double operator()(double x)
    {
      Vec4D px = 0.5*((p_ff->P1()+p_ff->P2())+x*(p_ff->P1()-p_ff->P2()));
      return p_ff->G(x)/(px*px);
    }
    double operator()() { return m_defval; }

  };// end of class IG2

}// end of namespace PHOTONS

YFS_Form_Factor::YFS_Form_Factor(const Particle_Vector& part,
                                 const double& ks) {
  p_ig1 = new IG1(this);
  p_ig2 = new IG2(this);
  p_gi1 = new Gauss_Integrator(p_ig1);
  p_gi2 = new Gauss_Integrator(p_ig2);

  double YSum = 0.;
  for (unsigned int j=0; j<part.size(); j++)
    for (unsigned int i=0; i<j; i++)
      YSum+=YFS_Form_Factor(part[i],part[j],ks).Get();

  m_Y = YSum;
}

YFS_Form_Factor::YFS_Form_Factor(const Particle * part1, const Particle * part2,
                                 const double& ks) {
  p_ig1 = new IG1(this);
  p_ig2 = new IG2(this);
  p_gi1 = new Gauss_Integrator(p_ig1);
  p_gi2 = new Gauss_Integrator(p_ig2);
  m_ks = ks;

  m_p1 = part1->Momentum();
  m_p2 = part2->Momentum();

  m_m1 = part1->FinalMass();
  m_m2 = part2->FinalMass();

  m_Z1 = part1->Flav().Charge();
  m_Z2 = part2->Flav().Charge();

  // choose such that E_2 >= E_1
  if (part2->Momentum()[0] >= part1->Momentum()[0]) {
    std::swap(m_p1,m_p2);
    std::swap(m_m1,m_m2);
    std::swap(m_Z1,m_Z2);
  }

  if (part1->ProductionBlob() == part2->ProductionBlob())  m_t1t2 = +1.;
  else if (part1->ProductionBlob() == part2->DecayBlob())  m_t1t2 = -1.;
  else if (part1->DecayBlob() == part2->ProductionBlob())  m_t1t2 = -1.;
  else if (part1->DecayBlob() == part2->DecayBlob())       m_t1t2 = +1.;
  else                                                     m_t1t2 = 0.;

  // roots of p_x^2, special case for (p1-p2)^2=0 needed only for W->lnu
  if (!(m_t1t2==-1 && abs((m_p1-m_p2).Abs2()) < 1E-6)) {
    m_x1  = - (m_p1.Abs2() - m_p2.Abs2()
                + 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                  / (m_p1-m_p2).Abs2();
    m_x2  = - (m_p1.Abs2() - m_p2.Abs2()
                - 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                  / (m_p1-m_p2).Abs2();
    if (abs(1.+m_x1) < 1E-10) {
      msg_Error()<<METHOD<<"() error: case should not appear !!!"<<endl;
    }
    if (abs(1.+m_x2) < 1E-10) {
      // expansion for (p1-p2)<~0 => x2>~-1
      double p22(m_p2.Abs2());
      double x(m_p1.Abs2()*m_p2.Abs2()/sqr(m_p1*m_p2));
      double r(m_p1*m_p2);
      m_x2  = - 1. + (2.*p22-r*(-0.5*x+0.125*x*x-0.0625*x*x*x))
                      / (m_p1-m_p2).Abs2();
    }
  }
  else m_x1=m_x2=0.;

  //roots of p_x^'2
  if (m_t1t2 == +1.){
    m_xx1 = - (m_p1.Abs2() - m_p2.Abs2()
                + 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                  / (m_p1+m_p2).Abs2();
    m_xx2 = - (m_p1.Abs2() - m_p2.Abs2()
                - 2.*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                  / (m_p1+m_p2).Abs2();
  }
  else m_xx1=m_xx2=0.;

  // calculate
  m_Y = Y();
}


YFS_Form_Factor::~YFS_Form_Factor() 
{
  delete p_gi1;
  delete p_gi2;
  delete p_ig1;
  delete p_ig2;
}

// private members

// Y = 2 alpha (Re B + B~)
//   = - alpha/pi Z_iZ_j theta_itheta_j
//     * [ log(E_iE_j/omega^2)
//         + (pipj)/2 * ( theta_itheta_j \int_-1^1 dx log(px'^2/lambda^2)/px'^2
//                        + \int_-1^1 dx log(px^2/lambda^2)/px^2 )
//         - (pipj)/2 * \int_-1^1 dx log(Ex^2/omega^2)/px^2
//         + 1/4 * \int_-1^1 dx log(px'^2/m1m2)
//         + G(1) + G(-1) - (pipj) * \int_-1^1 dx G(x)/px^2 ]             (C.44)
double YFS_Form_Factor::Y() {
  return -Photons::s_alpha/M_PI*m_Z1*m_Z2*m_t1t2
           *(log((m_p1[0]*m_p2[0])/(m_ks*m_ks))
             + 0.5*(m_p1*m_p2)*IntP1() - 0.5*(m_p1*m_p2)*IntE()
             + 0.25*IntP2()
             + G(1.) + G(-1.) - (m_p1*m_p2)*IntG());
}

// t1t2 * int dx ln (px'²/lambda²)/px'² + int dx ln (px²/lambda²)/px²     (C.47)
double YFS_Form_Factor::IntP1() {
  if (m_t1t2 == -1.) {
    return 0.;
  }
  else if (m_t1t2 == +1.) {
    double A(0.);
    if (m_xx1*m_xx2 >= 0.) 
      A = 8.*(M_PI*M_PI)/((m_xx2-m_xx1)*(m_p1+m_p2).Abs2());

    double B = 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
               *(log(abs(m_x1))*(DiLog((m_x1-1.)/m_x1)-DiLog((m_x1+1.)/m_x1))
                 -log(abs(m_x2))*(DiLog((m_x2-1.)/m_x2)-DiLog((m_x2+1.)/m_x2)));
    return (A+B);
  }
  else {
    return 0.;
  }
}

// int dx ln(Ex²/omega²)/px²                                              (C.61)
double YFS_Form_Factor::IntE() {
  // E1 = E2
  if (abs(m_p1[0]-m_p2[0]) < 1E-6) {
    return 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
           *log((m_p1[0]+m_p2[0])/(2*m_ks))
           *log(abs(((1.-m_x1)*(1.+m_x2))/((1.+m_x1)*(1.-m_x2))));
  }
  else {
    // (p1-p2)^2 < 0
    if ((m_p1-m_p2).Abs2() < -1E-6) {
      double xE   = -(m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
      double zeta = -(m_p1[0]-m_p2[0])/(2.*m_p1[0]);
      double y1   = 1.+zeta*(1.-m_x1);
      return 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
              * (log(m_p1[0]/m_ks)*log(abs((1.-m_x1)/(1.+m_x1)))
                 + log(abs(y1))*log(abs((1.-y1)/(1.-y1+2.*zeta)))
                 + DiLog(1.-((1.+2.*zeta)/y1)) - DiLog(1.-(1./y1))
                 - DiLog(-(1.+m_x2)/(xE-m_x2))
                 + DiLog((1.-m_x2)/(xE-m_x2))
                 - log((1./(2.*m_ks))*((1.+m_x2)*m_p1[0]+(1.-m_x2)*m_p2[0]))
                    *log(abs((1.-m_x2)/(1.+m_x2))));
    }
    // (p1-p2)^2 > 0
    else if ((m_p1-m_p2).Abs2() > 1E-6) {
      // m2^2 > 2(p1*p2)-m1^2
      if ((m_m2*m_m2) > (2.*(m_p1*m_p2)-m_m1*m_m1)) {
        double xi   = (m_p1[0]-m_p2[0])/(2.*m_p2[0]);
        double zeta = -(m_p1[0]-m_p2[0])/(2.*m_p1[0]);
        double y1   = 1.+zeta*(1.-m_x1);
        double y2   = 1.+xi*(1.+m_x2);
        return 8./((m_p1-m_p2).Abs2()*(m_x1-m_x2))
                * (log(m_p1[0]/m_ks)*log(abs((1.-m_x1)/(1.+m_x1)))
                   + log(abs(y1))*log(abs((1.-y1)/(1.-y1+2.*zeta)))
                   + DiLog(1.-((1.+2.*zeta)/y1)) - DiLog(1.-(1./y1))
                   - DiLog((y2)/(y2-1.-2.*xi)) + DiLog((y2)/(y2-1.))
                   + 0.5*(sqr(log(abs(y2/(y2-1.))))
                           - sqr(log(abs(y2/(y2-1.-2.*xi)))))
                   - log(abs(y2))*log(abs((y2+1.+2.*xi)/(y2+1.)))
                   - log(m_p2[0]/m_ks)*log(abs((1.-m_x2)/(1.+m_x2))));
      }
      // m12 > 2(p1*p2)-m2^2
      else if ((m_m1*m_m1) > (2.*(m_p1*m_p2)-m_m2*m_m2)) {
        msg_Error()<<METHOD<<"() error: case should not appear !!!"<<endl;
        return 0.;
      }
      else {
        msg_Error()<<METHOD<<"() error: case should not appear !!!"<<endl;
        return 0.;
      }
    }
    // (p1-p2)^2 = 0
    else {
      // m1 = m2
      if (abs(m_m1-m_m2) < 1E-6) {
        // E1 != E2
        if (abs(m_p1[0]-m_p2[0]) > 1E-6) {
          double xE = - (m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
          return 8./(m_p1+m_p2).Abs2()
                  *(2.*log(1./2.*(m_p2[0]-m_p1[0])/m_ks)
                    -xE*log((xE-1.)/(xE+1.))-log(xE*xE-1.)-2.);
        }
        else {
          msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
          return 0.;
        }
      }
      else {
        // E1 != E2
        if (abs(m_p1[0]-m_p2[0]) > 1E-6) {
          double xE = - (m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
          double xp = - (m_m1*m_m1+m_m2*m_m2)/(m_m1*m_m1-m_m2*m_m2);
          // xE = xp
          if (abs(xE-xp) < 1E-6)
            return 4./(m_m2*m_m2-m_m1*m_m1)
                    *(log((m_p2[0]-m_p1[0])/(2.*m_ks))*log(abs((xp+1.)/(xp-1.)))
                      - (1./2.)*(sqr(log(xp-1.))-sqr(log(xp+1.))));
          // xE > xp
          else if (xE > xp)
            return 4./(m_m2*m_m2-m_m1*m_m1)
                    *(log((m_p2[0]-m_p1[0])/(2.*m_ks))*log(abs((xp+1.)/(xp-1.)))
                      + log(xE-xp)*log(abs((xp+1.)/(xp-1.)))
                      + DiLog((xp-1.)/(xp-xE)) - DiLog((xp+1.)/(xp-xE)));
          // xE < xp
          else if (xp > xE) {
            double xi = (m_p1[0]-m_p2[0])/(2*m_p2[0]);
            double yp = 1.+xi*(1.+xp);
            return 4./(m_m2*m_m2-m_m1*m_m1)
                    *((log(m_p2[0]/m_ks)+log(abs(yp)))*log(abs((xp+1.)/(xp-1.)))
                      - (1./2.)*sqr(log(abs(yp/(yp-1.-2.*xi))))
                      + (1./2.)*sqr(log(abs(yp/(yp-1.))))
                      + DiLog(yp/(yp-1.)) - DiLog(yp/(yp-1.-2.*xi)));
          }
          else {
            msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
            return 0.;
          }
        }
        else {
          msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
          return 0.;
        }
      }
    }
  }
}

// int dx ln (px'²/m1m2)                                                  (C.55)
double YFS_Form_Factor::IntP2() {
  if(m_t1t2 == +1.)
    return 2.*log((m_p1+m_p2).Abs2()/(4.*m_m1*m_m2))
            + log(abs((1.-m_xx1*m_xx1)*(1.-m_xx2*m_xx2)))
            - m_xx1*log(abs((1.-m_xx1)/(1.+m_xx1)))
            - m_xx2*log(abs((1.-m_xx2)/(1.+m_xx2))) - 4.;
  else if (m_t1t2 == -1.) {
    // (p1-p2)^2 != 0
    if (abs((m_p1-m_p2).Abs2()) > 1E-6)
      return 2.*log(abs((m_p1-m_p2).Abs2())/(4.*m_m1*m_m2))
              + log(abs((1.-m_x1*m_x1)*(1.-m_x2*m_x2)))
              - m_x1*log(abs((1.-m_x1)/(1.+m_x1)))
              - m_x2*log(abs((1.-m_x2)/(1.+m_x2))) - 4.;
    // (p1-p2)^2 = 0
    else {
      // m1^2 != m2^2
      if (abs(m_m1*m_m1 - m_m2*m_m2) > 1E-6) {
        double xp = -(m_m1*m_m1+m_m2*m_m2)/(m_m1*m_m1 - m_m2*m_m2);
        return 2.*log(abs(m_m1*m_m1-m_m2*m_m2)/(2.*m_m1*m_m2))
                + log(abs(1.-xp*xp)) + xp*log(abs((1.+xp)/(1.-xp))) - 2.;
      }
      // m1^2 = m2^2
      else
        return 2.*log((m_p1+m_p2).Abs2()/(4.*m_m1*m_m2));
    }
  }
  else {
    msg_Error()<<METHOD<<"(): error: case should not appear !!!"<<endl;
    return 0.;
  }
}

// G(x)
double YFS_Form_Factor::G(double x) {
  Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
  double b  = CalculateBeta(px);
  double r(0.);
  if (b == 0.)      r = 1.-LOG_2;
  else if (b == 1.) r = 0.;
  else              r = ((1.-b)/(2.*b)*log((1.+b)/(1.-b))+log((1.+b)/2.));
  return r;
}

// Function for evaluating IntG() for dipole of different masses in its
// rest frame
double YFS_Form_Factor::GFunc(double x) {
  double xE = (m_p2[0]+m_p1[0])/(m_p2[0]-m_p1[0]);
  double xD = 2.*Vec3D(m_p1).Abs()/(m_p2[0]-m_p1[0]);

  double r1 = xE/(2.*xD*m_x1*(m_x1-m_x2))
                  *(-1./2.*sqr(log(abs(x-m_x1)))
                    + log(abs(m_x1/m_x2))*log(abs(x-m_x1))
                    + log(abs(x-m_x2))*log(abs((x-m_x1)/(m_x2-m_x1)))
                    + DiLog((m_x2-x)/(m_x2-m_x1)))
            + xE/(2.*xD*m_x2*(m_x1-m_x2))
                  *(-1./2.*sqr(log(abs(x-m_x2)))
                    - log(abs(m_x1/m_x2))*log(abs(x-m_x2))
                    + log(abs(x-m_x1))*log(abs((x-m_x2)/(m_x1-m_x2)))
                    + DiLog((m_x1-x)/(m_x1-m_x2)))
            + xE/(2.*xD*m_x1*m_x2)
                  *(sqr(log(abs(m_x1)))
                    - sqr(log(abs(m_x2)))
                    + DiLog(x/m_x1)
                    - DiLog(x/m_x2));

  double r2 = -(1.+xD)/(2.*xD)
                *(-sqr(log(abs(m_x1/m_x2*(m_x2-x)/(m_x1-x))))/(2.*(m_x1-m_x2)));

  double r3 = 1./(m_x1-m_x2)*(1./2.*sqr(log(abs(x-m_x2)))
                              + log(abs((1.-xD)/2.))*log(abs((x-m_x1)/(x-m_x2)))
                              - log(abs(x-xE))*log(abs((m_x1-x)/(m_x1-xE)))
                              + log(abs(x-xE))*log(abs((m_x2-x)/(m_x2-xE)))
                              + log(abs(x-m_x2))*log(abs((m_x1-x)/(m_x1-m_x2)))
                              - DiLog((x-xE)/(m_x1-xE))
                              + DiLog((x-xE)/(m_x2-xE))
                              + DiLog((x-m_x2)/(m_x1-m_x2)));

  double r4 = xE/(2.*xD*m_x1*(m_x1-m_x2))
                  *(1./2.*sqr(log(abs(x+m_x1)))
                    + log(abs(m_x2/m_x1))*log(abs(x+m_x1))
                    - log(abs(x+m_x2))*log(abs((m_x1+x)/(m_x1-m_x2)))
                    - DiLog((m_x2+x)/(m_x2-m_x1)))
            + xE/(2.*xD*m_x2*(m_x1-m_x2))
                  *(1./2.*sqr(log(abs(x+m_x2)))
                    - log(abs(m_x2/m_x1))*log(abs(x+m_x2))
                    - log(abs(x+m_x1))*log(abs((m_x2+x)/(m_x2-m_x1)))
                    - DiLog((m_x1+x)/(m_x1-m_x2)))
            + xE/(2.*xD*m_x1*m_x2)*(DiLog(-x/m_x1)-DiLog(-x/m_x2));

  double r5 = (1.-xD)/(2.*xD)
                *(-sqr(log(abs(m_x2/m_x1*(m_x1+x)/(m_x2+x))))/(2.*(m_x1-m_x2)));
  double r6 = 1./(m_x1-m_x2)*(-1./2.*sqr(log(abs(x+m_x1)))
                              - log(abs((1.-xD)/2.))*log(abs((x+m_x1)/(x+m_x2)))
                              + log(abs(x+xE))*log(abs((m_x1+x)/(m_x1-xE)))
                              - log(abs(x+xE))*log(abs((m_x2+x)/(m_x2-xE)))
                              + log(abs(x+m_x1))*log(abs((m_x2+x)/(m_x2-m_x1)))
                              + DiLog((xE+x)/(xE-m_x1))
                              - DiLog((xE+x)/(xE-m_x2))
                              + DiLog((m_x1+x)/(m_x1-m_x2)));

  return (r1+r2+r3+r4+r5+r6);
}

// int dx G(x)/px²                                                        (C.88)
double YFS_Form_Factor::IntG() {
  // if dipole in its CMS
  if ((Vec3D(m_p1)+Vec3D(m_p2)).Abs() < 1E-3) {
    // same mass or both nearly massless or both of nearly same beta
    if ((abs(m_m1-m_m2) < 1E-6) || 
        ((1.-CalculateBeta(m_p1) < 5E-3) && (1.-CalculateBeta(m_p2) < 5E-3)) ||
        ((CalculateBeta(m_p1)-CalculateBeta(m_p2))
          /(CalculateBeta(m_p1)+CalculateBeta(m_p2)) < 5E-3)) {
      double E = m_p1[0];
      double b = CalculateBeta(m_p1);
      double r = 1./(b*E*E)*(1./2.*sqr(log((1.+b)/2.)) + LOG_2*log(1.+b)
                             - 1./2.*sqr(LOG_2) - 1./2.*sqr(log(1.+b))
                             + DiLog((1.-b)/2.) - DiLog((1.+b)/2.)
                             + DiLog(b) - DiLog(-b));
      return r;
    }
    // (p1-p2)^2=0 and m2 >> m1 (leptonic W-decay)
    else if ((abs((m_p1-m_p2).Abs2()) < 1E-6) &&
             (m_p1.Abs2()/m_p2.Abs2() < 1E-3))
      return 2./m_p2.Abs2()*(3./12.*M_PI*M_PI+DiLog(-2.));
  }
  // numerical calculation
// #define USING__Explicit_Check
#ifdef USING__Explicit_Check
  unsigned int n1   = 5000;
  unsigned int n2   = 5000;
  double       sum1 = 0, sum2 = 0;
  for (unsigned int i=0; i<n1; i++) {
    double x1 = -1 + (0.1*i)/n1;
    double x2 =  1 - (0.1*i)/n1;
    Vec4D px1 = (1./2.)*((m_p1+m_p2)+x1*(m_p1-m_p2));
    Vec4D px2 = (1./2.)*((m_p1+m_p2)+x2*(m_p1-m_p2));
    sum1 = sum1 + 0.1/n1*G(x1)/(px1*px1);
    sum1 = sum1 + 0.1/n1*G(x2)/(px2*px2);
  }
#endif
  double csum=p_gi1->Integrate(0.9,1.0,1.0e-4);
#ifdef USING__Explicit_Check
  for (unsigned int i=0; i<=n2; i++) {
    double x  = -1 + 0.1 + (1.8*i)/n2;
    Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
    sum2 = sum2 + 1.8/n2*G(x)/(px*px);
  }
#endif
  double ccsum=p_gi2->Integrate(-0.9,0.9,1.0e-4);
#ifdef USING__Explicit_Check
  msg_Debugging()<<"YFS FF: sum 1 = "<<sum1<<" vs. "<<csum
		 <<", rel. diff. = "<<sum1/csum-1.0<<"\n";
  msg_Debugging()<<"YFS FF: sum 2 = "<<sum2<<" vs. "<<ccsum
		 <<", rel. diff. = "<<sum2/ccsum-1.0<<"\n";
#endif
  return csum+ccsum;
}

double YFS_Form_Factor::CalculateBeta(const Vec4D& p) {
  return (Vec3D(p).Abs()/p[0]);
}

