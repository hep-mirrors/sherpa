#include "YFS_Form_Factor.H"
#include "Model_Base.H"
#include "Message.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

YFS_Form_Factor::YFS_Form_Factor(Particle_Vector part, double ks) {
  double YSum = 0;
  for (unsigned int j=0; j<part.size(); j++) {
    for (unsigned int i=0; i<j; i++) {
      YFS_Form_Factor YFS(part.at(i),part.at(j),ks);
      YSum = YSum + YFS.Get();
    }
  }
  m_Y = YSum;
}

YFS_Form_Factor::YFS_Form_Factor(Particle * part1, Particle * part2, double ks) {
  m_ks = ks;

  // choose such that E_2 >= E_1
  if (part2->Momentum()[0] >= part1->Momentum()[0]) {
    p_part1 = new Particle(*part1);
    p_part2 = new Particle(*part2);
  }
  else {
    p_part1 = new Particle(*part2);
    p_part2 = new Particle(*part1);
  }

  m_p1 = p_part1->Momentum();
  m_p2 = p_part2->Momentum();

  m_m1 = p_part1->FinalMass();
  m_m2 = p_part2->FinalMass();

  m_Z1 = p_part1->Flav().Charge();
  m_Z2 = p_part2->Flav().Charge();

  if (part1->ProductionBlob() == part2->ProductionBlob())         m_t1t2 = +1;
  else if (part1->ProductionBlob() == part2->DecayBlob())         m_t1t2 = -1;
  else if (part1->DecayBlob() == part2->ProductionBlob())         m_t1t2 = -1;
  else if (part1->DecayBlob() == part2->DecayBlob())              m_t1t2 = +1;
  else                                                            m_t1t2 = 0;

  //roots of p_x^2
  m_x1  = - (m_p1.Abs2() - m_p2.Abs2() + 2*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                / (m_p1-m_p2).Abs2();
  m_x2  = - (m_p1.Abs2() - m_p2.Abs2() - 2*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                / (m_p1-m_p2).Abs2();

  //roots of p_x^'2
  if (m_t1t2 == +1){
    m_xx1 = - (m_p1.Abs2() - m_p2.Abs2() + 2*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                / (m_p1+m_p2).Abs2();
    m_xx2 = - (m_p1.Abs2() - m_p2.Abs2() - 2*sqrt((m_p1*m_p2)*(m_p1*m_p2)-(m_p1.Abs2()*m_p2.Abs2())))
                / (m_p1+m_p2).Abs2();
  }

  m_Y = Y();
  delete p_part1, delete p_part2;
}


YFS_Form_Factor::~YFS_Form_Factor() {
}

// private members

// Y = 2 alpha (Re B + B~)
double YFS_Form_Factor::Y() {
  double alpha = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  return (-alpha/M_PI*m_Z1*m_Z2*m_t1t2 * (log((m_p1[0]*m_p2[0])/(m_ks*m_ks))
          + (1./2.)*(m_p1*m_p2)*IntP1() - (1./2.)*(m_p1*m_p2)*IntE() + (1./4.)*IntP2()
          + G(1) + G(-1) - (m_p1*m_p2)*IntG()));
}

// t1t2 * int dx ln (px'²/lambda²)/px'² + int dx ln (px²/lambda²)/px²
double YFS_Form_Factor::IntP1() {
  if (m_t1t2 == -1) {
    return 0;
    }
  else if (m_t1t2 == +1) {
    double A;
    if (m_xx1*m_xx2 >= 0) 
      A = 8*(M_PI*M_PI)/((m_xx2-m_xx1)*(m_p1+m_p2).Abs2());
    else
      A = 0;

    double B = 8/((m_p1-m_p2).Abs2()*(m_x1-m_x2)) * 
                                      (  log(abs(m_x1))*(DiLog((m_x1-1)/m_x1)-DiLog((m_x1+1)/m_x1))
                                       - log(abs(m_x2))*(DiLog((m_x2-1)/m_x2)-DiLog((m_x2+1)/m_x2)));
    return (A+B);
  }
  else
    return 0;
}

// int dx ln(Ex²/omega²)/px²
double YFS_Form_Factor::IntE() {
  if (abs(m_p1[0]-m_p2[0]) < 1E-6) {
    double r = (8/((m_p1-m_p2).Abs2()*(m_x1-m_x2)))*log((m_p1[0]+m_p2[0])/(2*m_ks))*log(abs(((1-m_x1)*(1+m_x2))/((1+m_x1)*(1-m_x2))));
    return r;
    }
  else {
    if ((m_p1-m_p2).Abs2() < -1E-6) {
      double xE   = -(m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
      double zeta = -(m_p1[0]-m_p2[0])/(2*m_p1[0]);
      double y1   = 1+zeta*(1-m_x1);
      double r = (8/((m_p1-m_p2).Abs2()*(m_x1-m_x2)) * (log(m_p1[0]/m_ks)*log(abs((1-m_x1)/(1+m_x1)))
                                                 + log(abs(y1))*log(abs((1-y1)/(1-y1+2*zeta)))
                                                 + DiLog(1-((1+2*zeta)/y1)) - DiLog(1-(1/y1))
                                                 - DiLog(-(1+m_x2)/(xE-m_x2))
                                                 + DiLog((1-m_x2)/(xE-m_x2))
                                                 - log((1/(2*m_ks))*((1+m_x2)*m_p1[0]+(1-m_x2)*m_p2[0]))*log(abs((1-m_x2)/(1+m_x2)))));
      return r;
    }
    else if ((m_p1-m_p2).Abs2() > 1E-6) {
      if ((m_m2*m_m2) > (2*(m_p1*m_p2)-m_m1*m_m1)) {
        double xi   = (m_p1[0]-m_p2[0])/(2*m_p2[0]);
        double zeta = -(m_p1[0]-m_p2[0])/(2*m_p1[0]);
        double y1   = 1+zeta*(1-m_x1);
        double y2   = 1+xi*(1+m_x2);
        double r = (8/((m_p1-m_p2).Abs2()*(m_x1-m_x2)) * (log(m_p1[0]/m_ks)*log(abs((1-m_x1)/(1+m_x1)))
                                       + log(abs(y1))*log(abs((1-y1)/(1-y1+2*zeta)))
                                       + DiLog(1-((1+2*zeta)/y1)) - DiLog(1-(1/y1))
                                       - DiLog((y2)/(y2-1-2*xi)) + DiLog((y2)/(y2-1))
                                       + (1./2.)*(pow(log(y2/(y2-1)),2)-pow(log(y2/(y2-1-2*xi)),2))
                                       - log(abs(y2))*log(abs((y2+1+2*xi)/(y2+1)))
                                       - log(m_p2[0]/m_ks)*log(abs((1-m_x2)/(1+m_x2)))));
        return r;
      }
      else if ((m_m1*m_m1) > (2*(m_p1*m_p2)-m_m2*m_m2)) {
        msg_Out()<<"!!! error: case should not appear because of ordering of particles !!!"<<endl;
        return 0;
      }
      else {
        msg_Out()<<"!!! error: case should not appear !!!"<<endl;
        return 0;
      }
    }
    else {
      // (m_p1-m_p2).Abs2() = 0
      if (abs(m_m1-m_m2) < 1E-6) {
        if (abs(m_p1[0]-m_p2[0]) > 1E-6) {
          double xE = - (m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
          double r = 8/(m_p1+m_p2).Abs2() 
                      *(2*log(1./2.*(m_p2[0]-m_p1[0])/m_ks)-xE*log((xE-1)/(xE+1))-log(xE*xE-1)-2);
          return r;
        }
        else {
          msg_Out()<<"!!! error: case should not appear !!!"<<endl;
          return 0;
        }
      }
      else {
        if (abs(m_p1[0]-m_p2[0]) > 1E-6) {
          double xE = - (m_p1[0]+m_p2[0])/(m_p1[0]-m_p2[0]);
          double xp = - (m_m1*m_m1+m_m2*m_m2)/(m_m1*m_m1-m_m2*m_m2);
          double r;
          if (abs(xE-xp) < 1E-6) {
            r = 4/(m_m2*m_m2-m_m1*m_m1)*(log((m_p2[0]-m_p1[0])/(2*m_ks))*log(abs((xp+1)/(xp-1)))
                                         - (1./2.)*(pow(log(xp-1),2)-pow(log(xp+1),2)));
          }
          else if (xE > xp) {
            r = 4/(m_m2*m_m2-m_m1*m_m1)*(log((m_p2[0]-m_p1[0])/(2*m_ks))*log(abs((xp+1)/(xp-1)))
                                         + log(xE-xp)*log(abs((xp+1)/(xp-1)))
                                         + DiLog((xp-1)/(xp-xE)) - DiLog((xp+1)/(xp-xE)));
          }
          else if (xp > xE) {
            double xi = (m_p1[0]-m_p2[0])/(2*m_p2[0]);
            double yp = 1+xi*(1+xp);
            r = 4/(m_m2*m_m2-m_m1*m_m1)*(log(m_p2[0]/m_ks)*log(abs((xp+1)/(xp-1)))
                                         + log(abs(yp))*log(abs((xp+1)/(xp-1)))
                                         - (1./2.)*pow(log(abs(yp/(yp-1-2*xi))),2)
                                         + (1./2.)*pow(log(abs(yp/(yp-1))),2)
                                         + DiLog(yp/(yp-1)) - DiLog(yp/(yp-1-2*xi)));
          }
          else {
            msg_Out()<<"!!! error: case should not appear !!!"<<endl;
            return 0;
          }
          return r;
        }
        else {
          msg_Out()<<"!!! error: case should not appear !!!"<<endl;
          return 0;
        }
      }
    }
  }
}

// int dx ln (px²/m1m2)
double YFS_Form_Factor::IntP2() {
  if(m_t1t2 == +1) {
    double r = (2*log((m_p1+m_p2).Abs2()/(4*m_m1*m_m2)) + log(abs((1-m_xx1*m_xx1)*(1-m_xx2*m_xx2)))
          - m_xx1*log(abs((1-m_xx1)/(1+m_xx1))) - m_xx2*log(abs((1-m_xx2)/(1+m_xx2))) - 4);
    return r;
  }
  else if (m_t1t2 == -1) {
    if (abs((m_p1-m_p2).Abs2()) > 1E-6) {
      double r = (2.*log(abs((m_p1-m_p2).Abs2())/(4*m_m1*m_m2)) + log(abs((1-m_x1*m_x1)*(1-m_x2*m_x2)))
                - m_x1*log(abs((1-m_x1)/(1+m_x1))) - m_x2*log(abs((1-m_x2)/(1+m_x2))) - 4);
      return r;
    }
    else if (abs((m_p1-m_p2).Abs2()) < 1E-6) {
      if (abs(m_m1*m_m1 - m_m2*m_m2) > 1E-6) {
        double xp = -(m_m1*m_m1+m_m2*m_m2)/(m_m1*m_m1 - m_m2*m_m2);
        double r = 2.*log(abs(m_m1*m_m1-m_m2*m_m2)/(2*m_m1*m_m2)) + log(abs(1-xp*xp)) + xp*log(abs((1+xp)/(1-xp))) - 2;
        return r;
      }
      else {
        double r = 2*log((m_p1+m_p2).Abs2()/(4*m_m1*m_m2));;
        return r;
      }
    }
    else {
      msg_Out()<<"!!! error: case should not appear !!!"<<endl;
      return 0;
    }
  }
  else {
    msg_Out()<<"!!! error:case should not appear !!!"<<endl;
    return 0;
  }
}

// G(x)
double YFS_Form_Factor::G(double x) {
  Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
  double b  = CalculateBeta(px);
  double r;
  if (b == 0)       r = 1.-log(2.);
  else if (b == 1)  r = 0.;
  else              r = ((1-b)/(2*b)*log((1+b)/(1-b))+log((1+b)/2.));
  return r;
}

// Function for evaluating IntG() for dipole of different masses in its rest frame
double YFS_Form_Factor::GFunc(double x) {
  double xE = (m_p2[0]+m_p1[0])/(m_p2[0]-m_p1[0]);
  double xD = 2*Vec3D(m_p1).Abs()/(m_p2[0]-m_p1[0]);

  double r1 = xE/(2*xD*m_x1*(m_x1-m_x2))*(-1./2.*pow(log(abs(x-m_x1)),2)
                                          + log(abs(m_x1/m_x2))*log(abs(x-m_x1))
                                          + log(abs(x-m_x2))*log(abs((x-m_x1)/(m_x2-m_x1)))
                                          + DiLog((m_x2-x)/(m_x2-m_x1)))
            + xE/(2*xD*m_x2*(m_x1-m_x2))*(-1./2.*pow(log(abs(x-m_x2)),2)
                                          - log(abs(m_x1/m_x2))*log(abs(x-m_x2))
                                          + log(abs(x-m_x1))*log(abs((x-m_x2)/(m_x1-m_x2)))
                                          + DiLog((m_x1-x)/(m_x1-m_x2)))
            + xE/(2*xD*m_x1*m_x2)*(pow(log(abs(m_x1)),2)
                                   - pow(log(abs(m_x2)),2)
                                   + DiLog(x/m_x1)
                                   - DiLog(x/m_x2));
  double r2 = -(1+xD)/(2*xD)*(-pow(log(abs(m_x1/m_x2*(m_x2-x)/(m_x1-x))),2)/(2*(m_x1-m_x2)));
  double r3 = 1./(m_x1-m_x2)*(1./2.*pow(log(abs(x-m_x2)),2)
                              + log(abs((1-xD)/2.))*log(abs((x-m_x1)/(x-m_x2)))
                              - log(abs(x-xE))*log(abs((m_x1-x)/(m_x1-xE)))
                              + log(abs(x-xE))*log(abs((m_x2-x)/(m_x2-xE)))
                              + log(abs(x-m_x2))*log(abs((m_x1-x)/(m_x1-m_x2)))
                              - DiLog((x-xE)/(m_x1-xE))
                              + DiLog((x-xE)/(m_x2-xE))
                              + DiLog((x-m_x2)/(m_x1-m_x2)));
  double r4 = xE/(2*xD*m_x1*(m_x1-m_x2))*(1./2.*pow(log(abs(x+m_x1)),2)
                                          + log(abs(m_x2/m_x1))*log(abs(x+m_x1))
                                          - log(abs(x+m_x2))*log(abs((m_x1+x)/(m_x1-m_x2)))
                                          - DiLog((m_x2+x)/(m_x2-m_x1)))
            + xE/(2*xD*m_x2*(m_x1-m_x2))*(1./2.*pow(log(abs(x+m_x2)),2)
                                          - log(abs(m_x2/m_x1))*log(abs(x+m_x2))
                                          - log(abs(x+m_x1))*log(abs((m_x2+x)/(m_x2-m_x1)))
                                          - DiLog((m_x1+x)/(m_x1-m_x2)))
            + xE/(2*xD*m_x1*m_x2)*(DiLog(-x/m_x1)-DiLog(-x/m_x2));
  double r5 = (1-xD)/(2*xD)*(-pow(log(abs(m_x2/m_x1*(m_x1+x)/(m_x2+x))),2)/(2*(m_x1-m_x2)));
  double r6 = 1./(m_x1-m_x2)*(-1./2.*pow(log(abs(x+m_x1)),2)
                              - log(abs((1-xD)/2))*log(abs((x+m_x1)/(x+m_x2)))
                              + log(abs(x+xE))*log(abs((m_x1+x)/(m_x1-xE)))
                              - log(abs(x+xE))*log(abs((m_x2+x)/(m_x2-xE)))
                              + log(abs(x+m_x1))*log(abs((m_x2+x)/(m_x2-m_x1)))
                              + DiLog((xE+x)/(xE-m_x1))
                              - DiLog((xE+x)/(xE-m_x2))
                              + DiLog((m_x1+x)/(m_x1-m_x2)));

  return (r1+r2+r3+r4+r5+r6);
}

// int dx G(x)/px²
double YFS_Form_Factor::IntG() {
  // if dipole in its CMS
  if ((Vec3D(m_p1)+Vec3D(m_p2)).Abs() < 1E-3) {
    // same mass or both nearly massless or both of nearly same beta
    if ((abs(m_m1-m_m2) < 1E-6) || 
        ((1-CalculateBeta(m_p1) < 5E-3) && (1- CalculateBeta(m_p2) < 5E-3)) ||
        ((CalculateBeta(m_p1)-CalculateBeta(m_p2))/(CalculateBeta(m_p1)+CalculateBeta(m_p2)) < 5E-3)) {
      double E = m_p1[0];
      double b = CalculateBeta(m_p1);
      double r = 1./(b*E*E)*(1./2.*pow(log((1+b)/2.),2) + log(2.)*log(1+b)
                             - 1./2.*pow(log(2.),2) - 1./2.*pow(log(1+b),2)
                             + DiLog((1-b)/2.) - DiLog((1+b)/2.)
                             + DiLog(b) - DiLog(-b));
      return r;
    }
    // (p1-p2)^2=0 and m2 >> m1 (leptonic W-decay)
    else if ((abs((m_p1-m_p2).Abs2()) < 1E-6) && (m_p1.Abs2()/m_p2.Abs2() < 1E-3))  {
      double r = 2./m_p2.Abs2()*(3./12.*M_PI*M_PI+DiLog(-2));
      return r;
    }
//     else {
//       double r = 4/(m_p1-m_p2).Abs2()*(GFunc(1) - GFunc(0));
//       msg_Events()<<-GFunc(1)<<" "<<GFunc(0)<<" "<<4/(m_p1-m_p2).Abs2()<<endl;
//       msg_Events()<<"IntG(m1!=m2,rest frame): "<<-(m_p1*m_p2)*r<<endl;
// //       return r;
//     }
  }

  unsigned int n1   = 500;
  unsigned int n2   = 500;
  double       sum  = 0;
  for (unsigned int i=0; i<n1; i++) {
    double x1 = -1 + (0.1*i)/n1;
    double x2 =  1 - (0.1*i)/n1;
    Vec4D px1 = (1./2.)*((m_p1+m_p2)+x1*(m_p1-m_p2));
    Vec4D px2 = (1./2.)*((m_p1+m_p2)+x2*(m_p1-m_p2));
    sum = sum + 0.1/n1*G(x1)/(px1*px1);
    sum = sum + 0.1/n1*G(x2)/(px2*px2);
  }
  for (unsigned int i=0; i<=n2; i++) {
    double x  = -1 + 0.1 + (1.8*i)/n2;
    Vec4D  px = (1./2.)*((m_p1+m_p2)+x*(m_p1-m_p2));
    sum = sum + 1.8/n2*G(x)/(px*px);
  }
  return (sum);
}

double YFS_Form_Factor::CalculateBeta(Vec4D p) {
  return (Vec3D(p).Abs()/p[0]);
}

