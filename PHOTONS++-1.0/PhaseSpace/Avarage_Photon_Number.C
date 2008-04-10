#include "Avarage_Photon_Number.H"
#include "Model_Base.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Avarage_Photon_Number::Avarage_Photon_Number(Particle_Vector dip, double wmax, double wmin) {
  m_nbar     = 0;
  m_omegaMax = wmax;
  m_omegaMin = wmin;
  m_dipole   = dip;
  m_number1  = (int) 1E2;
  m_number2  = (int) 1E4;
  m_number3  = (int) 1E2;

  CalculateAvaragePhotonNumber();

#ifdef PHOTONS_DEBUG
  msg_Events()<<endl<<"nbar: "<<m_nbar<<" = ";
  for (unsigned int i=0; i<m_nbars.size(); i++) {
    msg_Events()<<m_nbars.at(i)<<" + ";
  }
  msg_Events()<<endl;
#endif
}

Avarage_Photon_Number::~Avarage_Photon_Number() {
}

void Avarage_Photon_Number::CalculateAvaragePhotonNumber() {
  double alpha   = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  double sum    = 0;
  for(unsigned int j=0; j<m_dipole.size(); j++) {
    for(unsigned int i=0; i<j; i++) {
      double Zi       = m_dipole.at(i)->Flav().Charge();
      double Zj       = m_dipole.at(j)->Flav().Charge();
      double titj     = 0;
      if ((m_dipole.at(i)->ProductionBlob() == m_dipole.at(j)->ProductionBlob()) && (m_dipole.at(i)->ProductionBlob() != NULL))       titj = +1;
      else if ((m_dipole.at(i)->DecayBlob() == m_dipole.at(j)->ProductionBlob()) && (m_dipole.at(i)->DecayBlob() != NULL))            titj = -1;
      else if ((m_dipole.at(i)->ProductionBlob() == m_dipole.at(j)->DecayBlob()) && (m_dipole.at(i)->ProductionBlob() != NULL))       titj = -1;
      else if ((m_dipole.at(i)->DecayBlob() == m_dipole.at(j)->DecayBlob()) && (m_dipole.at(i)->DecayBlob() != NULL))            titj = +1;
      double bi       = CalculateBeta(m_dipole.at(i)->Momentum());
      double bj       = CalculateBeta(m_dipole.at(j)->Momentum());
      Vec3D pi        = Vec3D(m_dipole.at(i)->Momentum());
      Vec3D pj        = Vec3D(m_dipole.at(j)->Momentum());
      double angle    = 1./2.*(M_PI - acos(pi*pj/(pi.Abs()*pj.Abs())));
      double integral = IntegrateInterferenceTerm(bi,bj,angle);
      double whole    = alpha/M_PI*Zi*Zj*titj*log(m_omegaMax/m_omegaMin)*(2 - 1./(2*M_PI)*(1-bi*bj*pow(sin(angle),2)+bi*bj*pow(cos(angle),2))*integral);
      m_nbars.push_back(whole);
      sum = sum + whole;
    }
  }
  m_nbar = sum;
}

double Avarage_Photon_Number::CalculateBeta(Vec4D p) {
  return Vec3D(p).Abs()/p[0];
}

double Avarage_Photon_Number::Function(double bi, double bj, double alpha, double phi) {
  double a1 = bi*sin(alpha)*sin(phi);
  double a2 = bj*sin(alpha)*sin(phi);
  double b1 = bi*cos(alpha);
  double b2 = bj*cos(alpha);
  double F = 2*(a1*(1+b1*b2)-a2*(1-b1*b1))/sqrt(1-a1*a1-b1*b1)*(M_PI/2+atan(a1/sqrt(1-a1*a1-b1*b1))) + 2*(a2*(1+b1*b2)-a1*(1-b2*b2))/sqrt(1-a2*a2-b2*b2)*(M_PI/2+atan(a2/sqrt(1-a2*a2-b2*b2)));
  double f = -(1-b1*b1)*a2*a2-(1-b2*b2)*a1*a1+2*(1+b1*b2)*a1*a2-(b1+b2)*(b1+b2);
  return (F/f);
}

double Avarage_Photon_Number::IntegralOne(double bi, double bj, double alpha) {
  // devide into peaked part (a->pi/2) and relatively flat part (0->a)
  double a = 1.5 - 1./50.*log(bi*bj);
  // integrate flat part by rectangle method
  double sumf = 0;
  for (unsigned int i=0; i<m_number1; i++) {
    double x = a/m_number1*i + a/(2.*m_number1);
    double F = Function(bi,bj,alpha,x)*(a/m_number1);
    sumf = sumf + F;
  }
  // integrate peaked part, for now with rectangle method
  double sump = 0;
  for (unsigned int i=0; i<m_number2; i++) {
    double x = a + (M_PI/2-a)/m_number2*i + (M_PI/2-a)/(2.*m_number2);
    double F = Function(bi,bj,alpha,x)*((M_PI/2-a)/m_number2);
    sump = sump + F;
  }
  return (sumf + sump);
}

double Avarage_Photon_Number::IntegralTwo(double bi, double bj, double alpha) {
  // integrate by rectangle method
  double sum = 0;
  for (unsigned int i=0; i<m_number3; i++) {
    double x = M_PI + (M_PI/2.*i)/m_number3 + M_PI/(4.*m_number3);
    double F = Function(bi,bj,alpha,x)*((M_PI/2)/m_number3);
    sum = sum + F;
  }
  return (sum);
}

double Avarage_Photon_Number::IntegrateInterferenceTerm(double bi, double bj, double alpha) {
  // calculates the integral dphi dtheta sin(theta)/((1-b1 sin(alpha)sin(theta)sin(phi)-b1
  // cos(alpha)cos(theta))(1-b2 sin(alpha)sin(theta)sin(phi)+b2 cos(alpha)cos(theta)))

  double term1 = -2*M_PI/((bi+bj)*cos(alpha))*log(((1-bi)*(1-bj))/((1+bi)*(1+bj)))
                  * 1/sqrt(1-(pow(2*bi*bj*cos(alpha),2)-pow(bi-bj,2))/(pow((bi+bj)*cos(alpha),2))*pow(sin(alpha),2));
  double term2 = -2*IntegralOne(bi,bj,alpha);
  double term3 = -2*IntegralTwo(bi,bj,alpha);

  return (term1+term2+term3);
}
