#include "PHOTONS++/Tools/Weight_Higher_Order_Corrections.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "MODEL/Main/Model_Base.H"
#include "PHOTONS++/MEs/PHOTONS_ME_Base.H"
#include "PHOTONS++/Main/Photons.H"

using namespace ATOOLS;
using namespace PHOTONS;
using namespace std;

Weight_Higher_Order_Corrections::Weight_Higher_Order_Corrections
(const Particle_Vector_Vector& pvv_old, const Particle_Vector_Vector& pvv_new,
 Dipole_Type::code dtype) : m_n(pvv_new[4].size()), p_pme(NULL) {
  DEBUG_FUNC(PHOTONS::Photons::s_useme);
  if (PHOTONS::Photons::s_useme)
    p_pme = PHOTONS_ME_Base::GetIRsubtractedME(pvv_old);
  if (p_pme) {
    msg_Debugging()<<"ME -> "<<p_pme->Name()<<std::endl;
    p_pme->FillMomentumArrays(pvv_new);
    CalculateWeightAndMaxWithME();
  }
  else {
    msg_Debugging()<<"ME -> none"<<std::endl;
    m_dtype       = dtype;
    m_newdipole   = pvv_new[2];
    m_olddipole   = pvv_old[2];
    m_softphotons = pvv_new[4];
    if (m_dtype == Dipole_Type::ff)  m_M = pvv_old[1][0]->FinalMass();
    if (m_dtype == Dipole_Type::fi)  m_M = pvv_old[0][0]->FinalMass();
    CalculateWeight();
    CalculateMax();
  }
}

Weight_Higher_Order_Corrections::~Weight_Higher_Order_Corrections() {
  if (p_pme) delete p_pme;
}

double Weight_Higher_Order_Corrections::RealCorrectionsOrder(int order) {
  if (order == 0)  return 0.;
  else if (order == 1) {
    msg_Debugging()<<"calc for "<<m_softphotons.size()<<" "
                   <<m_newdipole.size()<<std::endl;
    double sum = 0;
    for (unsigned int i=0; i<m_softphotons.size(); i++) {
      double summ(0.);
      for (unsigned int j=0; j<m_newdipole.size(); j++) {
        for (unsigned int k=0; k<j; k++) {
          double Z1(m_newdipole[k]->Flav().Charge());
          double Z2(m_newdipole[j]->Flav().Charge());
          double t1t2(0.);
          if (m_newdipole[j]->ProductionBlob()
                == m_newdipole[k]->ProductionBlob())
            t1t2 = +1.;
          else if (m_newdipole[j]->DecayBlob()
                     == m_newdipole[k]->ProductionBlob())
            t1t2 = -1.;
          else if (m_newdipole[j]->ProductionBlob()
                     == m_newdipole[k]->DecayBlob())
            t1t2 = -1.;
          else if (m_newdipole[j]->DecayBlob()
                     == m_newdipole[k]->DecayBlob())
            t1t2 = +1.;
          else
            t1t2 = 0.;
          summ = summ + Z1*Z2*t1t2*(Dmod(k,j,i)+Dmod(j,k,i));
        }
      }
      sum = sum + summ/Smod(i);
    }
    return -sum;
  }
  return 0.;
}

double Weight_Higher_Order_Corrections::VirtualCorrectionsOrder(int order) {
  return 0.;
}

double Weight_Higher_Order_Corrections::Dmod(unsigned int i, unsigned int j,
                                             unsigned int kk) {
  Vec4D pi = m_newdipole[i]->Momentum();
  Vec4D pj = m_newdipole[j]->Momentum();
  Vec4D k  = m_softphotons[kk]->Momentum();
  double D = 1./(pi*k)*(2.*(pi*pj)/(pi*k+pj*k)-(pi*pi)/(pi*k));
  if ((m_newdipole[i]->ProductionBlob() == m_newdipole[j]->ProductionBlob()) &&
      (m_newdipole[i]->ProductionBlob() != NULL )) {
    // emitter and spectator final state
      double y = (pi*k)/(pi*k+pj*k+pi*pj);
      double z = (pi*pj)/(pi*pj+pj*k);
      double p = (pi+pj+k)*(pi+pj+k);
      double P = p-pi*pi-pj*pj-k*k;
      double s = sqrt(sqr(2.*(pj*pj)+P-P*y)-4.*p*(pj*pj));
      double R = s/sqrt(Kallen(p,pi*pi,pj*pj));
      if (m_newdipole[i]->Flav().IntSpin() == 0)
        return 0.;
      else if (m_newdipole[i]->Flav().IntSpin() == 1)
        return 1./((pi*k)*R)*(2./(1.-z*(1.-y))-1.-z-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole[i]->Flav().IntSpin() == 2)
        return 1./(pi*k) * (2./(1.-z*(1.-y))+2./(1.-(1.-z)*(1.-y))+2.*z*(1.-z)
                            -4.-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole[i]->Flav().IntSpin() == 3)
        return 0.;
  }
  else if ((m_newdipole[i]->DecayBlob() == m_newdipole[j]->ProductionBlob()) &&
           (m_newdipole[i]->DecayBlob() != NULL)) {
    // emitter is initial state, spectator is final state
      double x = (pi*k+pj*k-pi*pj)/(pi*pj+pj*k);
      double z = (pi*pj)/(pi*pj+pj*k);
      double p = (pj+k-pi)*(pj+k-pi);
      double P = p-pi*pi-pj*pj-k*k;
      double R = sqrt(sqr(2.*(pj*pj)*x+P)-4.*p*x*x*(pj*pj))
                  /sqrt(Kallen(p,pi*pi,pj*pj));
      if (m_newdipole[i]->Flav().IntSpin() == 0)
        return 0.;
      else if (m_newdipole[i]->Flav().IntSpin() == 1)
        return 1./((pi*k)*x)*(2./(2.-x-z)-R*(1.+x)-x*(pi*pi)/(pi*k)) - D;
      else if (m_newdipole[i]->Flav().IntSpin() == 2)
//         return 1./((pi*k)*x)*(2./(2.-x-z)-2.+2.*x*(1.-x)+2.*(1.-x)/x
//                               -2.*(1.-z)*(pj*pj)/(z*x*p)-x*(pi*pi)/(pi*k)) - D;
        return 0.;
      else if (m_newdipole[i]->Flav().IntSpin() == 3)
        return 0.;
  }
  else if ((m_newdipole[i]->ProductionBlob() == m_newdipole[j]->DecayBlob()) &&
           (m_newdipole[i]->ProductionBlob() != NULL )) {
    // emitter is final state, spectator is initial state
      double x = (pi*pj+pj*k-pi*k)/(pi*pj+pj*k);
      double z = (pi*pj)/(pi*pj+pj*k);
      if (m_newdipole[i]->Flav().IntSpin() == 0)
        return 0.;
      else if (m_newdipole[i]->Flav().IntSpin() == 1)
        return 1./((pi*k)*x)*(2./(2.-x-z)-1.-z-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole[i]->Flav().IntSpin() == 2)
        return 1./((pi*k)*x)*(2./(2.-x-z)+2./(2.-x-(1.-z))+2.*z*(1.-z)-4.
                              -(pi*pi)/(pi*k)) - D;
      else if (m_newdipole[i]->Flav().IntSpin() == 3)
        return 0.;
  }
  else if ((m_newdipole[i]->DecayBlob() == m_newdipole[j]->DecayBlob()) &&
           (m_newdipole[i]->ProductionBlob() != NULL )) {
    // emitter is initial state, spectator is initial state
    return 0.;
  }
  else
    return 0.;
  return 0.;
}

double Weight_Higher_Order_Corrections::Smod(unsigned int kk) {
  double sum(0.);
  Vec4D k = m_softphotons[kk]->Momentum();
  for (unsigned int j=0; j<m_newdipole.size(); j++) {
    for (unsigned int i=0; i<j; i++) {
      Vec4D pi  = m_newdipole[i]->Momentum();
      Vec4D pj  = m_newdipole[j]->Momentum();
      double Z1 = m_newdipole[i]->Flav().Charge();
      double Z2 = m_newdipole[j]->Flav().Charge();
      double t1t2 = 0.;
      if (m_newdipole[i]->ProductionBlob() == m_newdipole[j]->ProductionBlob())
        t1t2 = +1.;
      else if (m_newdipole[i]->DecayBlob() == m_newdipole[j]->ProductionBlob())
        t1t2 = -1.;
      else if (m_newdipole[i]->ProductionBlob() == m_newdipole[j]->DecayBlob())
        t1t2 = -1.;
      else if (m_newdipole[i]->DecayBlob() == m_newdipole[j]->DecayBlob())
        t1t2 = +1.;
      else
        t1t2 = 0.;
      sum = sum + Z1*Z2*t1t2*((pi/(pi*k)-pj/(pj*k))*(pi/(pi*k)-pj/(pj*k)));
    }
  }
  return sum;
}

double Weight_Higher_Order_Corrections::Kallen(double x, double y, double z) {
  return ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z );
}

void Weight_Higher_Order_Corrections::CalculateWeight() {
  m_weight = 1. + VirtualCorrectionsOrder(1) + RealCorrectionsOrder(1)
                + VirtualCorrectionsOrder(2);
}

void Weight_Higher_Order_Corrections::CalculateMax() {
  m_maxweight = 1. + VirtualCorrectionsOrder(1)
                   + VirtualCorrectionsOrder(2);
}

void Weight_Higher_Order_Corrections::CalculateWeightAndMaxWithME() {
  DEBUG_FUNC("");
  double B_0_0 = p_pme->GetBeta_0_0();
  double B_0_1 = p_pme->GetBeta_0_1();
  double B_0_2 = p_pme->GetBeta_0_2();

  double Sum_1_1(0.), Sum_1_2(0.);
  for (unsigned int i=0; i<m_n; i++) {
    Sum_1_1 = Sum_1_1 + p_pme->GetBeta_1_1(i)/p_pme->Smod(i);
    Sum_1_2 = Sum_1_2 + p_pme->GetBeta_1_2(i)/p_pme->Smod(i);
  }

  double Sum_2_2(0.);
  for (unsigned int j=0; j<m_n; j++) {
    for (unsigned int i=0; i<j; i++) {
      Sum_2_2 = Sum_2_2 + p_pme->GetBeta_2_2(i,j)
                          /(p_pme->Smod(i)*p_pme->Smod(j));
    }
  }
  if (msg_LevelIsDebugging()) {
    std::string name(p_pme->Name());
    std::string space(name.size(),' ');
    msg_Out()<<name <<"  \\beta_0^0 = "<<B_0_0<<std::endl
             <<space<<"  \\beta_0^1 = "<<B_0_1<<std::endl
             <<space<<"  \\beta_0^2 = "<<B_0_2<<std::endl
             <<space<<"  \\sum_i[\\beta_1^1(i)/S(i)] = "<<Sum_1_1<<std::endl
             <<space<<"  \\sum_i[\\beta_1^2(i)/S(i)] = "<<Sum_1_2<<std::endl
             <<space<<"  \\sum_i\\sum_j[\\beta_2^2(i,j)/S(i)S(j)] = "<<Sum_2_2
             <<std::endl;
  }

  m_weight    = 1. + 1./B_0_0*( B_0_1 + B_0_2 + Sum_1_1 + Sum_1_2 + Sum_2_2 );
  m_maxweight = 1. + 1./B_0_0*( B_0_1 + B_0_2 );
}

