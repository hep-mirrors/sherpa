#include "PHOTONS++/Tools/Weight_Higher_Order_Corrections.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHOTONS;
using namespace std;

Weight_Higher_Order_Corrections::Weight_Higher_Order_Corrections(Particle_Vector_Vector pvv_old, Particle_Vector_Vector pvv_new, Dipole_Type::code dtype) {
  m_weight      = 1;
  m_maxweight   = 1;
  m_ME          = false;
//   check if exact matrix element exists
  ME_Browser ME(pvv_old);
  m_ME = ME.GetMEInfo();

  if ((m_ME == true) && (PHOTONS::Photons::s_useme == true)) {
    m_which = ME.WhichME();
    CalculateWeightAndMax(pvv_old,pvv_new);
  }
  else {
    m_dtype       = dtype;
    m_newdipole   = pvv_new.at(2);
    m_olddipole   = pvv_old.at(2);
    m_softphotons = pvv_new.at(4);
    if (m_dtype == Dipole_Type::ff)  m_M = pvv_old.at(1).at(0)->FinalMass();
    if (m_dtype == Dipole_Type::fi)  m_M = pvv_old.at(0).at(0)->FinalMass();

    CalculateWeight();
    CalculateMax();
  }
}

Weight_Higher_Order_Corrections::~Weight_Higher_Order_Corrections() {
}

double Weight_Higher_Order_Corrections::RealCorrectionsOrder(int order) {
  if (order == 0)  return 0;
  else if (order == 1) {
    double sum = 0;
    for (unsigned int i=0; i<m_softphotons.size(); i++) {
      double summ = 0;
      for (unsigned int j=0; j<m_newdipole.size(); j++) {
        for (unsigned int k=0; k<j; k++) {
          double Z1 = m_newdipole.at(k)->Flav().Charge();
          double Z2 = m_newdipole.at(j)->Flav().Charge();
          double t1t2 = 0;
          if (m_newdipole.at(j)->ProductionBlob() == m_newdipole.at(k)->ProductionBlob())
            t1t2 = +1;
          else if (m_newdipole.at(j)->DecayBlob() == m_newdipole.at(k)->ProductionBlob())
            t1t2 = -1;
          else if (m_newdipole.at(j)->ProductionBlob() == m_newdipole.at(k)->DecayBlob())
            t1t2 = -1;
          else if (m_newdipole.at(j)->DecayBlob() == m_newdipole.at(k)->DecayBlob())
            t1t2 = +1;
          else
            t1t2 = 0;
          summ = summ + Z1*Z2*t1t2*(Dmod(k,j,i)+Dmod(j,k,i));
        }
      }
      sum = sum + summ/Smod(i);
    }
    double r = -sum;
    return r;
  }
  else  return 0;
}

double Weight_Higher_Order_Corrections::VirtualCorrectionsOrder(int order) {
  double alpha   = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  if (order == 0)  return 0;
  else if (order == 1) {
    double r = 0;
    if (m_dtype == Dipole_Type::ff) {
      // only for decays Z -> fermion/antifermion
      if ((abs(m_M - 91.188) < 1E-6) && (m_newdipole.size() == 2) && (m_newdipole.at(0)->Flav().IsLepton() == true))
        r = alpha/M_PI * (2*log(m_M/m_newdipole.at(0)->FinalMass())+3./2.);
    }
    else if (m_dtype == Dipole_Type::fi) {
      // only for decays W -> fermion/antifermion
      if ((abs(m_M - 80.419) < 1E-6) && (m_newdipole.at(0)->Flav().Kfcode() == kf_Wplus) && (m_newdipole.size() == 2) && (m_newdipole.at(1)->Flav().IsLepton() == true))
        r = alpha/M_PI * (log(m_M/m_newdipole.at(1)->FinalMass())+1./2.);
    }
    return r;
  }
  else if (order == 2) {
    double r = 0;
    if (m_dtype == Dipole_Type::ff) {
      // only for decays Z -> fermion/antifermion
      if ((abs(m_M - 91.188) < 1E-6) && (m_newdipole.size() == 2) && (m_newdipole.at(0)->Flav().IsLepton() == true))
        r = 1./2.*pow(alpha/M_PI*2*log(m_M/m_newdipole.at(0)->FinalMass()),2);
    }
    return r;
  }
  else  return 0;
}

double Weight_Higher_Order_Corrections::Dmod(unsigned int i, unsigned int j, unsigned int kk) {
  Vec4D pi = m_newdipole.at(i)->Momentum();
  Vec4D pj = m_newdipole.at(j)->Momentum();
  Vec4D k  = m_softphotons.at(kk)->Momentum();
  double D = 1./(pi*k)*(2.*(pi*pj)/(pi*k+pj*k)-(pi*pi)/(pi*k));
  if ((m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->ProductionBlob()) && (m_newdipole.at(i)->ProductionBlob() != NULL )) {
    // emitter and spectator final state
      double y = (pi*k)/(pi*k+pj*k+pi*pj);
      double z = (pi*pj)/(pi*pj+pj*k);
      double p = (pi+pj+k)*(pi+pj+k);
      double P = p-pi*pi-pj*pj-k*k;
      double s = sqrt(pow(2.*(pj*pj)+P-P*y,2)-4.*p*(pj*pj));
      double R = s/sqrt(Kallen(p,pi*pi,pj*pj));
      if (m_newdipole.at(i)->Flav().IntSpin() == 0)
        return 0;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 1)
        return 1./((pi*k)*R)*(2./(1-z*(1-y))-1.-z-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 2)
        return 1./(pi*k) * (2./(1-z*(1-y))+2./(1-(1-z)*(1-y))+2*z*(1-z)-4-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 3)
        return 0;
  }
  else if ((m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->ProductionBlob()) && (m_newdipole.at(i)->DecayBlob() != NULL)) {
    // emitter is initial state, spectator is final state
      double x = (pi*k+pj*k-pi*pj)/(pi*pj+pj*k);
      double z = (pi*pj)/(pi*pj+pj*k);
      double p = (pj+k-pi)*(pj+k-pi);
      double P = p-pi*pi-pj*pj-k*k;
      double R = sqrt(pow(2.*(pj*pj)*x+P,2)-4.*p*x*x*(pj*pj))/sqrt(Kallen(p,pi*pi,pj*pj));
      if (m_newdipole.at(i)->Flav().IntSpin() == 0)
        return 0;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 1)
        return 1./((pi*k)*x)*(2./(2.-x-z)-R*(1+x)-x*(pi*pi)/(pi*k)) - D;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 2)
//         return 1./((pi*k)*x)*(2./(2-x-z)-2+2*x*(1-x)+2*(1-x)/x-2*(1-z)*(pj*pj)/(z*x*p)-x*(pi*pi)/(pi*k)) - D;
        return 0;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 3)
        return 0;
  }
  else if ((m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->DecayBlob()) && (m_newdipole.at(i)->ProductionBlob() != NULL )) {
    // emitter is final state, spectator is initial state
      double x = (pi*pj+pj*k-pi*k)/(pi*pj+pj*k);
      double z = (pi*pj)/(pi*pj+pj*k);
      if (m_newdipole.at(i)->Flav().IntSpin() == 0)
        return 0;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 1)
        return 1./((pi*k)*x)*(2./(2.-x-z)-1.-z-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 2)
        return 1./((pi*k)*x)*(2./(2.-x-z)+2./(2.-x-(1-z))+2*z*(1-z)-4-(pi*pi)/(pi*k)) - D;
      else if (m_newdipole.at(i)->Flav().IntSpin() == 3)
        return 0;
  }
  else if (m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->DecayBlob()) {
    // emitter is initial state, spectator is initial state
    return 0;
  }
  else
    return 0;
  return 0;
}

double Weight_Higher_Order_Corrections::Smod(unsigned int kk) {
  double sum    = 0;
  Vec4D k = m_softphotons.at(kk)->Momentum();
  for (unsigned int j=0; j<m_newdipole.size(); j++) {
    for (unsigned int i=0; i<j; i++) {
      Vec4D pi  = m_newdipole.at(i)->Momentum();
      Vec4D pj  = m_newdipole.at(j)->Momentum();
      double Z1 = m_newdipole.at(i)->Flav().Charge();
      double Z2 = m_newdipole.at(j)->Flav().Charge();
      double t1t2 = 0;
      if (m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->ProductionBlob())
        t1t2 = +1;
      else if (m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->ProductionBlob())
        t1t2 = -1;
      else if (m_newdipole.at(i)->ProductionBlob() == m_newdipole.at(j)->DecayBlob())
        t1t2 = -1;
      else if (m_newdipole.at(i)->DecayBlob() == m_newdipole.at(j)->DecayBlob())
        t1t2 = +1;
      else
        t1t2 = 0;
      sum = sum + Z1*Z2*t1t2*((pi/(pi*k)-pj/(pj*k))*(pi/(pi*k)-pj/(pj*k)));
    }
  }
  return sum;
}

double Weight_Higher_Order_Corrections::Kallen(double x, double y, double z) {
  return ( x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z );
}

void Weight_Higher_Order_Corrections::CalculateWeight() {
  m_weight = 1 + VirtualCorrectionsOrder(1) + RealCorrectionsOrder(1) + VirtualCorrectionsOrder(2);
}

void Weight_Higher_Order_Corrections::CalculateMax() {
  m_maxweight = 1 + VirtualCorrectionsOrder(1) + VirtualCorrectionsOrder(2);
}

PHOTONS_ME_Base * Weight_Higher_Order_Corrections::MESelector(Particle_Vector_Vector pvv_old, Particle_Vector_Vector pvv_new)
{
  if      (m_which == ME_List::Z_To_Lepton_Lepton) 
    return new  Z_To_Lepton_Lepton(pvv_old,pvv_new);
  else if (m_which == ME_List::W_To_Lepton_Neutrino)
    return new  W_To_Lepton_Neutrino(pvv_old,pvv_new);
  else if (m_which == ME_List::Vector_To_Fermion_Fermion)
    return new  Vector_To_Fermion_Fermion(pvv_old,pvv_new);
  else if (m_which == ME_List::Scalar_To_Fermion_Fermion)
    return new  Scalar_To_Fermion_Fermion(pvv_old,pvv_new);
  else if (m_which == ME_List::Tau_To_Lepton_Neutrinos)
    return new  Tau_To_Lepton_Neutrinos(pvv_old,pvv_new);
  else
    return NULL;
}

void Weight_Higher_Order_Corrections::CalculateWeightAndMax(Particle_Vector_Vector pvv_old, Particle_Vector_Vector pvv_new) {
  PHOTONS_ME_Base * me = MESelector(pvv_old,pvv_new);
  if (me == NULL) {
    m_weight    = 0;
    m_maxweight = 0;
    return;
  }
  else {
    double B_0_0 = me->GetBeta_0_0();
    double B_0_1 = me->GetBeta_0_1();
    double B_0_2 = me->GetBeta_0_2();

    double Sum_1_1 = 0;
    double Sum_1_2 = 0;
    for (unsigned int i=0; i<pvv_new.at(4).size(); i++) {
      Sum_1_1 = Sum_1_1 + me->GetBeta_1_1(i)/me->Smod(i);
      Sum_1_2 = Sum_1_2 + me->GetBeta_1_2(i)/me->Smod(i);
    }

    double Sum_2_2 = 0;
    for (unsigned int j=0; j<pvv_new.at(4).size(); j++) {
      for (unsigned int i=0; i<j; i++) {
        Sum_2_2 = Sum_2_2 + me->GetBeta_2_2(i,j)/(me->Smod(i)*me->Smod(i));
      }
    }

    m_weight    = 1 + 1./B_0_0*( B_0_1 + B_0_2 + Sum_1_1 + Sum_1_2 + Sum_2_2 );
    m_maxweight = 1 + 1./B_0_0*( B_0_1 + B_0_2 );
  }
  delete me;
}
