#include "PHOTONS++/MEs/Tau_To_Lepton_Neutrinos.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////                                                                  ////
////     CAUTION :                                                    ////
////                                                                  ////
////     pvv_zero contains m_chargedoutparticles at position 2        ////
////     --> final state lepton sits at position 0                    ////
////                                                                  ////
////     pvv_one contains m_newdipole at position 2                   ////
////     --> tau sits by default at m_newdipole.at(0)                 ////
////     --> final state lepton sits on postion 1                     ////
////                                                                  ////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

Tau_To_Lepton_Neutrinos::Tau_To_Lepton_Neutrinos(Particle_Vector_Vector pvv_zero, Particle_Vector_Vector pvv_one) : Dipole_FI(pvv_zero) {
  m_flavs[0] = (pvv_zero.at(0)).at(0)->Flav();    // tau
  m_flavs[1] = (pvv_zero.at(2)).at(0)->Flav();    // lepton
  // switch ordering if necessary
  m_switch = ((pvv_zero.at(3)).at(0)->Flav().Kfcode() == kf_nutau);
  // m_switch == true if first neutral final state particle is tau-neutrino
  if (m_switch == false) {
    m_flavs[2] = (pvv_zero.at(3)).at(0)->Flav();  // lepton-neutrino
    m_flavs[3] = (pvv_zero.at(3)).at(1)->Flav();  // tau-neutrino
  }
  else {
    m_flavs[2] = (pvv_zero.at(3)).at(1)->Flav();  // lepton-neutrino
    m_flavs[3] = (pvv_zero.at(3)).at(0)->Flav();  // tau-neutrino
  }
  for (unsigned int i=4; i<9; i++) {
    m_flavs[i] = Flavour(kf_photon);
  }

  m_cL = m_i*m_e/(m_sW*sqrt(2.));
  m_cR = 0.;

  m_pvv_zero = pvv_zero;
  m_pvv_one  = pvv_one;

  // operate in Fermi-Theory
  m_fermi = true;

  FillMomentumArrays();
}

Tau_To_Lepton_Neutrinos::~Tau_To_Lepton_Neutrinos() {
}

void Tau_To_Lepton_Neutrinos::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into -z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum = sum + m_olddipole.at(i)->Momentum();
  }
  Vec4D p1 = m_olddipole.at(0)->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,-1.));
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D vec = m_olddipole.at(i)->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_olddipole.at(i)->SetMomentum(vec);
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D vec = m_oldspectator.at(i)->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_oldspectator.at(i)->SetMomentum(vec);
  }
}

void Tau_To_Lepton_Neutrinos::FillMomentumArrays() {
  // m_moms0 - no photon
  Poincare boost(m_pvv_zero.at(0).at(0)->Momentum());
  Vec4D vec;
  // tau
  vec = m_pvv_zero.at(0).at(0)->Momentum();
  boost.Boost(vec);
  m_moms0[0] = vec;
  // lepton
  vec = m_pvv_zero.at(2).at(0)->Momentum();
  boost.Boost(vec);
  m_moms0[1] = vec;
  if (m_switch == false) {
    // lepton-neutrino
    vec = m_pvv_zero.at(3).at(0)->Momentum();
    boost.Boost(vec);
    m_moms0[2] = vec;
    // tau-neutrino
    vec = m_pvv_zero.at(3).at(1)->Momentum();
    boost.Boost(vec);
    m_moms0[3] = vec;
  }
  else {
    // lepton-neutrino
    vec = m_pvv_zero.at(3).at(1)->Momentum();
    boost.Boost(vec);
    m_moms0[2] = vec;
    // tau-neutrino
    vec = m_pvv_zero.at(3).at(0)->Momentum();
    boost.Boost(vec);
    m_moms0[3] = vec;
  }

  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated

  // not necessary if only one photon
  if (m_pvv_one.at(4).size() == 1) {
    m_moms1[0][0] = m_pvv_one.at(2).at(0)->Momentum();
    m_moms1[0][1] = m_pvv_one.at(2).at(1)->Momentum();
    if (m_switch == false) {
      m_moms1[0][2] = m_pvv_one.at(3).at(0)->Momentum();
      m_moms1[0][3] = m_pvv_one.at(3).at(1)->Momentum();
    }
    else {
      m_moms1[0][2] = m_pvv_one.at(3).at(1)->Momentum();
      m_moms1[0][3] = m_pvv_one.at(3).at(0)->Momentum();
    }
    m_moms1[0][4] = m_pvv_one.at(4).at(0)->Momentum();
  }
  else {
    Dipole_FI::DefineDipole();
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<m_pvv_one.at(4).size(); i++) {
      m_softphotons.push_back(m_pvv_one.at(4).at(i));
      m_K = CalculateMomentumSum(m_softphotons);
      CorrectMomenta();
      m_moms1[i][0] = m_newdipole.at(0)->Momentum();
      m_moms1[i][1] = m_newdipole.at(1)->Momentum();
      if (m_switch == true) {
        m_moms1[i][2] = m_newspectator.at(0)->Momentum();
        m_moms1[i][3] = m_newspectator.at(1)->Momentum();
      }
      else {
        m_moms1[i][2] = m_newspectator.at(1)->Momentum();
        m_moms1[i][3] = m_newspectator.at(0)->Momentum();
      }
      m_moms1[i][4] = m_softphotons.at(0)->Momentum();
      m_softphotons.clear();
    }
  }
}

double Tau_To_Lepton_Neutrinos::Smod(unsigned int kk) {
  m_moms = m_moms1[kk];
  double sum = 0;
  Vec4D k   = m_moms[4];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = -1;   // both fermions of same charge for both tau+- decays
  double Zj = -1;
  int    ti = -1;   // always one initial state, one final state
  int    tj = +1;
  sum = m_alpha/(4*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k))*(pi/(pi*k)-pj/(pj*k));
  return sum;
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4D  q  = m_moms[0]-m_moms[3];          // W propagator
  double mW = Flavour(kf_Wplus).HadMass();
  XYZFunc XYZ(4,m_moms,m_flavs,1,false);
  // Fermi-Theory
  if (m_fermi == true)
    return m_i/(mW*mW)*
             XYZ.Z(3,m_spins[3],0,m_spins[0],1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
  // full SM
  else
    return  -m_i/(q*q-mW*mW)*
            (XYZ.Z(3,m_spins[3],0,m_spins[0],1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR)
              -1./(mW*mW)*XYZ.X(3,m_spins[3],q,0,m_spins[0],m_cL,m_cR)
                          *XYZ.X(1,m_spins[1],q,2,m_spins[2],m_cL,m_cR));
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_0_1() {
  return 0;
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_0_2() {
  return 0;
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsP   = conj(Polarization_Vector(m_moms[4]).at(m_spins[4]));
  Vec4D qt     = m_moms[0]-m_moms[4];       // tau propagator momenta
  Vec4D ql     = m_moms[1]+m_moms[4];       // lepton propagator momenta
  double mt    = m_flavs[0].HadMass();       // tau mass/propagator pole
  double ml    = m_flavs[1].HadMass();       // lepton mass/propagator pole
  m_moms[5]    = m_moms[6] = qt;            // enter those into m_moms
  m_moms[7]    = m_moms[8] = ql;
  m_flavs[5]   = m_flavs[0];                // set to corresponding particle/antiparticle
  m_flavs[6]   = m_flavs[0].Bar();
  m_flavs[7]   = m_flavs[1];
  m_flavs[8]   = m_flavs[1].Bar();
  Vec4D q1     = m_moms[1]+m_moms[2];       // W propagators
  Vec4D q2     = m_moms[0]-m_moms[3];
  Vec4D k      = m_moms[4];
  double mW    = Flavour(kf_Wplus).HadMass();   // W mass/propagator pole
  XYZFunc XYZ(9,m_moms,m_flavs,1,false);
  m_flavs[5] = m_flavs[6] = m_flavs[7] = m_flavs[8] = Flavour(kf_none);
  // Fermi-Theory
  if (m_fermi == true) {
    Complex r1 = Complex(0.,0.);
    Complex r2 = Complex(0.,0.);
    Complex r3 = Complex(0.,0.);
    Complex r4 = Complex(0.,0.);
    for (unsigned int s=0; s<=1; s++) { // spin of pseudo-particle in propagator representation
      r1 = r1 + XYZ.X(5,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.Z(3,m_spins[3],5,s,1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
      r2 = r2 + XYZ.X(6,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.Z(3,m_spins[3],6,s,1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
      r3 = r3 + XYZ.Z(3,m_spins[3],0,m_spins[0],7,s,2,m_spins[2],m_cL,m_cR,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,7,s,1.,1.);
      r3 = r3 + XYZ.Z(3,m_spins[3],0,m_spins[0],8,s,2,m_spins[2],m_cL,m_cR,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,8,s,1.,1.);
    }
    // add prefactors to those terms
    r1 = 0.5*(1+mt/sqrt(qt*qt))*r1/(qt*qt-mt*mt);
    r2 = 0.5*(1-mt/sqrt(qt*qt))*r2/(qt*qt-mt*mt);
    r3 = 0.5*(1+ml/sqrt(ql*ql))*r3/(ql*ql-ml*ml);
    r4 = 0.5*(1-ml/sqrt(ql*ql))*r4/(ql*ql-ml*ml);
    return m_i*m_e/(mW*mW)*(r1+r2+r3+r4);
  }
  // full SM
  else {
    // calculate terms arising from emission off the charged leptons
    Complex r1 = Complex(0.,0.);
    Complex r2 = Complex(0.,0.);
    Complex r3 = Complex(0.,0.);
    Complex r4 = Complex(0.,0.);
    Complex r5 = Complex(0.,0.);
    Complex r6 = Complex(0.,0.);
    Complex r7 = Complex(0.,0.);
    Complex r8 = Complex(0.,0.);
    // emission off fermions
    for (unsigned int s=0; s<=1; s++) { // spin of pseudo-particle in propagator representation
      r1 = r1 + XYZ.X(5,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.Z(3,m_spins[3],5,s,1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
      r2 = r2 + XYZ.X(5,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.X(3,m_spins[3],q1,5,s,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],q1,2,m_spins[2],m_cL,m_cR);
      r3 = r3 + XYZ.X(6,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.Z(3,m_spins[3],6,s,1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
      r4 = r4 + XYZ.X(6,s,epsP,0,m_spins[0],1.,1.)
                *XYZ.X(3,m_spins[3],q1,6,s,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],q1,2,m_spins[2],m_cL,m_cR);
      r5 = r5 + XYZ.Z(3,m_spins[3],0,m_spins[0],7,s,2,m_spins[2],m_cL,m_cR,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,7,s,1.,1.);
      r6 = r6 + XYZ.X(3,m_spins[3],q2,0,m_spins[0],m_cL,m_cR)
                *XYZ.X(7,s,q2,2,m_spins[2],m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,7,s,1.,1.);
      r7 = r7 + XYZ.Z(3,m_spins[3],0,m_spins[0],8,s,m_spins[2],2,m_cL,m_cR,m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,8,s,1.,1.);
      r8 = r8 + XYZ.X(3,m_spins[3],q2,0,m_spins[0],m_cL,m_cR)
                *XYZ.X(8,s,q2,2,m_spins[2],m_cL,m_cR)
                *XYZ.X(1,m_spins[1],epsP,8,s,1.,1.);
    }
    // add prefactors to those terms
    r1 = -(1+mt/sqrt(qt*qt))*r1;
    r2 =  (1+mt/sqrt(qt*qt))*r2/(mW*mW);
    r3 = -(1-mt/sqrt(qt*qt))*r3;
    r4 =  (1-mt/sqrt(qt*qt))*r4/(mW*mW);
    r5 = -(1+ml/sqrt(ql*ql))*r5;
    r6 =  (1+ml/sqrt(ql*ql))*r6/(mW*mW);
    r7 = -(1-ml/sqrt(ql*ql))*r7;
    r8 =  (1-ml/sqrt(ql*ql))*r8/(mW*mW);
    // calculate terms arising from emission off the W propagator
    Complex Z = XYZ.Z(3,m_spins[3],0,m_spins[0],1,m_spins[1],2,m_spins[2],m_cL,m_cR,m_cL,m_cR);
    Complex X1 = XYZ.X(3,m_spins[3],q2,0,m_spins[0],m_cL,m_cR);
    Complex X2 = XYZ.X(3,m_spins[3],q1,0,m_spins[0],m_cL,m_cR);
    Complex X3 = XYZ.X(3,m_spins[3],epsP,0,m_spins[0],m_cL,m_cR);
    Complex X4 = XYZ.X(3,m_spins[3],q1-k,0,m_spins[0],m_cL,m_cR);
    Complex X5 = XYZ.X(1,m_spins[1],q1,2,m_spins[2],m_cL,m_cR);
    Complex X6 = XYZ.X(1,m_spins[1],q2,2,m_spins[2],m_cL,m_cR);
    Complex X7 = XYZ.X(1,m_spins[1],epsP,2,m_spins[2],m_cL,m_cR);
    Complex X8 = XYZ.X(1,m_spins[1],q2+k,2,m_spins[2],m_cL,m_cR);
    Complex s1 = Complex(0.,0.);
    Complex s2 = Complex(0.,0.);
    Complex s3 = Complex(0.,0.);
    Complex s4 = Complex(0.,0.);
    Complex s5 = Complex(0.,0.);
    Complex s6 = Complex(0.,0.);
    Complex s7 = Complex(0.,0.);
    Complex s8 = Complex(0.,0.);
    Complex s9 = Complex(0.,0.);
    s1 = - Z * ((q1+q2)*epsP);
    s2 = X3*X8;
    s3 = -X4*X7;
    s4 = 1./(mW*mW)*((q1+q2)*epsP)*(X1*X6 + X2*X5);
    s5 = -1./(mW*mW)*(q1*(q2+k))*X3*X5;
    s6 = 1./(mW*mW)*(q2*(q1-k))*X1*X7;
    s7 = -1./(mW*mW)*(q2*epsP)*X1*X8;
    s8 = 1./(mW*mW)*(q1*epsP)*X4*X5;
    s9 = -1./(mW*mW*mW*mW)*((q1*q2)*((q1+q2)*epsP)
                            -(q2*epsP)*(q1*(q2+k))
                            +(q1*epsP)*(q2*(q1-k)))*X1*X5;
    // combined result with prefactors
    Complex r =  (m_i*m_e)/(2*(qt*qt-mt*mt)*(q1*q1-mW*mW))*(r1+r2+r3+r4)
                +(m_i*m_e)/(2*(ql*ql-ml*ml)*(q2*q2-mW*mW))*(r5+r6+r7+r8)
                -(m_i*m_e)/((q1*q1-mW*mW)*(q2*q2-mW*mW))*(s1+s2+s3+s4+s5+s6+s7+s8+s9);
    return r;
  }
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_1_15(unsigned int i) {
  return 0;
}

Complex Tau_To_Lepton_Neutrinos::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0;
}

double Tau_To_Lepton_Neutrinos::GetBeta_0_0() {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin lepton-neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=1; m++) {       // spin tau-neutrino
        for (unsigned int k=0; k<=1; k++) {     // spin tau
          m_spins[0] = k;
          m_spins[1] = j;
          m_spins[2] = i;
          m_spins[3] = m;
          Complex M_0_0 = InfraredSubtractedME_0_0();
          sum = sum + (M_0_0*conj(M_0_0)).real();
        }
      }
    }
  }
  // spin avarage over initial state
  sum = (1./2.)*sum;
  return sum;
}

double Tau_To_Lepton_Neutrinos::GetBeta_0_1() {
  return 0;
}

double Tau_To_Lepton_Neutrinos::GetBeta_0_2() {
  return 0;
}

double Tau_To_Lepton_Neutrinos::GetBeta_1_1(unsigned int a) {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin lepton-neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=1; m++) {       // spin tau-neutrino
        for (unsigned int k=0; k<=1; k++) {     // spin tau
          for (unsigned int l=0; l<=1; l++) {   // spin gamma
            m_spins[0] = k;
            m_spins[1] = j;
            m_spins[2] = i;
            m_spins[3] = m;
            m_spins[4] = l;
            Complex M_1_05 = InfraredSubtractedME_1_05(a);
            sum = sum + (M_1_05*conj(M_1_05)).real();
          }
        }
      }
    }
  }
  // spin avarage over initial state
  sum = (1./2.)*sum;
  sum = 1./(2*pow(2*M_PI,3))*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double Tau_To_Lepton_Neutrinos::GetBeta_1_2(unsigned int i) {
  return 0;
}

double Tau_To_Lepton_Neutrinos::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0;
}
