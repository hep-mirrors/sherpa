#include "PHOTONS++/MEs/Scalar_To_Scalar_Lepton_Neutrino.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define B_1(A,B,C,M)        PV_Bubble_1(A,B,C,M)
#define C_11(A,B,C,D,E,F,M) PV_Triangle_11(A,B,C,D,E,F,M)
#define C_12(A,B,C,D,E,F,M) PV_Triangle_12(A,B,C,D,E,F,M)
#define C_21(A,B,C,D,E,F,M) PV_Triangle_21(A,B,C,D,E,F,M)
#define C_22(A,B,C,D,E,F,M) PV_Triangle_22(A,B,C,D,E,F,M)
#define C_23(A,B,C,D,E,F,M) PV_Triangle_23(A,B,C,D,E,F,M)
#define C_24(A,B,C,D,E,F,M) PV_Triangle_24(A,B,C,D,E,F,M)

using namespace PHOTONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Scalar_To_Scalar_Lepton_Neutrino::Scalar_To_Scalar_Lepton_Neutrino
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FF(pvv) {
  m_name = "Scalar_To_Scalar_Lepton_Neutrino";
  m_flavs[0]  = pvv[1][0]->Flav();    // IS scalar (neutral)
  m_masses[0] = pvv[1][0]->FinalMass();
  m_flavs[3]  = pvv[3][0]->Flav();    // neutrino (calcs with m_nu = 0)
  m_masses[3] = 0.;
  // switch ordering if necessary (FS charged scalar at pos 1)
  m_switch = pvv[2][0]->Flav().IsLepton();
  // m_switch == true if first charged final state particle is lepton
  if (m_switch == false) {
    m_flavs[1]  = pvv[2][0]->Flav();  // charged FS scalar
    m_masses[1] = pvv[2][0]->FinalMass();
    m_flavs[2]  = pvv[2][1]->Flav();  // lepton
    m_masses[2] = pvv[2][1]->FinalMass();
  }
  else {
    m_flavs[1]  = pvv[2][1]->Flav();  // charged FS scalar
    m_masses[1] = pvv[2][1]->FinalMass();
    m_flavs[2]  = pvv[2][0]->Flav();  // lepton
    m_masses[2] = pvv[2][0]->FinalMass();
  }
  for (unsigned int i=4; i<9; i++) {
    m_flavs[i] = Flavour(kf_photon);
    m_masses[i] = 0.;
  }

  m_cL = 1.;
  m_cR = 0.;

  // operate in Form-Factor-Model, point-like otherwise
  m_ffmodel = false;

  m_fpluszero      = 1.;
  m_fplusprimezero = 1.;

  for (unsigned int i=0; i<=1; i++)            // spin neutrino
    for (unsigned int j=0; j<=1; j++)          // spin lepton
      for (unsigned int m=0; m<=0; m++)        // spin FS scalar
        for (unsigned int k=0; k<=0; k++)      // spin IS scalar
          m_M00results[k][m][j][i].first = false;
}

Scalar_To_Scalar_Lepton_Neutrino::~Scalar_To_Scalar_Lepton_Neutrino() {
}

void Scalar_To_Scalar_Lepton_Neutrino::BoostOriginalPVVToMultipoleCMS() {
  // pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into -z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum = sum + m_olddipole[i]->Momentum();
  }
  Vec4D p1 = m_olddipole[0]->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,-1.));
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D vec = m_olddipole[i]->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_olddipole[i]->SetMomentum(vec);
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D vec = m_oldspectator[i]->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_oldspectator[i]->SetMomentum(vec);
  }
}

void Scalar_To_Scalar_Lepton_Neutrino::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  Poincare boost(m_pvv_zero[1][0]->Momentum());
  // m_moms0 - no photon
  Vec4D vec(0.,0.,0.,0.);
  // neutral IS scalar
  vec = m_pvv_zero[1][0]->Momentum();
  boost.Boost(vec);
  m_moms0[0] = vec;
  // neutrino
  vec = m_pvv_zero[3][0]->Momentum();
  boost.Boost(vec);
  m_moms0[3] = vec;
  if (m_switch == false) {
    // charged FS scalar
    vec = m_pvv_zero[2][0]->Momentum();
    boost.Boost(vec);
    m_moms0[1] = vec;
    // lepton
    vec = m_pvv_zero[2][1]->Momentum();
    boost.Boost(vec);
    m_moms0[2] = vec;
  }
  else {
    // charged FS scalar
    vec = m_pvv_zero[2][1]->Momentum();
    boost.Boost(vec);
    m_moms0[1] = vec;
    // lepton
    vec = m_pvv_zero[2][0]->Momentum();
    boost.Boost(vec);
    m_moms0[2] = vec;
  }

  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[1][0]->Momentum();
    m_moms1[0][3] = pvv_one[3][0]->Momentum();
    if (m_switch == false) {
      m_moms1[0][1] = pvv_one[2][0]->Momentum();
      m_moms1[0][2] = pvv_one[2][1]->Momentum();
    }
    else {
      m_moms1[0][1] = pvv_one[2][1]->Momentum();
      m_moms1[0][2] = pvv_one[2][0]->Momentum();
    }
    m_moms1[0][4] = pvv_one[4][0]->Momentum();
  }
  else {
    Dipole_FF::DefineDipole();
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<pvv_one[4].size(); i++) {
      m_softphotons.push_back(pvv_one[4][i]);
      m_K = CalculateMomentumSum(m_softphotons);
      CorrectMomenta();
      if (m_switch == false) {
        m_moms1[i][1] = m_newdipole[0]->Momentum();
        m_moms1[i][2] = m_newdipole[1]->Momentum();
      }
      else {
        m_moms1[i][1] = m_newdipole[1]->Momentum();
        m_moms1[i][2] = m_newdipole[0]->Momentum();
      }
      m_moms1[i][3] = m_newspectator[0]->Momentum();
      m_moms1[i][4] = m_softphotons[0]->Momentum();
      m_moms1[i][0] = m_moms1[i][1]+m_moms1[i][2]+m_moms1[i][3]+m_moms1[i][4];
      m_softphotons.clear();
    }
  }
}

double Scalar_To_Scalar_Lepton_Neutrino::Smod(unsigned int kk) {
  m_moms = m_moms1[kk];
  Vec4D k   = m_moms[4];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = -1.;   // both charged particles of opposite charge
  double Zj = +1.;
  int    ti = +1;    // always both final state
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_0_0() {
  // if result has already been calculated, return this
  if (m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].first)
    return m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].second;
  m_moms = m_moms0;
  double t = (m_moms[0]-m_moms[1]).Abs2();
  XYZFunc XYZ(4,m_moms,m_flavs,false);
  m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].first = true;
  if (m_ffmodel == true) {
    m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].second
      = m_sqrt2*m_GF
          *(Fplus(t)*XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1],2,m_spins[2],m_cR,m_cL)
           -Fminus(t)*XYZ.X(3,m_spins[3],m_moms[0]-m_moms[1],2,m_spins[2],m_cR,m_cL));
  }
  else {
    m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].second
      = m_sqrt2*m_GF*XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1],2,m_spins[2],m_cR,m_cL);
  }
  return m_M00results[m_spins[0]][m_spins[1]][m_spins[2]][m_spins[3]].second;
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_0_1() {
  if (m_ffmodel == true) return 0.;
  m_moms = m_moms0;
  // add pseudo-fermion to use XYZ-functions
  Flavour flavs[5]; Vec4D moms[5];
  for (unsigned int i(0);i<4;++i) { flavs[i]=m_flavs[i]; moms[i]=m_moms[i]; }
  flavs[4] = Flavour(kf_e).Bar(); moms[4]=m_moms[0];
  XYZFunc XYZ(5,moms,flavs,false);
  double t((m_moms[0]-m_moms[1]).Abs2());
  double mu2(t);
  Vec4D p(m_moms[0]),p1(m_moms[1]),p2(m_moms[2]),p3(m_moms[3]);
  double m12(sqr(m_masses[1])),m22(sqr(m_masses[2]));
  double m2(sqrt(m22));
  double p1p2(p1*p2), pp((p1+p2).Abs());
  DivArrC term1(0.,0.,0.,0.,0.,0.),term2(0.,0.,0.,0.,0.,0.),
          term3(0.,0.,0.,0.,0.,0.),term4(0.,0.,0.,0.,0.,0.);
  // ~ u_3 (p+p1)P_Lv_2
  term1 = 0.25*B_0(pp,m12,m22,mu2)
          +0.5*(m12-p1p2)*C_11(m12,m22,pp,0.,m12,m22,mu2)
          +(p1p2-0.25*m22)*C_12(m12,m22,pp,0.,m12,m22,mu2)
          +0.25*B_0(0.,m12,m12,mu2)
          +(D-2.)/8.*B_0(0.,m22,m22,mu2)
          +(D-2.)/10.*A_0(m22,mu2)/m22
          -0.25*(B_0(m22,0.,m22,mu2)-B_0(0.,m22,m22,mu2))
          +0.125*(2.*B_0(pp,m12,m22,mu2)
                  -B_0(0.,m12,m12,mu2)-B_0(0.,m22,m22,mu2));
  term1 *= XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1],2,m_spins[2],m_cR,m_cL);
  // ~ u_3 p p_1 P_R v_2
  term2 = 0.5*m2*C_12(m12,m22,pp,0,m12,m22,mu2);
  Complex temp(0.,0.);
  for (unsigned short int s(0);s<=1;++s)
    temp += XYZ.Y(3,m_spins[3],4,s,1.,1.)
            *XYZ.X(4,s,m_moms[1],2,m_spins[2],m_cL,m_cR);
  term2 *= temp;
  // ~ u_3 p_1 P_L v_2
  term3 = 0.25*B_0(pp,m12,m22,mu2)
          +0.25*B_1(pp,m12,m22,mu2)
          -p1p2*C_11(m12,m22,pp,0.,m12,m22,mu2)
          -0.5*p1p2*C_21(m12,m22,pp,0.,m12,m22,mu2)
          -0.5*m22*C_23(m12,m22,pp,0.,m12,m22,mu2);
  term3 *= XYZ.X(3,m_spins[3],m_moms[1],2,m_spins[2],m_cR,m_cL);
  // ~ u_3 P_R v_2
  term4 = -0.25*m2*B_1(pp,m12,m22,mu2)
          +m2*(0.5*m12*+p1p2)*C_12(m12,m22,pp,0.,m12,m22,mu2)
          +0.5*m2*m22*C_22(m12,m22,pp,0.,m12,m22,mu2)
          +0.5*m2*p1p2*C_23(m12,m22,pp,0.,m12,m22,mu2)
          +0.5*m2*C_24(m12,m22,pp,0.,m12,m22,mu2);
  term4 *= XYZ.Y(3,m_spins[3],2,m_spins[2],m_cL,m_cR);
  return m_sqrt2*m_GF*m_alpha/M_PI*(term1+term2+term3+term4).Finite();
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsP   = conj(Polarization_Vector(m_moms[4]).at(m_spins[4]));
  double mS    = m_masses[1];      // charged FS scalar mass/propagator pole
  double ml    = m_masses[2];      // lepton mass/propagator pole
  Vec4D ql     = m_moms[2]+m_moms[4];       // lepton propagator momenta
  m_moms[5]    = m_moms[6] = ql;            // enter those into m_moms
  m_flavs[5]   = m_flavs[2];                // set to corresponding particle/antiparticle
  m_flavs[6]   = m_flavs[2].Bar();
  Vec4D k      = m_moms[4];
  XYZFunc XYZ(7,m_moms,m_flavs,false);
  m_flavs[5] = m_flavs[6] = Flavour(kf_none);
  // Form-Factor-Model
  Complex r1(0.,0.), r2(0.,0.), r3(0.,0.), r4(0.,0.);
  if (m_ffmodel == true) {
    double t1 = (m_moms[2]+m_moms[3]).Abs2();
    double t2 = (m_moms[0]-m_moms[1]).Abs2();
    r1 = XYZ.X(3,m_spins[3],(2.*m_moms[0]-m_moms[2])*Fplus(t1)
                              +m_moms[2]*Fminus(t1),2,m_spins[2],m_cR,m_cL);
    r1 = m_sqrt2*m_e*m_GF*(m_moms[1]*epsP)/(m_moms[1]*m_moms[4])*r1;
    for (unsigned int s=0; s<=1; s++) {
      r2 = r2 + XYZ.X(3,m_spins[3],(2.*m_moms[0]-ql)*Fplus(t2)
                                    +ql*Fminus(t2),5,s,m_cR,m_cL)
                *XYZ.X(5,s,epsP,2,m_spins[2],1,1);
      r3 = r3 + XYZ.X(3,m_spins[3],(2.*m_moms[0]-ql)*Fplus(t2)
                                    +ql*Fminus(t2),6,s,m_cR,m_cL)
                *XYZ.X(6,s,epsP,2,m_spins[2],1,1);
    }
    r2 = -m_sqrt2*m_e*m_GF*(1.-ml/sqrt(ql.Abs2()))/(4.*m_moms[2]*m_moms[4])*r2;
    r3 = -m_sqrt2*m_e*m_GF*(1.+ml/sqrt(ql.Abs2()))/(4.*m_moms[2]*m_moms[4])*r3;
    r4 = XYZ.X(3,m_spins[3],epsP*(Fplus(t2)-Fminus(t2))
                            -2.*(2.*m_moms[0]-m_moms[2])*m_fplusprimezero
                               *((m_moms[0]-m_moms[1])*epsP)
                                                      ,2,m_spins[2],m_cR,m_cL);
    r4 = -m_sqrt2*m_e*m_GF*r4;
    return r1+r2+r3+r4;
  }
  // point-like
  else {
    r1 = XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1]+m_moms[4]
                                                      ,2,m_spins[2],m_cR,m_cL);
    r1 = ((2.*m_moms[1]+m_moms[4])*epsP)*r1;
    r1 = m_sqrt2*m_e*m_GF/((m_moms[1]+m_moms[4]).Abs2()-mS*mS)*r1;
    for (unsigned int s=0; s<=1; s++) {
      r2 = r2 + XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1],5,s,m_cR,m_cL)
                *XYZ.X(5,s,epsP,2,m_spins[2],1.,1.);
      r3 = r3 + XYZ.X(3,m_spins[3],m_moms[0]+m_moms[1],6,s,m_cR,m_cL)
                *XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    }
    r2 = (1.-ml/sqrt(ql.Abs2()))*r2;
    r3 = (1.+ml/sqrt(ql.Abs2()))*r3;
    r2 = -m_sqrt2*m_e*m_GF/(2.*((m_moms[2]+m_moms[4]).Abs2()-ml*ml))*r2;
    r3 = -m_sqrt2*m_e*m_GF/(2.*((m_moms[2]+m_moms[4]).Abs2()-ml*ml))*r3;
    return r1+r2+r3;
  }
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_1_15(unsigned int i) {
  return 0.;
}

Complex Scalar_To_Scalar_Lepton_Neutrino::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0.;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_0_0() {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=0; m++) {       // spin FS scalar
        for (unsigned int k=0; k<=0; k++) {     // spin IS scalar
          m_spins[0] = k;
          m_spins[1] = m;
          m_spins[2] = j;
          m_spins[3] = i;
          Complex M_0_0 = InfraredSubtractedME_0_0();
          sum = sum + (M_0_0*conj(M_0_0)).real();
        }
      }
    }
  }
  // spin avaraging gives factor 1
  return sum;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_0_1() {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=0; m++) {       // spin FS scalar
        for (unsigned int k=0; k<=0; k++) {     // spin IS scalar
          m_spins[0] = k;
          m_spins[1] = m;
          m_spins[2] = j;
          m_spins[3] = i;
          Complex M_0_0 = InfraredSubtractedME_0_0();
          Complex M_0_1 = InfraredSubtractedME_0_1();
          sum = sum + 2.*(M_0_0*conj(M_0_1)).real();
        }
      }
    }
  }
  // spin avaraging gives factor 1
  return sum;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_0_2() {
  return 0.;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_1_1(unsigned int a) {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=0; m++) {       // spin FS scalar
        for (unsigned int k=0; k<=0; k++) {     // spin IS scalar
          for (unsigned int l=0; l<=1; l++) {   // spin gamma
            m_spins[0] = k;
            m_spins[1] = m;
            m_spins[2] = j;
            m_spins[3] = i;
            m_spins[4] = l;
            Complex M_1_05 = InfraredSubtractedME_1_05(a);
            sum = sum + (M_1_05*conj(M_1_05)).real();
          }
        }
      }
    }
  }
  // spin avaraging gives factor 1
  sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_1_2(unsigned int i) {
  return 0.;
}

double Scalar_To_Scalar_Lepton_Neutrino::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0.;
}

double Scalar_To_Scalar_Lepton_Neutrino::Fplus(double t) {
  // linear form factor model
  return m_fpluszero*(1.+m_fplusprimezero*t);
}

double Scalar_To_Scalar_Lepton_Neutrino::Fminus(double t) {
  return 0.;
}

DECLARE_PHOTONS_ME_GETTER(Scalar_To_Scalar_Lepton_Neutrino_Getter,
                          "Scalar_To_Scalar_Lepton_Neutrino")

PHOTONS_ME_Base * Scalar_To_Scalar_Lepton_Neutrino_Getter::operator()
(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 0) &&
       (pvv[1].size() == 1) && pvv[1][0]->Flav().IsScalar()        &&
      (((pvv[2].size() == 2) && pvv[2][0]->Flav().IsScalar() &&
                                pvv[2][1]->Flav().IsLepton())  ||
       ((pvv[2].size() == 2) && pvv[2][1]->Flav().IsScalar() &&
                                pvv[2][0]->Flav().IsLepton())    ) &&
       (pvv[3].size() == 1) && pvv[3][0]->Flav().IsLepton())
    return new Scalar_To_Scalar_Lepton_Neutrino(pvv);
  return NULL;
}
