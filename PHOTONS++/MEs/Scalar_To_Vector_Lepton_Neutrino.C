#include "PHOTONS++/MEs/Scalar_To_Vector_Lepton_Neutrino.H"
#include "ATOOLS/Math/Poincare.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "HELICITIES/Main/Polarization_Tools.H"
#include "HELICITIES/Loops/Master_Integrals.H"
#include "HELICITIES/Loops/PV_Integrals.H"


using namespace PHOTONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

Scalar_To_Vector_Lepton_Neutrino::Scalar_To_Vector_Lepton_Neutrino
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FF(pvv) {
  m_name = "Scalar_To_Vector_Lepton_Neutrino";
  m_flavs[0]  = pvv[1][0]->Flav();    // neutral IS scalar
  m_masses[0] = pvv[1][0]->FinalMass();
  m_flavs[3]  = pvv[3][0]->Flav();    // neutrino (calcs with m_nu = 0)
  m_masses[3] = 0.;
  // switch ordering if necessary (FS charged scalar at pos 1)
  m_switch = pvv[2][0]->Flav().IsLepton();
  // m_switch == true if first charged final state particle is lepton
  if (m_switch == false) {
    m_flavs[1]  = pvv[2][0]->Flav();  // charged FS vector
    m_masses[1] = pvv[2][0]->FinalMass();
    m_flavs[2]  = pvv[2][1]->Flav();  // lepton
    m_masses[2] = pvv[2][1]->FinalMass();
  }
  else {
    m_flavs[1]  = pvv[2][1]->Flav();  // charged FS vector
    m_masses[1] = pvv[2][1]->FinalMass();
    m_flavs[2]  = pvv[2][0]->Flav();  // lepton
    m_masses[2] = pvv[2][0]->FinalMass();
  }
  for (unsigned int i=4; i<9; i++) {
    m_flavs[i]  = Flavour(kf_photon);
    m_masses[i] = 0.;
  }

  m_cL = 1.;
  m_cR = 0.;

  // operate in Form-Factor-Model, point-like otherwise
  m_ffmodel = false;

  m_fpluszero      = 1.;
  m_fplusprimezero = 1.;
}

Scalar_To_Vector_Lepton_Neutrino::~Scalar_To_Vector_Lepton_Neutrino() {
}

void Scalar_To_Vector_Lepton_Neutrino::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into -z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum += m_olddipole[i]->Momentum();
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

void Scalar_To_Vector_Lepton_Neutrino::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  Poincare boost(m_pvv_zero[1][0]->Momentum());
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
    // charged FS vector
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

double Scalar_To_Vector_Lepton_Neutrino::Smod(unsigned int kk) {
  m_moms = m_moms1[kk];
  Vec4D k   = m_moms[4];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = -1.;   // both charged particles of opposite charge
  double Zj = +1.;
  int    ti = +1;   // always both final state
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4C epsV = conj(Polarization_Vector(m_moms[1]).at(m_spins[1]));
  double t = (m_moms[0]-m_moms[1]).Abs2();
  Vec4C auxvec = Contraction(AuxiliaryTensor(m_moms[0],m_moms[1],t),2,epsV);

  XYZFunc XYZ(4,m_moms,m_flavs,1,false);

  if (m_ffmodel == true) {
    return 0;
  }
  else {
    return m_sqrt2*m_GF*XYZ.X(3,m_spins[3],auxvec,2,m_spins[2],m_cL,m_cR);
  }
  return 0;
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_0_1() {
  m_moms = m_moms0;
  return 0;
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_0_2() {
  return 0;
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsV   = conj(Polarization_Vector(m_moms[1]).at(m_spins[1]));
  Vec4C epsP   = conj(Polarization_Vector(m_moms[4]).at(m_spins[4]));
  double mV    = m_masses[1];               // charged FS vector mass/propagator pole
  double ml    = m_masses[2];               // lepton mass/propagator pole
  Vec4D qV     = m_moms[1]+m_moms[4];       // vector propagator momentum
  Vec4D ql     = m_moms[2]+m_moms[4];       // lepton propagator momentum
  m_moms[5]    = m_moms[6] = ql;            // enter those into m_moms
  m_flavs[5]   = m_flavs[2];                // set to corresponding particle/antiparticle
  m_flavs[6]   = m_flavs[2].Bar();
  Vec4D k      = m_moms[4];
  XYZFunc XYZ(7,m_moms,m_flavs,1,false);
  m_flavs[5] = m_flavs[6] = Flavour(kf_none);
  // Form-Factor-Model
  Complex r1 = 0;
  Complex r2 = 0;
  Complex r3 = 0;
  Complex r4 = 0;
  if (m_ffmodel == true) {
    return 0;
  }
  // point-like
  else {
    // A
    double tA = (m_moms[0]-m_moms[1]).Abs2();
    Vec4C auxvecA = Contraction(AuxiliaryTensor(m_moms[0],m_moms[1],tA),2,epsV);
    for (unsigned int s=0; s<=1; s++) {
      r1 = r1 + XYZ.X(3,m_spins[3],auxvecA,5,s,m_cL,m_cR)
                * XYZ.X(5,s,epsP,2,m_spins[2],1.,1.);
      r2 = r2 + XYZ.X(3,m_spins[3],auxvecA,6,s,m_cL,m_cR)
                * XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    }
    r1 = -(0.5*m_sqrt2*m_e*m_GF)/(ql.Abs2()-ml*ml) * (1-ml/ql.Abs()) * r1;
    r2 = -(0.5*m_sqrt2*m_e*m_GF)/(ql.Abs2()-ml*ml) * (1+ml/ql.Abs()) * r2;
    // B
    double tB = (m_moms[0]-qV).Abs2();
    Vec4C auxvecB =  Contraction(AuxiliaryTensor(m_moms[0],qV,tB),2,epsP)*((m_moms[1]+2*k)*epsV)
                    +Contraction(AuxiliaryTensor(m_moms[0],qV,tB),2,(m_moms[1]-k))*(epsV*epsP)
                    +Contraction(AuxiliaryTensor(m_moms[0],qV,tB),2,epsV)*((-2*m_moms[1]-k)*epsP)
                    -1./(mV*mV)*Contraction(AuxiliaryTensor(m_moms[0],qV,tB),2,qV)
                      *(  (qV*epsP)*((m_moms[1]+2*k)*epsV)
                         +(epsV*epsP)*(qV*(m_moms[1]-k))
                         +(qV*epsV)*((-2*m_moms[1]-k)*epsP) );
    r3 = (m_sqrt2*m_e*m_GF)/(qV.Abs2()-mV*mV)
          * XYZ.X(3,m_spins[3],auxvecB,2,m_spins[2],m_cL,m_cR);
    return r1+r2+r3;
  }
  return 0;
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_1_15(unsigned int i) {
  return 0;
}

Complex Scalar_To_Vector_Lepton_Neutrino::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0;
}

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_0_0() {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=2; m++) {       // spin FS vector
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

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_0_1() {
  return 0;
//   double sum = 0;
//   for (unsigned int i=0; i<=1; i++) {           // spin neutrino
//     for (unsigned int j=0; j<=1; j++) {         // spin lepton
//       for (unsigned int m=0; m<=2; m++) {       // spin FS vector
//         for (unsigned int k=0; k<=0; k++) {     // spin IS scalar
//           m_spins[0] = k;
//           m_spins[1] = m;
//           m_spins[2] = j;
//           m_spins[3] = i;
//           Complex M_0_0 = InfraredSubtractedME_0_0();
//           Complex M_0_1 = InfraredSubtractedME_0_1();
//           sum = sum + (M_0_0*conj(M_0_1)+M_0_1*conj(M_0_0)).real();
//         }
//       }
//     }
//   }
//   // spin avaraging gives factor 1
//   return sum;
}

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_0_2() {
  return 0;
}

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_1_1(unsigned int a) {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin neutrino
    for (unsigned int j=0; j<=1; j++) {         // spin lepton
      for (unsigned int m=0; m<=2; m++) {       // spin FS vector
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

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_1_2(unsigned int i) {
  return 0;
}

double Scalar_To_Vector_Lepton_Neutrino::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0;
}

double Scalar_To_Vector_Lepton_Neutrino::FormFactorAp(double t) {
  return -0.0851648 - 0.00291091*t - 0.0000786814*t*t;
}

double Scalar_To_Vector_Lepton_Neutrino::FormFactorAm(double t) {
  return 0.;
}

double Scalar_To_Vector_Lepton_Neutrino::FormFactorF(double t) {
  return 4.8048 + 0.0914977*t + 0.00174239*t*t;
}

double Scalar_To_Vector_Lepton_Neutrino::FormFactorG(double t) {
  return 0.104396 + 0.00400149*t + 0.000117909*t*t;
}

Lorentz_Ten2D Scalar_To_Vector_Lepton_Neutrino::AuxiliaryTensor(Vec4D p1, Vec4D p2, double t) {
  // convert t to dimensionless quantity tt=t/[MeV]
  double tt = t*(1000.*1000.);
  return  FormFactorG(tt)*EpsilonTensorContraction(p1,p2)
         -FormFactorF(tt)*MetricTensor()
         -FormFactorAp(tt)*BuildTensor(p1+p2,p1-p2)
         -FormFactorAm(tt)*BuildTensor(p1-p2,p1+p2);
}


DECLARE_PHOTONS_ME_GETTER(Scalar_To_Vector_Lepton_Neutrino_Getter,
                          "Scalar_To_Vector_Lepton_Neutrino")

PHOTONS_ME_Base * Scalar_To_Vector_Lepton_Neutrino_Getter::operator()
(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 0) &&
       (pvv[1].size() == 1) && pvv[1][0]->Flav().IsScalar()        &&
      (((pvv[2].size() == 2) && pvv[2][0]->Flav().IsVector() &&
                                pvv[2][1]->Flav().IsLepton())  ||
       ((pvv[2].size() == 2) && pvv[2][1]->Flav().IsVector() &&
                                pvv[2][0]->Flav().IsLepton())    ) &&
       (pvv[3].size() == 1) && pvv[3][0]->Flav().IsLepton())
    return new Scalar_To_Vector_Lepton_Neutrino(pvv);
  return NULL;
}
