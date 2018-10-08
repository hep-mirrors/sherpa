#include "PHOTONS++/MEs/Collinear_Approximation_FI.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

// Collinear approximation for final initial dipoles
Collinear_Approximation_FI::Collinear_Approximation_FI
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FI(pvv) {
  Data_Reader reader(" ",";","#","=");
  m_name = "Collinear_Approximation_FI";
  m_no_weight = 0;
  m_flavs[0]  = pvv[0][0]->Flav();
  m_masses[0] = pvv[0][0]->FinalMass();
  m_flavs[1]  = pvv[2][0]->Flav();
  m_masses[1] = pvv[2][0]->FinalMass();
  m_flavs[2]  = pvv[3][0]->Flav();
  m_masses[2] = pvv[3][0]->FinalMass();
}

Collinear_Approximation_FI::~Collinear_Approximation_FI() {
}

void Collinear_Approximation_FI::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into +z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum += m_olddipole[i]->Momentum();
  }
  Vec4D p1 = m_olddipole[0]->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,-1.)); // in Dipole_FI W is rotated into -z-direction!!!
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

void Collinear_Approximation_FI::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  Poincare boost(m_pvv_zero[0][0]->Momentum());
  Vec4D vec;
  vec = m_pvv_zero[0][0]->Momentum();
  boost.Boost(vec);
  m_moms0[0] = vec;
  vec = m_pvv_zero[2][0]->Momentum();
  boost.Boost(vec);
  m_moms0[1] = vec;
  vec = m_pvv_zero[3][0]->Momentum();
  boost.Boost(vec);
  m_moms0[2] = vec;
  vec = m_pvv_zero[0][0]->Momentum();
  // Store fully rescaled final state
  m_momsFull[0] = pvv_one[2][0]->Momentum();
  m_momsFull[1] = pvv_one[2][1]->Momentum();
  m_momsFull[2] = pvv_one[3][0]->Momentum();
  for (size_t i = 0; i < pvv_one[4].size(); i++) {
    m_momsFull[i+3] = pvv_one[4][i]->Momentum();
  }


  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated
  Dipole_FI::DefineDipole();

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[2][0]->Momentum();
    m_moms1[0][1] = pvv_one[2][1]->Momentum();
    m_moms1[0][2] = pvv_one[3][0]->Momentum();
    m_moms1[0][3] = pvv_one[4][0]->Momentum();
  }
  else if (pvv_one[4].size() > 1) {
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<pvv_one[4].size(); i++) {
      m_softphotons.push_back(pvv_one[4][i]);
      m_K = CalculateMomentumSum(m_softphotons);
      DetermineQAndKappa();
      CorrectMomenta();
      if (m_u == -1.) {
	// if any u = -1. is found, momenta are not corrected properly since root could not be found
	// return no weights in this case!
	m_no_weight = 1;
	continue;
      }
      m_moms1[i][0] = m_newdipole[0]->Momentum();
      m_moms1[i][1] = m_newdipole[1]->Momentum();
      m_moms1[i][2] = m_newspectator[0]->Momentum();
      m_moms1[i][3] = m_softphotons[0]->Momentum();
      m_softphotons.clear();
    }
  }
}

double Collinear_Approximation_FI::SmodFull(unsigned int kk) {
  // Smod calculated using fully dressed momenta 
  // Used to cancel the Smod in the dipole weight

  // Index kk denotes which photon is used in calculation. Note that m_momsFull contains all photons
  // starting at position 3, hence index of photon kk is 3+kk.
  m_moms = m_momsFull;
  Vec4D k   = m_moms[3+kk];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = -1.;
  double Zj = -1.;
  int    ti = -1;
  int    tj = +1;
  // Note the way this is set up, SmodFull does not need m_alpha/4pi^2 since there's none in
  // GetBeta_1_1 either
  return Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double Collinear_Approximation_FI::Smod(unsigned int kk) {
  return 0.;
}

double Collinear_Approximation_FI::Smod(unsigned int a, unsigned int b, unsigned int c) {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_0_0() {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_0_1(const int& ewmode) {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_1_05(unsigned int i, const double& xi) {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_1_15(unsigned int a, const double& xi) {
  return 0.;
}

Complex Collinear_Approximation_FI::InfraredSubtractedME_2_1(unsigned int a, unsigned int b) {
  return 0.;
}

double Collinear_Approximation_FI::GetBeta_0_0() {
  return 0.;
}

double Collinear_Approximation_FI::GetBeta_0_1(const int& ewmode) {
  return 0.;
}

double Collinear_Approximation_FI::GetBeta_0_2() {
  return 0.;
}
  
double Collinear_Approximation_FI::GetBeta_1_1(unsigned int a) {
  return Dmod(0,1,a)+Dmod(1,0,a);
}

double Collinear_Approximation_FI::Dmod(unsigned int i, unsigned int j, unsigned int kk) {
  m_moms = m_moms1[kk];
  Vec4D pi = m_moms[i];
  Vec4D pj = m_moms[j];
  Vec4D k  = m_moms[3];
  double D = 1./(pi*k)*(2.*(pi*pj)/(pi*k+pj*k)-(pi*pi)/(pi*k));
  if ((m_newdipole[i]->DecayBlob() == m_newdipole[j]->ProductionBlob()) &&
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
      // Approximation not justified for large mass vectors, e. g. W-bosons
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
}
 
 double Collinear_Approximation_FI::GetBeta_1_2(unsigned int a) {
   return 0.;
}

double Collinear_Approximation_FI::GetBeta_2_2(unsigned int a, unsigned int b) {
  return 0.;
}
 
void Collinear_Approximation_FI::Print_Info() {
  return;
}

double Collinear_Approximation_FI::Kallen(double x, double y, double z) {
  return ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z );
}

DECLARE_PHOTONS_ME_GETTER(Collinear_Approximation_FI,
                          "Collinear_Approximation_FI")

PHOTONS_ME_Base *ATOOLS::Getter<PHOTONS_ME_Base,Particle_Vector_Vector,
				Collinear_Approximation_FI>::
operator()(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv[0].size() == 1) && (pvv[0][0]->Flav().Mass() != 0.) &&
       (pvv[1].size() == 0) &&
       (pvv[2].size() == 1) && 
       (pvv[3].size() == 1))
    return new Collinear_Approximation_FI(pvv);
  return NULL;
}
