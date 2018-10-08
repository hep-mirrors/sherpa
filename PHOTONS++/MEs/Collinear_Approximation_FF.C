#include "PHOTONS++/MEs/Collinear_Approximation_FF.H"
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

// Collinear Approximation for FF dipoles
Collinear_Approximation_FF::Collinear_Approximation_FF
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FF(pvv) {
  Data_Reader reader(" ",";","#","=");
  m_name = "Collinear_Approximation_FF";
  m_no_weight = 0;
  m_flavs[0]  = pvv[1][0]->Flav();
  m_masses[0] = pvv[1][0]->FinalMass();
  // switch ordering if necessary
  m_switch = pvv[2][0]->Flav().IsAnti();
  if (m_switch == false) {
    m_flavs[1] = pvv[2][0]->Flav(); m_masses[1] = pvv[2][0]->FinalMass();
    m_flavs[2] = pvv[2][1]->Flav(); m_masses[2] = pvv[2][1]->FinalMass();
  }
  else {
    m_flavs[2] = pvv[2][0]->Flav(); m_masses[2] = pvv[2][0]->FinalMass();
    m_flavs[1] = pvv[2][1]->Flav(); m_masses[1] = pvv[2][1]->FinalMass();
  }
}

Collinear_Approximation_FF::~Collinear_Approximation_FF() {
}

void Collinear_Approximation_FF::BoostOriginalPVVToMultipoleCMS() {
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
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,1.));
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

void Collinear_Approximation_FF::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  m_moms0[0] = m_pvv_zero[1][0]->Momentum();
  if (m_switch == false) {
    m_moms0[1] = m_pvv_zero[2][0]->Momentum();
    m_moms0[2] = m_pvv_zero[2][1]->Momentum();
  }
  else {
    m_moms0[2] = m_pvv_zero[2][0]->Momentum();
    m_moms0[1] = m_pvv_zero[2][1]->Momentum();
  }

  // Store fully rescaled final state
  m_momsFull[0] = pvv_one[1][0]->Momentum();
  if (m_switch == false) {
    m_momsFull[1] = pvv_one[2][0]->Momentum();
    m_momsFull[2] = pvv_one[2][1]->Momentum();
  }
  else {
    m_momsFull[2] = pvv_one[2][0]->Momentum();
    m_momsFull[1] = pvv_one[2][1]->Momentum();
  }
  for (size_t i = 0; i < pvv_one[4].size(); i++) {
    m_momsFull[i+3] = pvv_one[4][i]->Momentum();
  }

  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated
  // Need to define dipole first and only once!
  Dipole_FF::DefineDipole();

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[1][0]->Momentum();
    if (m_switch == false) {
      m_moms1[0][1] = pvv_one[2][0]->Momentum();
      m_moms1[0][2] = pvv_one[2][1]->Momentum();
    }
    else {
      m_moms1[0][2] = pvv_one[2][0]->Momentum();
      m_moms1[0][1] = pvv_one[2][1]->Momentum();
    }
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
      if (m_switch == false) {
        m_moms1[i][1] = m_newdipole[0]->Momentum();
        m_moms1[i][2] = m_newdipole[1]->Momentum();
      }
      else {
        m_moms1[i][2] = m_newdipole[0]->Momentum();
        m_moms1[i][1] = m_newdipole[1]->Momentum();
      }
      m_moms1[i][3] = m_softphotons[0]->Momentum();
      m_moms1[i][0] = m_moms1[i][1]+m_moms1[i][2]+m_moms1[i][3];
      m_softphotons.clear();
    }
  }
}

double Collinear_Approximation_FF::SmodFull(unsigned int kk) {
  // Smod calculated using fully dressed momenta 
  // Used to cancel the Smod in the dipole weight

  // Index kk denotes which photon is used in calculation. Note that m_momsFull contains all photons
  // starting at position 3, hence index of photon kk is 3+kk.
  m_moms = m_momsFull;
  Vec4D k   = m_moms[3+kk];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  // Note the way this is set up, SmodFull does not need m_alpha/4pi^2 since there's none in
  // GetBeta_1_1 either
  return Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double Collinear_Approximation_FF::Smod(unsigned int kk) {
  return 0.;
}

double Collinear_Approximation_FF::Smod(unsigned int a, unsigned int b, unsigned int c) {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_0_0() {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_0_1(const int& ewmode) {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_1_05(unsigned int i, const double& xi) {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_1_15(unsigned int a, const double& xi) {
  return 0.;
}

Complex Collinear_Approximation_FF::InfraredSubtractedME_2_1(unsigned int a, unsigned int b) {
  return 0.;
}

double Collinear_Approximation_FF::GetBeta_0_0() {
  return 0.;
}

double Collinear_Approximation_FF::GetBeta_0_1(const int& ewmode) {
  return 0.;
}

double Collinear_Approximation_FF::GetBeta_0_2() {
  return 0.;
}
  
double Collinear_Approximation_FF::GetBeta_1_1(unsigned int a) {
  return Dmod(0,1,a)+Dmod(1,0,a);
}

double Collinear_Approximation_FF::Dmod(unsigned int i, unsigned int j, unsigned int kk) {
  m_moms = m_moms1[kk];
  Vec4D pi = m_moms[1+i];
  Vec4D pj = m_moms[1+j];
  Vec4D k  = m_moms[3];
  double D = 1./(pi*k)*(2.*(pi*pj)/(pi*k+pj*k)-(pi*pi)/(pi*k));
  // emitter and spectator final state
  double y = (pi*k)/(pi*k+pj*k+pi*pj);
  double z = (pi*pj)/(pi*pj+pj*k);
  double p = (pi+pj+k)*(pi+pj+k);
  double P = p-pi*pi-pj*pj-k*k;
  double s = sqrt(sqr(2.*(pj*pj)+P-P*y)-4.*p*(pj*pj));
  double R = s/sqrt(Kallen(p,pi*pi,pj*pj));
  if (m_newdipole[i]->Flav().IntSpin() == 0) {
    return 0.;
  }
  else if (m_newdipole[i]->Flav().IntSpin() == 1) {
    return 1./((pi*k)*R)*(2./(1.-z*(1.-y))-1.-z-(pi*pi)/(pi*k)) - D;
  }
  else if (m_newdipole[i]->Flav().IntSpin() == 2) {
    return 1./(pi*k) * (2./(1.-z*(1.-y))+2./(1.-(1.-z)*(1.-y))+2.*z*(1.-z)
			-4.-(pi*pi)/(pi*k)) - D;
  }
  else if (m_newdipole[i]->Flav().IntSpin() == 3) {
    return 0.;
  }
}
 
 double Collinear_Approximation_FF::GetBeta_1_2(unsigned int a) {
   return 0.;
}

double Collinear_Approximation_FF::GetBeta_2_2(unsigned int a, unsigned int b) {
  return 0.;
}
 
void Collinear_Approximation_FF::Print_Info() {
  return;
}

double Collinear_Approximation_FF::Kallen(double x, double y, double z) {
  return ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z );
}

DECLARE_PHOTONS_ME_GETTER(Collinear_Approximation_FF,
                          "Collinear_Approximation_FF")

PHOTONS_ME_Base *ATOOLS::Getter<PHOTONS_ME_Base,Particle_Vector_Vector,
				Collinear_Approximation_FF>::
operator()(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 0) &&
       (pvv[1].size() == 1) && (pvv[1][0]->Flav().Mass() != 0.) &&
       (pvv[2].size() == 2) &&
       (pvv[3].size() == 0))
    return new Collinear_Approximation_FF(pvv);
  return NULL;
}
