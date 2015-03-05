#include "PHOTONS++/Tools/Dipole_FF.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/PhaseSpace/Avarage_Photon_Number.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Dipole_FF::Dipole_FF(const Particle_Vector_Vector& pvv) {
  m_pvv                 = pvv;
  m_dtype               = Dipole_Type::ff;
  m_chargedinparticles  = pvv[0];
  m_neutralinparticles  = pvv[1];
  m_chargedoutparticles = pvv[2];
  m_neutraloutparticles = pvv[3];
  m_M                   = m_neutralinparticles[0]->Momentum().Mass();
  m_Q                   = Vec4D(0.,0.,0.,0.);
  m_QN                  = Vec4D(0.,0.,0.,0.);
  m_kappaC              = Vec3D(0.,0.,0.);
  m_kappaN              = Vec3D(0.,0.,0.);
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++) {
    m_mC.push_back(m_chargedoutparticles[i]->FinalMass());
  }
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++) {
    m_mN.push_back(m_neutraloutparticles[i]->FinalMass());
  }
  double sum = 0.;
  for (unsigned int i=0; i<m_mC.size(); i++) {
    sum = sum + m_mC[i];
  }
  for (unsigned int i=0; i<m_mN.size(); i++) {
    sum = sum + m_mN[i];
  }
  m_omegaMax  = Min(m_omegaMax,
                    Photons::s_reducemaxenergy*(m_M/2.)*( m_M/sum - sum/m_M ));
  if (m_omegaMax<0.) m_omegaMax = m_omegaMin;
  // set running alpha_QED to squared mass of incomming parton
  // -> taken at maximal scale
  Photons::SetAlphaQED(sqr(m_M));
}

Dipole_FF::~Dipole_FF() {
  DeleteAll(m_olddipole);
  DeleteAll(m_newdipole);
  DeleteAll(m_oldspectator);
  DeleteAll(m_newspectator);
}

void Dipole_FF::AddRadiation() {
  DEBUG_FUNC(ProcessName());
  DefineDipole();
  IdPairNbarVector nbars;
  // calculate avarage photon number in Lab frame
  if (PHOTONS::Photons::s_ircutoffframe==1) {
      Avarage_Photon_Number avnum(m_olddipole,m_omegaMax,m_omegaMin);
      m_nbar = avnum.GetNBar();
      nbars  = avnum.GetNBars();
  }
  Vec4D sumdip = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sumdip = sumdip + m_olddipole[i]->Momentum();
  }
  // boost into dipole CMS and rotate p_1 into z-axis
  Poincare boost(sumdip);
  boost.Boost(sumdip);
  Poincare rotate;
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D mom = m_olddipole[i]->Momentum();
    boost.Boost(mom);
    if (i==0)  rotate = Poincare(mom,Vec4D(0.,0.,0.,1.));
    rotate.Rotate(mom);
    m_olddipole[i]->SetMomentum(mom);
    m_Q = m_Q + mom;
  }
  // also transform neutral particles' momenta into Q-CMS
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D mom = m_oldspectator[i]->Momentum();
    boost.Boost(mom);
    rotate.Rotate(mom);
    m_oldspectator[i]->SetMomentum(mom);
    m_QN = m_QN + mom;
  }
  // calculate avarage photon number in Multipole_CMS
  double beta1(0.), beta2(0.);
  if (PHOTONS::Photons::s_ircutoffframe==0) {
    if (m_olddipole.size() == 2) {
      beta1 = CalculateBeta(m_olddipole[0]->Momentum());
      beta2 = CalculateBeta(m_olddipole[1]->Momentum());
      CalculateAvaragePhotonNumber(beta1,beta2);
    }
    else {
      Avarage_Photon_Number avnum(m_olddipole,m_omegaMax,m_omegaMin);
      m_nbar = avnum.GetNBar();
      nbars  = avnum.GetNBars();
    }
  }
  CheckAvaragePhotonNumberForNumericalErrors();
  msg_Debugging()<<"nbar: "<<m_nbar<<endl;
  // correction weights -> acception/rejection
  if (m_nbar > 0.) {
    bool genreject(true);
    while (genreject) {
      ResetVariables();
      // reject event if too many hard photons defy momentum conservation
      bool photoncheck(false);
      while (!photoncheck) {
        DeleteAllPhotons();
        if (m_olddipole.size() == 2) GeneratePhotons(beta1,beta2);
        else                         GeneratePhotons(nbars);
        m_K = CalculateMomentumSum(m_softphotons);
        DetermineKappa();
        if (m_n != 0)  photoncheck = CheckIfExceedingPhotonEnergyLimits();
        else photoncheck = true;
      }
      msg_Debugging()<<"n:    "<<m_n<<endl;
      if (m_n != 0) {
        CorrectMomenta();     // m_newdipole, m_newspectator, m_K are in P-CMS;
        if (m_u < 0.|| m_u > 1.) continue;
        CalculateWeights();
      }
      // optionally increase maximum
      m_genmaxweight *= Photons::s_increasemaxweight;
      if (ran->Get()*m_genmaxweight < m_genweight)   genreject = false;
      msg_Debugging()<<"-> "<<(genreject?"reject":"accept")<<std::endl;
      // accept new particle momenta if event accepted
      if (!genreject) {
      // if accepted rewrite momenta into Q-CMS,
      //   also transform the photon momenta
      // if no photons then P-CMS = Q-CMS
        if (m_n != 0) {
          // first boost to p-CMS (p = P + PN + K = Q + QN), then to Q-CMS
          Poincare boostPtop(m_P + m_PN + m_K);
          Poincare boostQtop(m_Q + m_QN);
          // now boost from P-CMS to Q-CMS via p-CMS
          for (unsigned int i=0; i<m_newdipole.size(); i++) {
            Vec4D mom = m_newdipole[i]->Momentum();
            boostPtop.Boost(mom);
            boostQtop.BoostBack(mom);
            m_newdipole[i]->SetMomentum(mom);
          }
          for (unsigned int i=0; i<m_newspectator.size(); i++) {
            Vec4D mom = m_newspectator[i]->Momentum();
            boostPtop.Boost(mom);
            boostQtop.BoostBack(mom);
            m_newspectator[i]->SetMomentum(mom);
          }
          for (unsigned int i=0; i<m_softphotons.size(); i++) {
            Vec4D k = m_softphotons[i]->Momentum();
            boostPtop.Boost(k);
            boostQtop.BoostBack(k);
            m_softphotons[i]->SetMomentum(k);
          }
          m_P   = CalculateMomentumSum(m_newdipole);
          m_PN  = CalculateMomentumSum(m_newspectator);
          m_K   = CalculateMomentumSum(m_softphotons);
        }
      }
    }
    // if no photon added
    if (m_n == 0) m_success = true;
    // if any photons added
    else {
      CheckMomentumConservationInQCMS(boost,rotate);
      // if momentum conserved, i.e. event successfully generated,
      // boost back to original system
      if (m_success) {
        m_photonsadded = true;
        for (unsigned int i=0; i<m_newdipole.size(); i++) {
          Vec4D mom = m_newdipole[i]->Momentum();
          rotate.RotateBack(mom);
          boost.BoostBack(mom);
          m_newdipole[i]->SetMomentum(mom);
        }
        for (unsigned int i=0; i<m_newspectator.size(); i++) {
          Vec4D mom = m_newspectator[i]->Momentum();
          rotate.RotateBack(mom);
          boost.BoostBack(mom);
          m_newspectator[i]->SetMomentum(mom);
        }
        for (unsigned int i=0; i<m_n; i++) {
          Vec4D k = m_softphotons[i]->Momentum();
          rotate.RotateBack(k);
          boost.BoostBack(k);
          m_softphotons[i]->SetMomentum(k);
        }
        ReturnMomenta();
      }
      else {
        DeleteAll(m_softphotons);
      }
    }
  }
}

void Dipole_FF::CalculateAvaragePhotonNumber(const double& b1,
                                             const double& b2) {
  double Z1      = m_olddipole[0]->Flav().Charge();
  double Z2      = m_olddipole[1]->Flav().Charge();
  m_nbar  = - Photons::s_alpha/M_PI * Z1*Z2 * log(m_omegaMax/m_omegaMin)
              *((1.+b1*b2)/(b1+b2)*log(((1.+b1)*(1.+b2))/((1.-b1)*(1.-b2)))-2.);
}

bool Dipole_FF::CheckIfExceedingPhotonEnergyLimits() {
  double sum = 0.;
  double nN = m_mN.size();
  double kappa2 = m_kappaN.Sqr();
  double Kmkappa2 = (Vec3D(m_K)-nN*m_kappaN).Sqr();
  for (unsigned int i=0; i<m_mC.size(); i++) sum+=m_mC[i];
  for (unsigned int i=0; i<m_mN.size(); i++) sum+=sqrt(sqr(m_mN[i])+kappa2);
  if (m_K[0] < (sqrt(m_M*m_M + Kmkappa2) - sum)) return true;
  else                                           return false;
}

void Dipole_FF::CheckMomentumConservationInQCMS
(const Poincare& boost,const Poincare& rotate) {
  // in Q-CMS
  if ((m_u > 1) || (m_u < 0)) {
    msg_Error()<<"u: "<<m_u<<" not in [0,1] ... all photons deleted ..."<<endl;
    m_success = false;
    return;
  }
  m_success = false;
  Vec4D p(m_neutralinparticles[0]->Momentum());
  boost.Boost(p);
  rotate.Rotate(p);
  // after radiation  in Q-CMS
  Vec4D P   = Vec4D(0.,0.,0.,0.);
  Vec4D PN  = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_newdipole.size(); i++)
    P += m_newdipole[i]->Momentum();
  for (unsigned int i=0; i<m_newspectator.size(); i++)
    PN += m_newspectator[i]->Momentum();
  // difference
  Vec4D diff  = p - P - PN - m_K;
  double accu = m_accu*p[0];
  if ((abs(diff[0]) < accu) && (abs(diff[1]) < accu) &&
      (abs(diff[2]) < accu) && (abs(diff[3]) < accu))
    m_success = true;
  else {
    msg_Error()<<"momentum not conserved! residual is: "<<diff
               <<" accuracy is: "<<accu<<endl;
    for (unsigned int i=0;i<m_olddipole.size();i++)
      msg_Debugging()<<*m_olddipole[i]<<endl;
    for (unsigned int i=0;i<m_oldspectator.size();i++)
      msg_Debugging()<<*m_oldspectator[i]<<endl;
    for (unsigned int i=0;i<m_newdipole.size();i++)
      msg_Debugging()<<*m_newdipole[i]<<endl;
    for (unsigned int i=0;i<m_newspectator.size();i++)
      msg_Debugging()<<*m_newspectator[i]<<endl;
    for (unsigned int i=0;i<m_n;i++)
      msg_Debugging()<<*m_softphotons[i]<<endl;
    msg_Error()<<"all photons deleted..."<<endl;
  }
}

void Dipole_FF::CorrectMomenta() {
  DetermineU();
  // reconstruct momenta in P-CMS
  if ((m_u >= 0.) && (m_u <= 1.)) {
    // charged final state momenta
    for (unsigned int i=0; i<m_olddipole.size(); i++) {
      Vec3D mom = m_u *Vec3D(m_olddipole[i]->Momentum());
      double mom0 = sqrt(sqr(m_mC[i])+mom.Sqr());
      m_newdipole[i]->SetMomentum(Vec4D(mom0,mom));
      m_P += m_newdipole[i]->Momentum();
    }
    // neutral final state momenta
    for (unsigned int i=0; i<m_oldspectator.size(); i++) {
      Vec3D mom = m_u *Vec3D(m_oldspectator[i]->Momentum())-m_kappaN;
      double mom0 = sqrt(sqr(m_mN[i])+mom.Sqr());
      m_newspectator[i]->SetMomentum(Vec4D(mom0,mom));
      m_PN += m_newspectator[i]->Momentum();
    }
  }
}

void Dipole_FF::DefineDipole() {
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++) {
    m_olddipole.push_back(new Particle(*m_chargedoutparticles[i]));
    m_olddipole[i]->SetProductionBlob(m_chargedoutparticles[i]->ProductionBlob());
    m_olddipole[i]->SetDecayBlob(m_chargedoutparticles[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++) {
    m_oldspectator.push_back(new Particle(*m_neutraloutparticles[i]));
    m_oldspectator[i]->SetProductionBlob(m_neutraloutparticles[i]->ProductionBlob());
    m_oldspectator[i]->SetDecayBlob(m_neutraloutparticles[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    m_newdipole.push_back(new Particle(*m_olddipole[i]));
    m_newdipole[i]->SetProductionBlob(m_olddipole[i]->ProductionBlob());
    m_newdipole[i]->SetDecayBlob(m_olddipole[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    m_newspectator.push_back(new Particle(*m_oldspectator[i]));
    m_newspectator[i]->SetProductionBlob(m_oldspectator[i]->ProductionBlob());
    m_newspectator[i]->SetDecayBlob(m_oldspectator[i]->DecayBlob());
  }
}

double Dipole_FF::Func(const double& M2, const std::vector<double>& mC2,
                       const std::vector<double>& mN2,
                       const std::vector<Vec3D>& q, const double& u) {
  double sum = 0.;
  double nN = m_mN.size();
  for (unsigned int i=0; i<mC2.size(); i++)
    sum+=sqrt(mC2[i]+u*u*(q[i]*q[i]));
  for (unsigned int i=0; i<mN2.size(); i++)
    sum+=sqrt(mN2[i]+(u*q[i+mC2.size()]-m_kappaN).Sqr());
  return sqrt(M2+(u*Vec3D(m_QN)+Vec3D(m_K)-nN*m_kappaN).Sqr()) - m_K[0] - sum;
}

void Dipole_FF::ResetVariables() {
  DeleteAllPhotons();
  for (unsigned int i=0; i<m_olddipole.size(); i++)
    m_newdipole[i]->SetMomentum(m_olddipole[i]->Momentum());
  for (unsigned int i=0; i<m_oldspectator.size(); i++)
    m_newspectator[i]->SetMomentum(m_oldspectator[i]->Momentum());
  m_u   = 1.;
  m_K   = Vec4D(0.,0.,0.,0.);
  m_P   = Vec4D(0.,0.,0.,0.);
  m_PN  = Vec4D(0.,0.,0.,0.);
  m_kappaC = Vec3D(0.,0.,0.);
  m_kappaN = Vec3D(0.,0.,0.);
  m_genweight     = 1.;
  m_genmaxweight  = 1.;
}

void Dipole_FF::ReturnMomenta() {
  for(unsigned int i=0; i<m_newdipole.size(); i++)
    m_chargedoutparticles[i]->SetMomentum(m_newdipole[i]->Momentum());
  for(unsigned int i=0; i<m_newspectator.size(); i++)
    m_neutraloutparticles[i]->SetMomentum(m_newspectator[i]->Momentum());
}

void Dipole_FF::DetermineKappa() {
  double nN = m_mN.size();
  m_kappaC = Vec3D(0.,0.,0.);
  switch (Photons::s_ffrecscheme) {
  case 1:
    m_kappaN = 1./(nN+1.)*Vec3D(m_K);
    break;
  case 2:
    if (nN>0) {
      m_kappaN = 1./nN*Vec3D(m_K);
      break;
    }
  default:
    m_kappaN = Vec3D(0.,0.,0.);
    break;
  }
}
