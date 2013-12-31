#include "PHOTONS++/Tools/Dipole_FI.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Particle.H"
#include "MODEL/Main/Model_Base.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/PhaseSpace/Avarage_Photon_Number.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Dipole_FI::Dipole_FI(const Particle_Vector_Vector& pvv) {
  m_pvv                 = pvv;
  m_dtype               = Dipole_Type::fi;
  m_chargedinparticles  = pvv[0];
  m_neutralinparticles  = pvv[1];
  m_chargedoutparticles = pvv[2];
  m_neutraloutparticles = pvv[3];
  m_M                   = m_chargedinparticles[0]->Momentum().Mass();
  m_Q                   = Vec4D(0.,0.,0.,0.);
  m_QN                  = Vec4D(0.,0.,0.,0.);
  m_kappaC              = Vec3D(0.,0.,0.);
  m_kappaN              = Vec3D(0.,0.,0.);
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++)
    m_mC.push_back(m_chargedoutparticles[i]->FinalMass());
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++)
    m_mN.push_back(m_neutraloutparticles[i]->FinalMass());
  m_omegaMax = Min(m_omegaMax,
                   Photons::s_reducemax * DetermineMaximumPhotonEnergy());
  // set running alpha_QED to squared mass of incomming parton
  // -> taken at maximal scale
  Photons::SetAlphaQED(sqr(m_M));
}

Dipole_FI::~Dipole_FI() {
  DeleteAll(m_olddipole);
  DeleteAll(m_newdipole);
  DeleteAll(m_oldspectator);
  DeleteAll(m_newspectator);
}

void Dipole_FI::AddRadiation() {
  DEBUG_FUNC(ProcessName());
  DefineDipole();
  Vec4D sumdip = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sumdip = sumdip + m_olddipole[i]->Momentum();
  }
  // boost into dipole CMS and initial state particle into -z
  Poincare boost(sumdip);
  boost.Boost(sumdip);
  Poincare rotate;
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D mom = m_olddipole[i]->Momentum();
    boost.Boost(mom);
    if (i==0)  rotate = Poincare(mom/mom.PSpat2(),Vec4D(0.,0.,0.,-1.));
    rotate.Rotate(mom);
    m_olddipole[i]->SetMomentum(mom);
    if (i==0)  m_p = m_olddipole[0]->Momentum();
    else       m_Q = m_Q + m_olddipole[i]->Momentum();
  }
  // also transform neutral particles' momenta into (p+Q)-CMS
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D mom = m_oldspectator[i]->Momentum();
    boost.Boost(mom);
    rotate.Rotate(mom);
    m_oldspectator[i]->SetMomentum(mom);
    m_QN = m_QN + mom;
  }
  // calculate avarage photon number (only for two particle dipole)
  double beta1(0.), beta2(0.);
  IdPairNbarVector nbars;
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
  CheckAvaragePhotonNumberForNumericalErrors();
  msg_Debugging()<<"nbar: "<<m_nbar<<endl;
  // correction weights -> acception/rejection
  if (m_nbar > 0) {
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
        CorrectMomenta();
        if (m_u < 0.|| m_u > 1.) continue;
        // m_newdipole, m_newspectator, m_K, m_P and m_PN are in (p+P)-CMS
        // m_plddipole, m_oldspectator. m_Q and m_QN are in (p+Q)-CMS
        CalculateWeights();
      }
      // if no photon generated, event always accepted
      else break;
      if (ran->Get()*m_genmaxweight < m_genweight)   genreject = false;
      msg_Debugging()<<"-> "<<(genreject?"reject":"accept")<<std::endl;
      // accept new particle momenta if event accepted
      if (!genreject) {
      // if accepted rewrite momenta into (p+Q)-CMS, also transform the photon momenta
      // if no photons then (p+P)-CMS = (p+Q)-CMS
        if (m_n != 0) {
          // first boost to p-CMS (p = P + PN + K = Q + QN), then to (p+Q)-CMS
          Poincare boostPtop(m_P + m_PN + m_K);
          Poincare boostQtop(m_Q + m_QN);
          //now boost from (p+P)-CMS to (p+Q)-CMS via p-CMS
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
        }
      }
    }
    // if no photon added
    if (m_n == 0) m_success = true;
    // if any photons added
    if (m_n != 0) {
      CheckMomentumConservationInQCMS(boost,rotate);
      // if momentum conserved, i.e. event successfully generated, boost back to original system
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

void Dipole_FI::CalculateAvaragePhotonNumber(const double& b1,
                                             const double& b2) {
  double Z1      = m_olddipole[0]->Flav().Charge();
  double Z2      = m_olddipole[1]->Flav().Charge();
  m_nbar  = Photons::s_alpha/M_PI * Z1*Z2 * log(m_omegaMax/m_omegaMin)
            *((1.+b1*b2)/(b1+b2)*log(((1.+b1)*(1.+b2))/((1.-b1)*(1.-b2))) - 2.);
}

bool Dipole_FI::CheckIfExceedingPhotonEnergyLimits() {
  double sum = 0.;
  double nC  = m_mC.size();
  double kappaC2 = m_kappaC.Sqr();
  double kappaN2 = m_kappaN.Sqr();
  for (unsigned int i=0; i<m_mC.size(); i++) sum+=sqrt(m_mC[i]*m_mC[i]+kappaC2);
  for (unsigned int i=0; i<m_mN.size(); i++) sum+=sqrt(m_mN[i]*m_mN[i]+kappaN2);
  if (m_K[0] < sqrt(m_M*m_M + nC*nC*kappaC2) - sum) return true;
  else                                             return false;
}

void Dipole_FI::CheckMomentumConservationInQCMS
(const Poincare& boost,const Poincare& rotate) {
  // in (p+Q)-CMS
  if ((m_u > 1) || (m_u < 0)) {
    msg_Error()<<"u: "<<m_u<<" not in [0,1] ... all photons deleted ..."<<endl;
    m_success = false;
    return;
  }
  m_success = false;
  Vec4D p   = m_newdipole[0]->Momentum();
  // after radiation in Q-CMS
  m_P = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=1; i<m_newdipole.size(); i++)
    m_P = m_P + m_newdipole[i]->Momentum();
  m_PN  = CalculateMomentumSum(m_newspectator);
  m_K   = CalculateMomentumSum(m_softphotons);
  // difference
  Vec4D diff  = p - m_P - m_PN - m_K;
  double accu = m_accu*m_newdipole[0]->Momentum()[0];
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

void Dipole_FI::CorrectMomenta() {
  DetermineU();
  double nC = m_mC.size();
  if ((m_u >= 0.) && (m_u <= 1.)) {
    // initial state charged momentum
    Vec3D p = (-m_u)*Vec3D(m_Q)+nC*m_kappaC;
    m_newdipole[0]->SetMomentum(Vec4D(sqrt(m_M*m_M + p*p), p));
    // final state charged momenta
    for (unsigned int i=1; i<m_newdipole.size(); i++) {
      Vec3D q = m_u*Vec3D(m_olddipole[i]->Momentum())-m_kappaC;
      m_newdipole[i]->SetMomentum(Vec4D(sqrt(m_mC[i-1]*m_mC[i-1] + q*q), q));
      m_P = m_P + m_newdipole[i]->Momentum();
    }
    // final state neutral momenta
    for (unsigned int i=0; i<m_newspectator.size(); i++) {
      Vec3D q = m_u*Vec3D(m_oldspectator[i]->Momentum())-m_kappaN;
      m_newspectator[i]->SetMomentum(Vec4D(sqrt(m_mN[i]*m_mN[i] + q*q), q));
      m_PN = m_PN + m_newspectator[i]->Momentum();
    }
  }
}

void Dipole_FI::DefineDipole() {
    m_olddipole.push_back(new Particle(*m_chargedinparticles[0]));
    m_olddipole.at(0)->SetProductionBlob(m_chargedinparticles[0]->ProductionBlob());
    m_olddipole.at(0)->SetDecayBlob(m_chargedinparticles[0]->DecayBlob());
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++) {
    m_olddipole.push_back(new Particle(*m_chargedoutparticles[i]));
    m_olddipole[i+1]->SetProductionBlob(m_chargedoutparticles[i]->ProductionBlob());
    m_olddipole[i+1]->SetDecayBlob(m_chargedoutparticles[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++) {
    m_oldspectator.push_back(new Particle(*m_neutraloutparticles[i]));
    m_oldspectator[i]->SetProductionBlob(m_neutraloutparticles[i]->ProductionBlob());
    m_oldspectator[i]->SetDecayBlob(m_neutraloutparticles[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    m_newdipole.push_back(new Particle(*m_olddipole[i]));
    m_newdipole.at(i)->SetProductionBlob(m_olddipole[i]->ProductionBlob());
    m_newdipole.at(i)->SetDecayBlob(m_olddipole[i]->DecayBlob());
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    m_newspectator.push_back(new Particle(*m_oldspectator[i]));
    m_newspectator[i]->SetProductionBlob(m_oldspectator[i]->ProductionBlob());
    m_newspectator[i]->SetDecayBlob(m_oldspectator[i]->DecayBlob());
  }
}

double Dipole_FI::DetermineMaximumPhotonEnergy() {
  double sum(0.);
  unsigned int nC(m_mC.size()), nN(m_mN.size()), n(nC+nN);
  std::vector<double> mCN2;
  for (unsigned int i=0; i<nC; ++i) {
    sum+=m_mC[i];
    mCN2.push_back(sqr(m_mC[i]));
  }
  for (unsigned int i=0; i<nN; ++i) {
    sum+=m_mN[i];
    mCN2.push_back(sqr(m_mN[i]));
  }
  if (n != m_mC.size()+m_mN.size()) {
    msg_Out()<<METHOD<<"error while determining maximum photon energy\n";
    return m_omegaMin;
  }
  unsigned int count = 0;
  double fac         = 1./(2.*nC+nN);
  double x0          = 0.5*(m_M-sum);
  double xNminus1    = x0;
  double xN          = 0.;
  double F_xNminus1  = 0.;
  while (abs(xN - xNminus1) > 1E-6) {
    if (count==500) {
      msg_Out()<<"failed to determine maximum photon energy... set to IR cut-off..."<<endl;
      return m_omegaMin;
    }
    xNminus1 = xN;
    F_xNminus1 = sqr(fac*xNminus1);
    sum = 0.;
    for (unsigned int i=0; i<n; ++i)  sum += sqrt(mCN2[i]+F_xNminus1);
    xN = sqrt(m_M*m_M + nC*nC*F_xNminus1) - sum;
    count++;
  }
  if (xN<0.) return m_omegaMin;
  return xN;
}

double Dipole_FI::Func(const double& M2, const std::vector<double>& mC2,
                       const std::vector<double>& mN2,
                       const std::vector<Vec3D>& q, const double& u) {
  double sum  = 0.;
  int    nC   = m_mC.size();
  // first entry in q belongs to initial state dipole constituent
  for (unsigned int i=0; i<mC2.size(); i++)
    sum+=sqrt(mC2[i]+(u*q[1+i]-m_kappaC).Sqr());
  for (unsigned int i=0; i<mN2.size(); i++)
    sum+=sqrt(mN2[i]+(u*q[1+nC+i]-m_kappaN).Sqr());
  return sqrt(M2+(u*Vec3D(m_Q)-nC*m_kappaC).Sqr()) - sum - m_K[0];
}

void Dipole_FI::ResetVariables() {
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

void Dipole_FI::ReturnMomenta() {
  for(unsigned int i=1; i<m_newdipole.size(); i++)
    m_chargedoutparticles[i-1]->SetMomentum(m_newdipole[i]->Momentum());
  for(unsigned int i=0; i<m_newspectator.size(); i++)
    m_neutraloutparticles[i]->SetMomentum(m_newspectator[i]->Momentum());
}

void Dipole_FI::DetermineKappa() {
  int nC      = m_mC.size();
  int nN      = m_mN.size();
  switch (Photons::s_firecscheme) {
  case 1:
    m_kappaN = Vec3D(0.,0.,0.);
    m_kappaC = 1./(2.*nC)*Vec3D(m_K);
  case 2:
    if (nN>0) {
      m_kappaN = 1./nN*Vec3D(m_K);
      m_kappaC = Vec3D(0.,0.,0.);
      break;
    }
  default:
    m_kappaC = m_kappaN = 1./(2.*nC+nN)*Vec3D(m_K);
    break;
  }
}
