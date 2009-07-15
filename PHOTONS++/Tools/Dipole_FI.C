#include "PHOTONS++/Tools/Dipole_FI.H"
#include "PHOTONS++/Tools/Generate_One_Photon.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Dipole_FI::Dipole_FI(Particle_Vector_Vector pvv) {
  m_pvv                 = pvv;
  m_dtype               = Dipole_Type::fi;
  m_chargedinparticles  = pvv.at(0);
  m_neutralinparticles  = pvv.at(1);
  m_chargedoutparticles = pvv.at(2);
  m_neutraloutparticles = pvv.at(3);
  m_M                   = m_chargedinparticles.at(0)->FinalMass();
  m_Q                   = Vec4D(0.,0.,0.,0.);
  m_QN                  = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++) {
    if (m_chargedoutparticles.at(i)->FinalMass() == 0.105)
      m_mC.push_back(m_chargedoutparticles.at(i)->FinalMass() - 1E-9);
    else
      m_mC.push_back(m_chargedoutparticles.at(i)->FinalMass());
  }
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++) {
    m_mN.push_back(m_neutraloutparticles.at(i)->FinalMass());
  }
  m_omegaMax = DetermineMaximumPhotonEnergy();
}

Dipole_FI::~Dipole_FI() {
  DeleteAll(m_olddipole);
  DeleteAll(m_newdipole);
  DeleteAll(m_oldspectator);
  DeleteAll(m_newspectator);
}

void Dipole_FI::AddRadiation() {
  DefineDipole();
  Vec4D sumdip = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sumdip = sumdip + m_olddipole.at(i)->Momentum();
  }
  // boost into dipole CMS and initial state particle into -z
  Poincare boost(sumdip);
  boost.Boost(sumdip);
  Poincare rotate;
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D mom = m_olddipole.at(i)->Momentum();
    boost.Boost(mom);
    if (i==0)  rotate = Poincare(mom,Vec4D(0.,0.,0.,-1.));
    rotate.Rotate(mom);
    m_olddipole.at(i)->SetMomentum(mom);
    if (i==0)  m_p = m_olddipole.at(0)->Momentum();
    else       m_Q = m_Q + m_olddipole.at(i)->Momentum();
  }
  // also transform neutral particles' momenta into (p+Q)-CMS
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D mom = m_oldspectator.at(i)->Momentum();
    boost.Boost(mom);
    rotate.Rotate(mom);
    m_oldspectator.at(i)->SetMomentum(mom);
    m_QN = m_QN + mom;
  }
  // calculate avarage photon number (only for two particle dipole)
  double beta1 = 0;
  double beta2 = 0;
  std::vector<double> nbars;
  if (m_olddipole.size() == 2) {
    beta1 = CalculateBeta(m_olddipole.at(0)->Momentum());
    beta2 = CalculateBeta(m_olddipole.at(1)->Momentum());
    CalculateAvaragePhotonNumber(beta1,beta2);
  }
  else {
//     msg_Info()<<"multipole calculation"<<endl;
    Avarage_Photon_Number avnum(m_olddipole,m_omegaMax,m_omegaMin);
    m_nbar = avnum.GetNBar();
    nbars  = avnum.GetNBars();
  }
  CheckAvaragePhotonNumberForNumericalErrors();
#ifdef PHOTONS_DEBUG
  msg_Info()<<"nbar: "<<m_nbar<<endl;
#endif
  // correction weights -> acception/rejection
  if (m_nbar > 0) {
    bool genreject = true;
    while (genreject == true) {
      ResetVariables();
      // reject event if too many hard photons defy momentum conservation
      bool photoncheck = false;
      while (photoncheck == false) {
        DeleteAllPhotons();
        if (m_olddipole.size() == 2)  GeneratePhotons(beta1,beta2);
        else                          GeneratePhotons(nbars);
        m_K = CalculateMomentumSum(m_softphotons);
        if (m_n != 0)  photoncheck = CheckIfExceedingPhotonEnergyLimits();
        else photoncheck = true;
      }
      if (m_n != 0) {
        CorrectMomenta();
        // m_newdipole, m_newspectator, m_K, m_P and m_PN are in (p+P)-CMS
        // m_plddipole, m_oldspectator. m_Q and m_QN are in (p+Q)-CMS
        CalculateWeights();
      }
      if (ran.Get()*m_genmaxweight < m_genweight)   genreject = false;
      // accept new particle momenta if event accepted
      if (genreject == false) {
      // if accepted rewrite momenta into (p+Q)-CMS, also transform the photon momenta
      // if no photons then (p+P)-CMS = (p+Q)-CMS
        if (m_n > 0) {
          // first boost to p-CMS (p = P + PN + K = Q + QN), then to (p+Q)-CMS
          Poincare boostPtop(m_P + m_PN + m_K);
          Poincare boostQtop(m_Q + m_QN);
          //now boost from (p+P)-CMS to (p+Q)-CMS via p-CMS
          for (unsigned int i=0; i<m_newdipole.size(); i++) {
            Vec4D mom = m_newdipole.at(i)->Momentum();
            boostPtop.Boost(mom);
            boostQtop.BoostBack(mom);
            m_newdipole.at(i)->SetMomentum(mom);
          }
          for (unsigned int i=0; i<m_newspectator.size(); i++) {
            Vec4D mom = m_newspectator.at(i)->Momentum();
            boostPtop.Boost(mom);
            boostQtop.BoostBack(mom);
            m_newspectator.at(i)->SetMomentum(mom);
          }
          for (unsigned int i=0; i<m_softphotons.size(); i++) {
            Vec4D k = m_softphotons.at(i)->Momentum();
            boostPtop.Boost(k);
            boostQtop.BoostBack(k);
            m_softphotons.at(i)->SetMomentum(k);
          }
        }
      }
    }
    // if no photon added
    if (m_n == 0) m_success = true;
    // if any photons added
    if (m_n > 0) {
      CheckMomentumConservationInQCMS();
      // if momentum conserved, i.e. event successfully generated, boost back to original system
      if (m_success == true) {
        m_photonsadded = true;
        for (unsigned int i=0; i<m_newdipole.size(); i++) {
          Vec4D mom = m_newdipole.at(i)->Momentum();
          rotate.RotateBack(mom);
          boost.BoostBack(mom);
          m_newdipole.at(i)->SetMomentum(mom);
        }
        for (unsigned int i=0; i<m_newspectator.size(); i++) {
          Vec4D mom = m_newspectator.at(i)->Momentum();
          rotate.RotateBack(mom);
          boost.BoostBack(mom);
          m_newspectator.at(i)->SetMomentum(mom);
        }
        for (unsigned int i=0; i<m_n; i++) {
          Vec4D k = m_softphotons.at(i)->Momentum();
          rotate.RotateBack(k);
          boost.BoostBack(k);
          m_softphotons.at(i)->SetMomentum(k);
        }
        ReturnMomenta();
      }
      else {
        DeleteAll(m_softphotons);
      }
    }
  }
}

void Dipole_FI::CalculateAvaragePhotonNumber(const double b1, const double b2) {
  double alpha   = MODEL::s_model->ScalarConstant("alpha_QED(0)");
  double Z1      = m_olddipole.at(0)->Flav().Charge();
  double Z2      = m_olddipole.at(1)->Flav().Charge();
  m_nbar  = alpha/M_PI * Z1*Z2 * log(m_omegaMax/m_omegaMin) * ((1+b1*b2)/(b1+b2)*log(((1+b1)*(1+b2))/((1-b1)*(1-b2))) - 2);
}

bool Dipole_FI::CheckIfExceedingPhotonEnergyLimits() {
  double sum = 0;
  double nC   = m_mC.size();
  double nN   = m_mN.size();
  for (unsigned int i=0; i<m_mC.size(); i++) {
    sum = sum + sqrt(m_mC.at(i)*m_mC.at(i)+(1/(2*nC+nN)*Vec3D(m_K))*(1/(2*nC+nN)*Vec3D(m_K)));
  }
  for (unsigned int i=0; i<m_mN.size(); i++) {
    sum = sum + sqrt(m_mN.at(i)*m_mN.at(i)+(1/(2*nC+nN)*Vec3D(m_K))*(1/(2*nC+nN)*Vec3D(m_K)));
  }
  if (m_K[0] < sqrt(m_M*m_M + (nC/(2*nC+nN)*Vec3D(m_K))*(nC/(2*nC+nN)*Vec3D(m_K))) - sum)   return true;
  else {
    return false;
  }
}

void Dipole_FI::CheckMomentumConservationInQCMS() {
  // in (p+Q)-CMS
  if ((m_u > 1) || (m_u < 0)) msg_Out()<<"u: "<<m_u<<endl;
  m_success = false;
  Vec4D p   = m_newdipole.at(0)->Momentum();
  // after radiation in Q-CMS
  m_P = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=1; i<m_newdipole.size(); i++) {
    m_P = m_P + m_newdipole.at(i)->Momentum();
  }
  m_PN  = CalculateMomentumSum(m_newspectator);
  m_K   = CalculateMomentumSum(m_softphotons);
  // difference
  Vec4D diff  = p - m_P - m_PN - m_K;
  double accu = m_accu*m_newdipole.at(0)->Momentum()[0];
  if ((abs(diff[0]) < accu) && (abs(diff[1]) < accu) &&
      (abs(diff[2]) < accu) && (abs(diff[3]) < accu))
    m_success = true;
  else {
    msg_Out()<<"momentum not conserved! residual is: "<<diff<<" accuracy is: "<<accu<<endl;
    for (unsigned int i=0;i<m_olddipole.size();i++) {
      msg_Out()<<*m_olddipole.at(i)<<endl;
    }
    for (unsigned int i=0;i<m_oldspectator.size();i++) {
      msg_Out()<<*m_oldspectator.at(i)<<endl;
    }
    for (unsigned int i=0;i<m_newdipole.size();i++) {
      msg_Out()<<*m_newdipole.at(i)<<endl;
    }
    for (unsigned int i=0;i<m_newspectator.size();i++) {
      msg_Out()<<*m_newspectator.at(i)<<endl;
    }
    for (unsigned int i=0;i<m_n;i++) {
      msg_Out()<<*m_softphotons.at(i)<<endl;
    }
  }
}

void Dipole_FI::CorrectMomenta() {
  DetermineU();
  double nC = m_mC.size();
  double nN = m_mN.size();
  if ((m_u >= 0) && (m_u <= 1)) {
    // initial state charged momentum
    m_newdipole.at(0)->SetMomentum(Vec4D(sqrt(m_M*m_M + (m_u*Vec3D(m_Q)-(nC/(2*nC+nN)*Vec3D(m_K)))*(m_u*Vec3D(m_Q)-(nC/(2*nC+nN)*Vec3D(m_K)))), (-m_u)*Vec3D(m_Q)+(nC/(2*nC+nN)*Vec3D(m_K))));
    // final state charged momenta
    for (unsigned int i=1; i<m_newdipole.size(); i++) {
      m_newdipole.at(i)->SetMomentum(Vec4D(sqrt(m_mC.at(i-1)*m_mC.at(i-1) + (m_u*Vec3D(m_olddipole.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K)))*(m_u*Vec3D(m_olddipole.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K)))), m_u*Vec3D(m_olddipole.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K))));
      m_P = m_P + m_newdipole.at(i)->Momentum();
    }
    // final state neutral momenta
    double nN = m_newspectator.size();
    for (unsigned int i=0; i<m_newspectator.size(); i++) {
      m_newspectator.at(i)->SetMomentum(
          Vec4D(sqrt(m_mN.at(i)*m_mN.at(i) + 
                   (m_u*Vec3D(m_oldspectator.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K)))
                   *(m_u*Vec3D(m_oldspectator.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K)))),
                m_u*Vec3D(m_oldspectator.at(i)->Momentum())-(1./(2*nC+nN)*Vec3D(m_K))));
      m_PN = m_PN + m_newspectator.at(i)->Momentum();
    }
  }
}

void Dipole_FI::DefineDipole() {
    m_olddipole.push_back(new Particle(*m_chargedinparticles.at(0)));
    m_olddipole.at(0)->SetProductionBlob(m_chargedinparticles.at(0)->ProductionBlob());
    m_olddipole.at(0)->SetDecayBlob(m_chargedinparticles.at(0)->DecayBlob());
  for (unsigned int i=0; i<m_chargedoutparticles.size(); i++) {
    m_olddipole.push_back(new Particle(*m_chargedoutparticles.at(i)));
    m_olddipole.at(i+1)->SetProductionBlob(m_chargedoutparticles.at(i)->ProductionBlob());
    m_olddipole.at(i+1)->SetDecayBlob(m_chargedoutparticles.at(i)->DecayBlob());
  }
  for (unsigned int i=0; i<m_neutraloutparticles.size(); i++) {
    m_oldspectator.push_back(new Particle(*m_neutraloutparticles.at(i)));
    m_oldspectator.at(i)->SetProductionBlob(m_neutraloutparticles.at(i)->ProductionBlob());
    m_oldspectator.at(i)->SetDecayBlob(m_neutraloutparticles.at(i)->DecayBlob());
  }
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    m_newdipole.push_back(new Particle(*m_olddipole.at(i)));
    m_newdipole.at(i)->SetProductionBlob(m_olddipole.at(i)->ProductionBlob());
    m_newdipole.at(i)->SetDecayBlob(m_olddipole.at(i)->DecayBlob());
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    m_newspectator.push_back(new Particle(*m_oldspectator.at(i)));
    m_newspectator.at(i)->SetProductionBlob(m_oldspectator.at(i)->ProductionBlob());
    m_newspectator.at(i)->SetDecayBlob(m_oldspectator.at(i)->DecayBlob());
  }
}

double Dipole_FI::DetermineMaximumPhotonEnergy() {
  int    i        = 0;
  double sum      = 0;
  for (unsigned int i=0; i<m_mC.size(); i++) {
    sum = sum + m_mC.at(i);
  }
  for (unsigned int i=0; i<m_mN.size(); i++) {
    sum = sum + m_mN.at(i);
  }
  double x0       = (1./2.)*(m_M-sum);
  double xNminus1 = x0;
  double xN       = 0;
  double nC       = m_mC.size();
  double nN       = m_mN.size();
  while (abs(xN - xNminus1) > 1E-6) {
    if (i != 0) xNminus1 = xN;
    sum = 0;
    for (unsigned int j=0; j<nC; j++) {
      sum = sum + sqrt(m_mC.at(j)*m_mC.at(j) + pow(1./(2*nC+nN)*xNminus1,2));
    }
    for (unsigned int j=0; j<nN; j++) {
      sum = sum + sqrt(m_mN.at(j)*m_mN.at(j) + pow(1./(2*nC+nN)*xNminus1,2));
    }
    xN = sqrt(m_M*m_M + pow(nC/(2*nC+nN)*xNminus1,2)) - sum;
    i++;
    if (i==500) {
      msg_Out()<<"failed to determine maximum photon energy... set to IR cut-off..."<<endl;
      return m_omegaMin;
    }
  }
  return xN;
}

double Dipole_FI::Func(double M2, std::vector<double> mC2, std::vector<double> mN2, std::vector<Vec3D> q, double u) {
    double sum  = 0;
    int nC      = m_mC.size();
    int nN      = m_mN.size();
    // first entry in q belongs to initial state dipole constituent
    for (unsigned int i=0; i<mC2.size(); i++) {
      sum = sum + sqrt(mC2.at(i) + (u*q.at(1+i)-1./(2.*nC+nN)*Vec3D(m_K))*(u*q.at(1+i)-1./(2.*nC+nN)*Vec3D(m_K)));
    }
    for (unsigned int i=0; i<mN2.size(); i++) {
      sum = sum + sqrt(mN2.at(i) + (u*q.at(1+nC+i)-1./(2.*nC+nN)*Vec3D(m_K))*(u*q.at(1+nC+i)-1./(2.*nC+nN)*Vec3D(m_K)));
    }
    double r = sqrt(M2+(u*Vec3D(m_Q)-nC/(2.*nC+nN)*Vec3D(m_K))*(u*Vec3D(m_Q)-nC/(2.*nC+nN)*Vec3D(m_K))) - sum - m_K[0];
    return r;
}

void Dipole_FI::ResetVariables() {
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    m_newdipole.at(i)->SetMomentum(m_olddipole.at(i)->Momentum());
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    m_newspectator.at(i)->SetMomentum(m_oldspectator.at(i)->Momentum());
  }
  m_u   = 1;
  m_K   = Vec4D(0.,0.,0.,0.);
  m_P   = Vec4D(0.,0.,0.,0.);
  m_PN  = Vec4D(0.,0.,0.,0.);
  m_genweight     = 1;
  m_genmaxweight  = 1;
}

void Dipole_FI::ReturnMomenta() {
  for(unsigned int i=1; i<m_newdipole.size(); i++) {
    m_chargedoutparticles.at(i-1)->SetMomentum(m_newdipole.at(i)->Momentum());
  }
  for(unsigned int i=0; i<m_newspectator.size(); i++) {
    m_neutraloutparticles.at(i)->SetMomentum(m_newspectator.at(i)->Momentum());
  }
}
