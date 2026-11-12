#include "YFS/Main/Dipole.H"
#include <cstdlib>
#include "YFS/Main/FSR.H"
#include "YFS/Main/ISR.H"
#include "YFS/Main/Emission.H"

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H" 
#include "PHASIC++/Channels/Channel_Elements.H"



using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;


double delf = 0;
double deli = 0;



Dipole::Dipole(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, 
              ATOOLS::Vec4D_Vector const &born, dipoletype::code ty, const double alpha):
  m_type(ty), m_alp(alpha)
{
  if ((mom.size() != fl.size()) || fl.size() != 2 || mom.size() != 2 || born.size()!=2) {
    msg_Out()<<"Dipole type is  = "<<ty<<std::endl
             <<" mom.size() = "<<mom.size()<<std::endl
             <<" fl.size() = "<<fl.size()<<std::endl
             <<" born.size() = "<<born.size()<<std::endl
             <<"Flavours = "<<fl<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for dipoletype");
  }
  Clean();
  // todo get alpha from YFS_BASE
  m_irfinite = false;
  m_alpi = m_alp/M_PI;
  m_sp = (mom[0]+mom[1]).Abs2();
  m_Qi = fl[0].Charge();
  m_Qj = fl[1].Charge();
  // if(fl[0].IsBoson() || fl[1].IsBoson()) m_irfinite = true;// Case for on shell ww
  m_QiQj = m_Qi*m_Qj;
  // Direct 2-element init avoids push_back reallocations (built every event).
  m_flavs = {fl[0], fl[1]};
  m_masses = {fl[0].Mass(), fl[1].Mass()};
  m_charges = {fl[0].Charge(), fl[1].Charge()};
  m_names = {fl[0].IDName(), fl[1].IDName()};
  m_momenta = {mom[0], mom[1]};
  m_oldmomenta = m_momenta;
  m_newmomenta = m_momenta;
  // Final-type dipoles never read m_ghost, so skip filling it.
  if (ty != dipoletype::code::final) m_ghost = m_momenta;
  m_bornmomenta = {born[0], born[1]};
  m_eikmomentum = m_bornmomenta;
  if (ty == dipoletype::code::initial) {
    m_thetai = m_thetaj = 1;
  }
  else if (ty == dipoletype::code::final) {
    m_thetai = m_thetaj = -1;
  }
  else if (ty == dipoletype::code::ifi) {
    m_thetai = -1;
    m_thetaj = 1;
  }
  if ((m_momenta.size() != m_oldmomenta.size()) || m_newmomenta.size() != 2 || m_bornmomenta.size() != 2) {
    THROW(fatal_error, "Incorrect dipole size in YFS");
  }
  m_thetaij = m_thetai*m_thetaj;
  m_theta.push_back(m_thetai);
  m_theta.push_back(m_thetaj);
  m_Q.push_back(m_Qi);
  m_Q.push_back(m_Qj);
  CalculateGamma();
}


Dipole::~Dipole() {
  Clean();
}



void Dipole::PrintInfo() {
   std::cout << " Dipole Type is "<<m_type
      << "\n Dipole components are "
      << m_names[0] << " " << m_names[1] << std::endl;
  for (int i = 0; i < 2; ++i)
  {
    std::cout << "Mass of " << m_names[i] << " = " << m_masses[i] << std::endl
        << "Charge of " << m_names[i] << " = " << m_charges[i] << std::endl
        << "Momentum of " << m_names[i] << " = " << m_momenta[i] << std::endl;
  }
  std::cout << "Invarinat mass " << " = " << (m_momenta[0]+m_momenta[1]).Mass() << std::endl;
  std::cout << "Number of Photons " << " = " << (m_dipolePhotons).size() << std::endl
            << "with four momentum :"<<std::endl;
  for(const auto &k: m_dipolePhotons){
    std::cout<<k<<std::endl;
  }
  if(m_type==dipoletype::final){
    std::string isres = (m_resonance)?"Yes":"No";
    std::cout << "Is Resonance: "<< isres << std::endl;
  }
}


void Dipole::Boost() {
  if (Type() == dipoletype::initial) {
    m_dipolePhotonsEEX=m_dipolePhotons;
    m_eikmomentum = m_bornmomenta;
    if (m_dipolePhotons.size() == 0) {
      DEBUG_FUNC("No ISR Photons, skipping boost");
      for (int i = 0; i < 2; ++i) m_newmomenta[i]=m_bornmomenta[i];
      return;
    }
    Vec4D Q;
    Q = m_bornmomenta[0] + m_bornmomenta[1] - m_photonSum;
    // if(Q.Abs2() > )
    double sp = Q * Q;
    double zz = sqrt(sp) / 2.;
    double z = zz * sqrt((sp - sqr(m_masses[0] - m_masses[1])) * (sp - sqr(m_masses[0] + m_masses[1]))) / sp;
    double m1 = m_masses[0];
    double m2 = m_masses[1];
    // m_newmomenta[0] = {zz, 0, 0, z};
    // m_newmomenta[1] = {zz, 0, 0, -z};
    double signz = m_bornmomenta[0][3]>0?1:-1;
    double lamCM = 0.5*sqrt(Lambda(Q.Abs2(),m1*m1,m2*m2)/Q.Abs2());
    double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
    double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
    m_newmomenta[0] = {E1, 0, 0, lamCM};
    m_newmomenta[1] = {E2, 0, 0, -lamCM};
    m_ranPhi = ran->Get()*2.*M_PI;
    // sqr(1.+2.*t/s)
    double s = (m_newmomenta[0]+m_newmomenta[1]).Abs2();
    double t = (m_newmomenta[0]-m_newmomenta[0]).Abs2();
    m_ranTheta = acos(1.+2.*t/s);
    ATOOLS::Poincare poin(Q);
    Poincare pRot(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
    for (int i = 0; i < 2; ++i) {
      pRot.RotateBack(m_newmomenta[i]);
      poin.BoostBack(m_newmomenta[i]);
    }
    m_sp = (m_newmomenta[0]+m_newmomenta[1]).Abs2();
  }
  else if (Type() == dipoletype::final) {
    if (m_dipolePhotons.size() == 0) {
      // Legs keep lab momenta (m_newmomenta) but the ghosts were made in the pair rest
      // frame and FSR::YFS_FORM reads them as lab vectors, so the lab and rest Piatek
      // terms would stop cancelling for a boosted pair.
      if (m_ghost.size() != 0) {
        p_Pboost.emplace(m_momenta[0] + m_momenta[1]);
        for (auto &r : m_ghost) Boost(r);
        p_Pboost.reset();
      }
      return;
    }
    if (m_dipolePhotons.size() != m_Nphotons){
      msg_Error()<<"Wrong Photon multiplicity in Boost \n"
                 <<"Photon vector size: "<<m_dipolePhotons.size()<<std::endl
                 <<"Photons Generated: "<<m_Nphotons<<std::endl;
    }
    if(!IsResonance()){
      msg_Error()<<"Trying to boost a non-resonant dipole"<<std::endl;
    }
    // Check that the final state fermions
    // are in their own restframe;
    Vec4D Q = m_momenta[0]+m_momenta[1];
    if(!IsEqual(0,Q.PSpat())){
      msg_Error()<<"Dipole is in the wrong frame\n";
    }
    if(m_ghost.size()!=0){
      Q = m_ghost[0]+m_ghost[1];
      if(!IsEqual(0,Q.PSpat())){
        msg_Error()<<"Dipole ghost is in the wrong frame";
      }
    }
    m_ranTheta = acos(1.-2.*ran->Get());
    m_ranPhi = ran->Get()*2.*M_PI;
    Vec4D qqk = m_momenta[0] + m_momenta[1] + m_photonSum;
    p_Pboost.emplace(qqk);
    m_eikmomentum = m_bornmomenta;
    for (size_t i = 0; i < 2; ++i)
    {
      Boost(m_momenta[i]);
      m_newmomenta[i]=m_momenta[i];
      if(m_ghost.size()!=0){
        Boost(m_ghost[i]);
      }
    }
    m_photonSum*=0.;
    // m_dipolePhotonsEEX.clear();
    for (auto &k : m_dipolePhotons) {
      m_dipolePhotonsEEX.push_back(k);
      Boost(k);
      m_photonSum+=k;
    }
    p_Pboost.reset();
  }
}

void Dipole::BoostNLO(ATOOLS::Vec4D &p) {
  p_Pboost->Boost(p);
  p_rotate.RotateBack(p);
  // RandomRotate(p);
  p_boost.BoostBack(p);
}

void Dipole::Boost(ATOOLS::Vec4D &p) {
  p_Pboost->Boost(p);
  p_rotate.RotateBack(p);
  p_boost.BoostBack(p);
}

bool Dipole::BoostNLO() {
  if (Type() == dipoletype::initial) {
    m_dipolePhotonsEEX=m_dipolePhotons;
    m_eikmomentum = m_bornmomenta;
    if (m_dipolePhotons.size() == 0) {
      DEBUG_FUNC("No ISR Photons, skipping boost");
      for (int i = 0; i < 2; ++i) m_newmomenta[i]=m_bornmomenta[i];
      return true;
    }
    Vec4D Q;
    Q = m_bornmomenta[0] + m_bornmomenta[1] - m_photonSum;
    // if(Q.Abs2() > )
    double sp = Q * Q;
    double zz = sqrt(sp) / 2.;
    double z = zz * sqrt((sp - sqr(m_masses[0] - m_masses[1])) * (sp - sqr(m_masses[0] + m_masses[1]))) / sp;
    double m1 = m_masses[0];
    double m2 = m_masses[1];
    // m_newmomenta[0] = {zz, 0, 0, z};
    // m_newmomenta[1] = {zz, 0, 0, -z};
    double signz = m_bornmomenta[0][3]>0?1:-1;
    double lamCM = 0.5*sqrt(Lambda(Q.Abs2(),m1*m1,m2*m2)/Q.Abs2());
    double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
    double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
    m_newmomenta[0] = {E1, 0, 0, signz*lamCM};
    m_newmomenta[1] = {E2, 0, 0, -signz*lamCM};
    m_ranPhi = ran->Get()*2.*M_PI;
    // sqr(1.+2.*t/s)
    double s = (m_newmomenta[0]+m_newmomenta[1]).Abs2();
    double t = (m_newmomenta[0]-m_newmomenta[0]).Abs2();
    m_ranTheta = acos(1.+2.*t/s);
    ATOOLS::Poincare poin(Q);
    Poincare pRot(m_bornmomenta[0], Vec4D(0., 0., 0., signz*1.));
    for (int i = 0; i < 2; ++i) {
      pRot.Rotate(m_newmomenta[i]);
      poin.BoostBack(m_newmomenta[i]);
    }
    m_sp = (m_newmomenta[0]+m_newmomenta[1]).Abs2();
    m_eikmomentum=m_newmomenta;
  }
  else if (Type() == dipoletype::final) {
    if (m_dipolePhotons.size() == 0) return true;
    m_photonSum*=0;
    for(auto &k: m_dipolePhotons) m_photonSum+=k;
    // if (m_dipolePhotons.size() != 1){
    //   msg_Error()<<"Wrong Photon multiplicity in BoostNLO \n"
    //              <<"Photon vector size: "<<m_dipolePhotons.size()<<std::endl
    //              <<"Photons Generated: "<<m_Nphotons<<std::endl;
    // }
    // if(!IsEqual(m_dipolePhotons[0],m_photonSum)){
    //   msg_Error()<<"Wrong photon momentum in "<<METHOD<<std::endl;
    // }
    // Check that the final state fermions
    // are in their own restframe;
    // m_ranTheta = acos(0.99999*(1.-2.*ran->Get()));
    m_ranPhi = ran->Get()*2.*M_PI;
    double s = (m_momenta[0]+m_momenta[1]).Abs2();
    double t = (m_momenta[0]-m_momenta[0]).Abs2();
    Vec4D Q = m_bornmomenta[0]+m_bornmomenta[1];
    m_ranTheta = acos((m_bornmomenta[0]+m_bornmomenta[1]).CosTheta());
    Poincare boost(m_bornmomenta[0]+m_bornmomenta[1]);
    // Poincare boost(m_newmomenta[0]+m_newmomenta[1]);
    // boost.Boost(m_photonSum);
    double x = 1./(1-m_photonSum.E());
    double y = 1./(1. + m_photonSum.E()/m_photonscale + 0.25*m_photonSum*m_photonSum/m_photonscale/m_photonscale);
    double sprim =(Q).Abs2()*y;
    Vec4D preboostk = m_photonSum;
    // if(IsBad(sprim)) return  false;
    // double m1 = m_momenta[0].Mass();
    // double m2 = m_momenta[1].Mass();
    // Vec4D rref = Q-m_photonSum;
    MakePair(sqrt(sprim), m_momenta[0], m_momenta[1]);
    // PHASIC::CE.Isotropic2Momenta(rref, m1*m1, m2*m2,m_momenta[0], m_momenta[1],ran->Get(), ran->Get());
    Vec4D qqk = m_momenta[0] + m_momenta[1] + m_photonSum;
    p_Pboost.emplace(qqk);
    Vec4D ref = m_bornmomenta[0];
    boost.Boost(ref);
    Poincare rot(ref, Vec4D(0,0,0,1));
    SetBoost(boost);
    SetRotate(rot);
    for (size_t i = 0; i < 2; ++i)
    {
      BoostNLO(m_momenta[i]);
      // p_boost.Boost(m_momenta[i]);
      // p_rotate.Rotate(m_momenta[i]);
      m_newmomenta[i]=m_momenta[i];
    }
    // m_eikmomentum = m_momenta;
    m_photonSum*=0.;
    for (auto &k : m_dipolePhotons) {
      BoostNLO(k);
      // p_Pboost->Boost(k);
      // // p_rotate.Rotate(k);
      // p_boost.BoostBack(k);
      m_photonSum+=k;
    }
    p_Pboost.reset();
    for (int i = 0; i < 2; ++i)
    {
      for(int j = 0; j < 4; ++j){
        double k = m_newmomenta[i][j];
        if(IsBad(k)){
         msg_Error()<<"NLO Boost Failed"<<std::endl;
         return false; 
        }
      }
    }
    // if(IsBad(m_newmomenta[0]) || IsBad(m_newmomenta[1]) ){
    //   msg_Error()<<"NLO Boost Failed"<<std::endl;
    //   return false;
    // }
    return true;
  }
  return true;
}


void Dipole::MakePair(double cms, Vec4D &p1, Vec4D &p2) {
  double E = cms / 2.;
  double s = sqr(cms);
  Vec4D P = {cms,0,0,0};
  double mass1 = p1.Mass();
  double mass2 = p2.Mass();
  double beta2 = (s - sqr(mass1 - mass2)) * (s - sqr(mass1 + mass2)) / (s * s);
  double beta =  sqrt(beta2);
  double eta1 = (s + sqr(mass1) - sqr(mass2)) / s;
  double eta2 = (s - sqr(mass1) + sqr(mass2)) / s;
  // p1 = {E * eta1, 0, 0, beta * E};
  // p2 = {E * eta2, 0, 0, -beta * E};
  double lamCM = 0.5*sqrt(Lambda(s,mass1*mass1,mass2*mass2)/s);
  double E1 = lamCM*sqrt(1+mass1*mass1/sqr(lamCM));
  double E2 = lamCM*sqrt(1+mass2*mass2/sqr(lamCM));
  p1 = {E1, 0, 0, lamCM};
  p2 =  {E2, 0, 0, -lamCM};
  Poincare boost(p1+p2);
  boost.Boost(m_photonSum);
  // p2 = P-p1;
}


void Dipole::RandomRotate(Vec4D &p){
    double x = p[1], y = p[2], z = p[3];

  // --- First: Rotate around X-axis (θ = m_ranTheta)
  // Affects Y and Z
  double cx = cos(m_ranTheta), sx = sin(m_ranTheta);
  double y1 = cx * y - sx * z;
  double z1 = sx * y + cx * z;

  // --- Then: Rotate around Z-axis (ϕ = m_ranPhi)
  // Affects X and new Y
  double cz = cos(m_ranPhi), sz = sin(m_ranPhi);
  double x1 = cz * x - sz * y1;
  double y2 = sz * x + cz * y1;

  // Set rotated components back into the vector
  p[1] = x1;
  p[2] = y2;
  p[3] = z1; // from the X rotation
}

void Dipole::BoostToQFM(bool boostback) {
  m_QFrame = m_bornmomenta[0] + m_bornmomenta[1];
  ATOOLS::Poincare poin(m_QFrame);
  for (auto &p : m_momenta) {
    if(boostback) poin.BoostBack(p);
    else poin.Boost(p);
  }
  // Recalcuate betas in this frame
  CalculateGamma();
}


// 1 - beta of a leg as m^2/(E(E+P)): no cancellation and never <= 0.
// 1 - |p|/E is ~1e-11 for an electron at LEP2 energies, so rounding can give
// |p| >= E, a NaN form factor and an endless Poisson loop in FSR::NPhotons.
static double OneMinusBeta(const Vec4D &p, double m)
{
  const double E(p.E()), m2(m*m);
  const double P(sqrt(Max(0., E*E - m2)));
  return m2/(E*(E + P));
}

/*
  gamma and gamma' of a FINAL-state dipole are Lorentz invariants of the
  pair, so the velocities have to be taken in the pair's own rest frame.
*/
void Dipole::CalculateGamma(){
  Vec4D q1(m_bornmomenta[0]), q2(m_bornmomenta[1]);
  if (m_type == dipoletype::final) {
    const Vec4D Q(q1 + q2);
    if (Q.Abs2() > 0. && Q[0] > 0. && !IsZero(Q.PSpat()/Q[0], 1e-14)) {
      Poincare rest(Q);
      rest.Boost(q1);
      rest.Boost(q2);
    }
  }
  const double omb1(OneMinusBeta(q1, m_masses[0]));
  const double omb2(OneMinusBeta(q2, m_masses[1]));
  m_b1 = 1. - omb1;
  m_b2 = 1. - omb2;
  double logarg = (1+m_b1)*(1+m_b2);
  logarg /= omb1*omb2;
  m_gamma  = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg)-2);
  m_gammap = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg));

  m_gamma  *= m_alpi*std::abs(ChargeNorm());
  m_gammap *= m_alpi*std::abs(ChargeNorm());
  if(Type()==dipoletype::final)   delf = 0.5*m_gamma;
  if(Type()==dipoletype::initial) deli = 0.5*m_gamma;
  if(RealOnly()) delf=deli=0;
}

void Dipole::AddPhotonsToDipole(ATOOLS::Vec4D_Vector &Photons) {
  m_photonSum *= 0;
  if (m_dipolePhotons.size() != 0) {
    msg_Debugging() << "Warning: Dipole still contains Photons, deleting old and adding new\n ";
    m_dipolePhotons.clear();
  }
  if (Photons.size() == 0) {
    DEBUG_FUNC("No Photons for this dipole" << this);
    return;
  }
  else {
    for (auto &k : Photons) AddPhotonToDipole(k);
  }
  DEBUG_FUNC("Photons added to this dipole " << this << "\n " << m_dipolePhotons);
}

void Dipole::AddPhotonToDipole(ATOOLS::Vec4D &k){
  m_dipolePhotons.push_back(k);
  m_photonSum +=k;
}

ATOOLS::Vec4D Dipole::Sum() {
  ATOOLS::Vec4D sum;
  for (auto m : m_bornmomenta) sum += m;
  return sum;
}

double Dipole::Mass() {
  return (m_flavs[0].Mass() + m_flavs[1].Mass()) / 2.;
}

void Dipole::AddToGhosts(ATOOLS::Vec4D &p) {
  if (m_ghost.size() > 2) {
    msg_Error() << "Too many four momentum in FSR for boosting" << std::endl;
  }
  m_ghost.push_back(p);
}

// EEX, Beta1/2/3, VirtualEEX, Hard and xi live in Dipole_EEX.C.

void Dipole::Clean(){
  m_masses.clear();
  m_charges.clear();
  m_names.clear();
  m_flavs.clear();
  m_momenta.clear();
  m_oldmomenta.clear();
  m_newmomenta.clear();
  m_bornmomenta.clear();
  m_ghost.clear();
  m_dipolePhotons.clear();
  m_photonSum*=0;
}

bool Dipole::IsDecayAllowed(){
  if(m_flavs[0].IsNeutrino() || m_flavs[1].IsNeutrino()){
    int diff = fabs(m_flavs[0].Kfcode() -m_flavs[1].Kfcode());
    if(diff==1) return true;
    else return false;
    // if(m_flavs[1])
  }
  else{
    if(m_flavs[0] == m_flavs[1].Bar() ) return true;
    else return false;
  }
}


double Dipole::Eikonal(const Vec4D &k,const Vec4D &p1,const Vec4D &p2) {
  // No extra sign for like-charge pairs: the dipole decomposition of the YFS
  // radiation function gives every unordered pair the coefficient
  // Q_iQ_j theta_i theta_j, and the mass terms resum by charge conservation.
  // Same convention as the Eikonal(k) overload. See YFS/Tools/IFI_Budget.C.
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double Dipole::EikonalMassless(const Vec4D &k,const Vec4D &p1, const Vec4D &p2) {
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (-2.*p1*p2 / ((p1 * k)*(p2 * k)));
}


double Dipole::Eikonal(const Vec4D &k) {
  Vec4D p1 = m_eikmomentum[0];
  Vec4D p2 = m_eikmomentum[1];
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}


double Dipole::EikonalInterferance(const Vec4D &k) {
  Vec4D p1 = m_eikmomentum[0];
  Vec4D p2 = m_eikmomentum[1];
  return -m_QiQj*m_thetaij*m_alp / (2 * M_PI * M_PI) * (p1*p2 / (p1 * k)/(p2 * k));
}

void Dipole::SetFlavLab(int i, int j){
  m_leftfl = i;
  m_rightfl = j;
}

std::ostream& YFS::operator<<(std::ostream &out, const Dipole &Dip) {
  out << " Dipole Type is "<<Dip.m_type
      << "\n Dipole components are "
      << Dip.m_names[0] << " " << Dip.m_names[1] << std::endl;
  for (int i = 0; i < 2; ++i)
  {
    out << "Mass of " << Dip.m_names[i] << " = " << Dip.m_masses[i] << std::endl
        << "Charge of " << Dip.m_names[i] << " = " << Dip.m_charges[i] << std::endl
        << "Momentum of " << Dip.m_names[i] << " = " << Dip.m_momenta[i] << std::endl
        << "Born Momentum of " << Dip.m_names[i] << " = " << Dip.m_bornmomenta[i] << std::endl;
  }
  out << "Invarinat mass " << " = " << (Dip.m_momenta[0]+Dip.m_momenta[1]).Mass() << std::endl
      <<"Sum of Photons = "<< Dip.m_photonSum << std::endl
      << "Q+sum_i K_i = "<< Dip.m_photonSum+Dip.m_momenta[0]+Dip.m_momenta[1]<<std::endl
      << "Born Qi+Qj = "<< Dip.m_bornmomenta[0]+Dip.m_bornmomenta[0]<<std::endl
      << "Left ID " << Dip.Left() << std::endl
      << "Right ID " <<Dip.Right() << std::endl
      << "Left Fl " << Dip.GetFlav(0) << std::endl
      << "Right Fl " << Dip.GetFlav(1) << std::endl
      << "Mass of photon-fermion system = "
      << (Dip.m_photonSum+Dip.m_newmomenta[0]+Dip.m_newmomenta[1]).Mass()<<std::endl;
  for (int i = 0; i < Dip.m_dipolePhotons.size(); ++i)
  {
    out<<" Photon["<<i<<"] = "<<Dip.m_dipolePhotons[i]<<std::endl;
  }
  if(Dip.m_type==dipoletype::final){
    std::string isres = (Dip.m_resonance)?"Yes":"No";
    out << "Is Resonance: "<< isres << std::endl;
  }
  return out;
}

std::ostream &YFS::operator<<(std::ostream &ostr,const dipoletype::code &it)
{
  if      (it==dipoletype::initial)  return ostr<<"Inital";
  else if (it==dipoletype::final)     return ostr<<"Final";
  else if (it==dipoletype::ifi)     return ostr<<"Initial-Final";
  return ostr<<"UNKNOWN";
}



EmissionResult Dipole::GenerateEmissions(ISR *isr, FSR *fsr,
                                         double born, double v,
                                         Vec4D_Vector &me_acc) {
  EmissionResult res;

  if (m_type == dipoletype::ifi) return res;   // carries interference, does not radiate

  if (m_type == dipoletype::initial) {
    if (!isr) THROW(fatal_error, "No ISR generator for an initial-initial dipole");
    isr->NPhotons();
    isr->GeneratePhotonMomentum();
    isr->Weight();
    SetBorn(born);
    const Vec4D_Vector k(isr->GetPhotons());
    for (const Vec4D &g : k) res.photons.push_back(Photon(g, this));
    res.me_photons = res.photons;   // ISR photons are not hidden
    res.weight = isr->GetWeight();
    res.photon_sum = isr->GetPhotonSum();
    Vec4D_Vector kk(k);
    AddPhotonsToDipole(kk);
    Boost();
    return res;
  }

  // dipoletype::final
  if (!fsr) THROW(fatal_error, "No FSR generator for a final-final dipole");
  fsr->Reset();
  BoostToQFM(0);
  SetBorn(born);
  fsr->SetV(v);
  if (!fsr->Initialize(*this)) { res.fail = EmissionResult::Failure::initialize; return res; }
  if (!fsr->MakeFSR())         { res.fail = EmissionResult::Failure::makefsr;    return res; }

  res.photon_sum = fsr->GetPhotonSum();

  Vec4D_Vector k(fsr->GetPhotons());
  if (!fsr->F()) { res.fail = EmissionResult::Failure::masswgt; res.weight = 0.; return res; }
  for (const Vec4D &g : k) {
    me_acc.push_back(g);
    res.me_photons.push_back(Photon(g, this));
  }
  AddPhotonsToDipole(k);
  SetMEPhotons(me_acc);
  Boost();
  if (!fsr->YFS_FORM()) { res.fail = EmissionResult::Failure::formfactor; return res; }
  fsr->HidePhotons();

  Vec4D_Vector hidden(fsr->GetPhotons());
  for (const Vec4D &g : hidden) res.photons.push_back(Photon(g, this));
  AddPhotonsToDipole(hidden);
  fsr->Weight();
  res.weight = fsr->GetWeight();
  return res;
}
