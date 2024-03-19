#include "YFS/Tools/Dipole.H"

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "ATOOLS/Org/Run_Parameter.H" 



using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;

double SqLam(double x,double y,double z)
{
  return abs(x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z);
  // double arg(sqr(s-s1-s2)-4.*s1*s2);
  // if (arg>0.) return sqrt(arg)/s;
  // return 0.;
}


Dipole::Dipole(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born, dipoletype::code ty):
  m_type(ty)
{
  if ((mom.size() != fl.size()) || fl.size() != 2 || mom.size() != 2 || born.size()!=2) {
    msg_Out()<<"Dipole type is  =  "<<ty<<std::endl
             <<" mom.size() =  "<<mom.size()<<std::endl
             <<" fl.size() =  "<<fl.size()<<std::endl
             <<" born.size() =  "<<born.size()<<std::endl;
    THROW(fatal_error, "Incorrect dipole size in YFS for dipoletype");
  }
  Clean();
  // todo get alpha from YFS_BASE
  Scoped_Settings s{ Settings::GetMainSettings()["YFS"] };
  bool use_model_alpha = s["USE_MODEL_ALPHA"].Get<bool>();
  if(use_model_alpha) m_alp = s_model->ScalarConstant("alpha_QED");
  else m_alp  = 1./s["1/ALPHAQED(0)"].Get<double>(); 
  m_alpi = m_alp/M_PI;
  if (use_model_alpha) m_rescale_alpha = 1;
  else m_rescale_alpha = (*aqed)(0) / s_model->ScalarConstant("alpha_QED");
  // m_QiQj = 1;
  m_sp = (mom[0]+mom[1]).Abs2();
  m_Qi = fl[0].Charge();
  m_Qj = fl[1].Charge();
  m_QiQj = m_Qi*m_Qj;
  if(IsEqual(fl[0],fl[1])) m_sameflav = 1;
  else m_sameflav = 0;
  for (auto &v : fl)
  {
    m_masses.push_back(v.Mass());
    m_charges.push_back(v.Charge());
    m_names.push_back(v.IDName());
    m_flavs.push_back(v);
  }
  for (auto &v : mom) {
    m_momenta.push_back(v);
    m_oldmomenta.push_back(v);
    m_newmomenta.push_back(v);
    m_eikmomentum.push_back(v);
    m_beams.push_back(v);
    m_ghost.push_back(v);
  }
  for (auto &v : born) m_bornmomenta.push_back(v);
  if (ty == dipoletype::code::initial) {
    if(fl[0].IsAnti()) m_thetai = -1;
    else m_thetai = 1;
    if(fl[1].IsAnti()) m_thetaj = 1;
      else m_thetaj = -1;
    if(IsEqual(m_Qi,m_Qj)){
      m_thetai = m_thetaj = -1;
    }
    // m_thetai = m_thetaj = -1;
    for (int i = 0; i < 2; ++i) m_beams.push_back(m_bornmomenta[i]);
  }
  else if (ty == dipoletype::code::final) {
    if(fl[0].IsAnti()) m_thetai = 1;
    else m_thetai = -1;
    if(fl[1].IsAnti()) m_thetaj = -1;
      else m_thetaj = 1;
    if(IsEqual(m_Qi,m_Qj)){
      m_thetai = m_thetaj = 1;
    }
    // m_thetai = m_thetaj = 1;

  }
  else if (ty == dipoletype::code::ifi) {
    m_thetai = -1;
    m_thetaj = 1;
  }
  if ((m_momenta.size() != m_oldmomenta.size()) || m_newmomenta.size() != 2 || m_bornmomenta.size() != 2) {
    THROW(fatal_error, "Incorrect dipole size in YFS");
  }
  if (ty == dipoletype::code::final) {
    m_ghost.clear();
    // p_boost  = new Poincare(m_bornmomenta[0] + m_bornmomenta[1]);
    // p_rotate = new Poincare(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
  }
  m_thetaij = m_thetai*m_thetaj;
  m_theta.push_back(m_thetai);
  m_theta.push_back(m_thetaj);
  m_Q.push_back(m_Qi);
  m_Q.push_back(m_Qj);
  CalculateBeta();
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
  if(m_type==dipoletype::final){
    std::string isres = (m_resonance)?"Yes":"No";
    std::cout << "Is Resonance: "<< isres << std::endl;
  }
}


void Dipole::Boost() {
  if (Type() == dipoletype::initial) {
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
    double lamCM = 0.5*sqrt(SqLam(Q.Abs2(),m1*m1,m2*m2)/Q.Abs2());
    double E1 = lamCM*sqrt(1+m1*m1/sqr(lamCM));
    double E2 = lamCM*sqrt(1+m2*m2/sqr(lamCM));
    m_newmomenta[0] = {E1, 0, 0, lamCM};
    m_newmomenta[1] = {E2, 0, 0, -lamCM};
    ATOOLS::Poincare poin(Q);
    Poincare pRot(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
    for (int i = 0; i < 2; ++i) {
      pRot.Rotate(m_newmomenta[i]);
      poin.BoostBack(m_newmomenta[i]);
    }
  }
  else if (Type() == dipoletype::final) {
    if (m_dipolePhotons.size() == 0) return;
    if (m_dipolePhotons.size() != m_Nphotons){
      msg_Error()<<"Wrong Photon multiplicity in Boost \n";
    }
    // Check that the final state fermions
    // are in their own restframe;
    Vec4D Q = m_momenta[0]+m_momenta[1];
    // if(!IsEqual(0,Q.PSpat())){
    //   msg_Error()<<"Dipole is in the wrong frame\n";
    // }
    if(m_ghost.size()!=0){
      Q = m_ghost[0]+m_ghost[1];
      // if(!IsEqual(0,Q.PSpat())){
        // msg_Error()<<"Dipole ghost is in the wrong frame";
      // }
    }

    Vec4D qqk = m_momenta[0] + m_momenta[1] + m_photonSum;
    p_Pboost = new Poincare(qqk);
    p_boost  = new Poincare(m_beams[0] + m_beams[1]);

    p_rotate = new Poincare(m_beams[0], Vec4D(0., 0., 0., 1.));
    for (size_t i = 0; i < 2; ++i)
    {
      Boost(m_momenta[i]);
      m_newmomenta[i]=m_momenta[i];
      if(m_ghost.size()!=0){
        Boost(m_ghost[i]);
      }
    }
    m_photonSum*=0.;
    for (auto &k : m_dipolePhotons) {
      Boost(k);
      m_photonSum+=k;
    }
    if (p_rotate) delete p_rotate;
    if (p_Pboost) delete p_Pboost;
    if (p_boost) delete p_boost;
  }
}

void Dipole::Boost(ATOOLS::Vec4D &p) {
  p_Pboost->Boost(p);
  p_rotate->RotateBack(p);
  p_boost->BoostBack(p);
}

void Dipole::BoostLab(){
  Poincare p_boost(m_beams[0] + m_beams[1]);
  p_boost.BoostBack(m_newmomenta[0]);
  p_boost.BoostBack(m_newmomenta[1]);
  for(auto &k : m_dipolePhotons) p_boost.BoostBack(k);
  // if (p_boost) delete p_boost;
}

void Dipole::BoostToCMS(Vec4D_Vector &k, bool boostback){
  Vec4D CMSFrame=m_bornmomenta[0] + m_bornmomenta[1];
  ATOOLS::Poincare poin(m_QFrame);
  for (auto &p : k) {
    if(boostback) poin.BoostBack(p);
    else poin.Boost(p);
  }
}


void Dipole::BoostToQFM(bool boostback) {
  m_QFrame = m_bornmomenta[0] + m_bornmomenta[1];
  ATOOLS::Poincare poin(m_QFrame);
  for (auto &p : m_momenta) {
    if(boostback) poin.BoostBack(p);
    else poin.Boost(p);
  }
  // Recalcuate betas in this frame
  CalculateBeta();
}


void Dipole::CalculateBeta(){

  m_b1 = (Vec3D(m_bornmomenta[0]).Abs() / m_bornmomenta[0].E());
  m_b2 = (Vec3D(m_bornmomenta[1]).Abs() / m_bornmomenta[1].E());
  double logarg = (1+m_b1)*(1+m_b2);
  logarg /= (1-m_b1)*(1-m_b2);
  m_gamma  = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg)-2);
  m_gammap = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg));

  m_gamma  *= m_alpi*abs(ChargeNorm());
  m_gammap *= m_alpi*abs(ChargeNorm());
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
    for (auto &k : Photons) AddPhotonsToDipole(k);
  }
  DEBUG_FUNC("Photons added to this dipole " << this << "\n " << m_dipolePhotons);
}

void Dipole::AddPhotonsToDipole(ATOOLS::Vec4D &k){
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

double Dipole::EEX(const Vec4D &k){
  double p1p2 = m_bornmomenta[0]*m_bornmomenta[1];
  double a = k*m_bornmomenta[0]/p1p2;
  double b = k*m_bornmomenta[1]/p1p2;
  double ap = a/(1.+a+b);
  double bp = b/(1.+a+b);
    // double V = 1+m_gamma/2.;
  if (Type() == dipoletype::initial) {
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-b));
  }
  else if (Type() == dipoletype::final) {
    return 0.5*Eikonal(k)*(sqr(1-ap)+sqr(1-bp));
  }
  // else if (Type() == dipoletype::ifi) {
  //   return -2*Eikonal(k)*((sqr(1-ap)+sqr(1-b))+(sqr(1-a)+sqr(1-bp)));
  // }
  return 0;
}

void Dipole::Clean(){
  m_masses.clear();
  m_charges.clear();
  m_names.clear();
  m_flavs.clear();
  m_momenta.clear();
  m_oldmomenta.clear();
  m_newmomenta.clear();
  m_bornmomenta.clear();
  m_beams.clear();
  m_ghost.clear();
  m_dipolePhotons.clear();
  m_photonSum*=0;
  m_theta.clear();
  m_Q.clear();
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


double Dipole::Eikonal(Vec4D k, Vec4D p1, Vec4D p2) {
  // Vec4D eik;
  // if(Type()==dipoletype::initial){
  //   eik+=m_thetai*m_Qi*p1/(p1*k);
  //   eik+=m_thetaj*m_Qj*p2/(p2*k);
  // }
  // else if(Type()==dipoletype::final){
  //   eik+=m_thetai*m_Qi*p1/(p1*k);
  //   eik+=m_thetaj*m_Qj*p2/(p2*k);
  // }
  // return m_alp/(4 * M_PI * M_PI)*eik*eik;
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}

double Dipole::EikonalMassless(Vec4D k, Vec4D p1, Vec4D p2) {
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (-2.*p1*p2 / ((p1 * k)*(p2 * k)));
}


double Dipole::Eikonal(Vec4D k) {
  Vec4D p1 = m_bornmomenta[0];
  Vec4D p2 = m_bornmomenta[1];
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}


double Dipole::EikonalInterferance(const Vec4D &k) {
  Vec4D p1 = m_bornmomenta[0];
  Vec4D p2 = m_bornmomenta[1];
  return -m_QiQj*m_thetaij*m_alp / (2 * M_PI * M_PI) * (p1*p2 / (p1 * k)/(p2 * k));
}



std::ostream& YFS::operator<<(std::ostream &out, const Dipole &Dip) {
  out << " Dipole Type is "<<Dip.m_type
      << "\n Dipole components are "
      << Dip.m_names[0] << " " << Dip.m_names[1] << std::endl;
  for (int i = 0; i < 2; ++i)
  {
    out << "Mass of " << Dip.m_names[i] << " = " << Dip.m_masses[i] << std::endl
        << "Charge of " << Dip.m_names[i] << " = " << Dip.m_charges[i] << std::endl
        << "Momentum of " << Dip.m_names[i] << " = " << Dip.m_momenta[i] << std::endl;
  }
  out << "Invarinat mass " << " = " << (Dip.m_momenta[0]+Dip.m_momenta[1]).Mass() << std::endl
      <<"Sum of Photons = "<< Dip.m_photonSum << std::endl
      << "Q+sum_i K_i = "<< Dip.m_photonSum+Dip.m_momenta[0]+Dip.m_momenta[1]<<std::endl
      << "Mass of photon-fermion system = "
      << (Dip.m_photonSum+Dip.m_newmomenta[0]+Dip.m_newmomenta[1]).Mass()<<std::endl;
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

