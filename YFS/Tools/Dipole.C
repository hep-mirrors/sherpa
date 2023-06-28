#include "YFS/Tools/Dipole.H"

#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Model_Base.H"
#include "MODEL/Main/Running_AlphaQED.H"



using namespace ATOOLS;
using namespace MODEL;
using namespace YFS;

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
  // if(use_model_alpha) m_alp = s_model->ScalarConstant("alpha_QED");
  m_alp  = (*aqed)(0); 
  m_alpi = m_alp/M_PI;
  if (use_model_alpha) m_rescale_alpha = 1;
  else m_rescale_alpha = (*aqed)(0) / s_model->ScalarConstant("alpha_QED");
  m_QiQj = 1;
  m_sp = (mom[0]+mom[1]).Abs2();
  for (auto &v : fl)
  {
    m_masses.push_back(v.Mass());
    m_charges.push_back(v.Charge());
    m_QiQj *= v.Charge();
    m_names.push_back(v.IDName());
    m_flavs.push_back(v);
  }
  for (auto &v : mom) {
    m_momenta.push_back(v);
    m_oldmomenta.push_back(v);
    m_newmomenta.push_back(v);
    // m_bornmomenta.push_back(v);
    m_beams.push_back(v);
  }
  for (auto &v : born) m_bornmomenta.push_back(v);
  if (ty == dipoletype::code::initial) {
    m_thetai = 1;
    m_thetaj = -1;
    m_thetaij = 1;
    for (int i = 0; i < 2; ++i) m_beams.push_back(m_bornmomenta[i]);
  }
  else if (ty == dipoletype::code::final) m_thetai = m_thetaj = m_thetaij = 1;
  else if (ty == dipoletype::code::ifi) {
    m_thetai = -1;
    m_thetaj = 1;
    m_thetaij = -1;
  }
  if ((m_momenta.size() != m_oldmomenta.size()) || m_newmomenta.size() != 2 || m_bornmomenta.size() != 2) {
    THROW(fatal_error, "Incorrect dipole size in YFS");
  }
  if (ty == dipoletype::code::final) {
    m_ghost.clear();
    // p_boost  = new Poincare(m_bornmomenta[0] + m_bornmomenta[1]);
    // p_rotate = new Poincare(m_bornmomenta[0], Vec4D(0., 0., 0., 1.));
  }
  CalculateBeta();
}


Dipole::~Dipole() {
  Clean();
}



void Dipole::PrintInfo() {
  std::cout << "Dipole components are (";

  for (size_t i = 0; i < m_masses.size(); ++i)
  {
    std::cout << "Mass of " << m_names[i] << " = " << m_masses[i] << std::endl
              << "Charge of " << m_names[i] << " = " << m_charges[i] << std::endl
              << "Momentum of " << m_names[i] << " = " << m_momenta[i] << std::endl;
  }

  // std::cout<<out;
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
    m_newmomenta[0] = {zz, 0, 0, z};
    m_newmomenta[1] = {zz, 0, 0, -z};
    ATOOLS::Poincare poin(Q);
    Poincare pRot(m_newmomenta[0], Vec4D(0., 0., 0., 1.));
    for (int i = 0; i < 2; ++i) {
      pRot.Rotate(m_newmomenta[i]);
      poin.BoostBack(m_newmomenta[i]);
    }
  }
  else if (Type() == dipoletype::final) {
    // if (m_dipolePhotons.size() == 0) return;
    if (m_dipolePhotons.size() != m_Nphotons){
      msg_Error()<<"Wrong Photon multiplicity in Boost \n";
    }
    // Check that the final state fermions
    // are in their own restframe;
    Vec4D Q = m_momenta[0]+m_momenta[1];
    if(!IsEqual(0,Q.PSpat())){
      msg_Error()<<"Dipole is in the wrong frame\n";
    }
    Q = m_ghost[0]+m_ghost[1];
    if(!IsEqual(0,Q.PSpat())){
      msg_Error()<<"Dipole ghost is in the wrong frame";
    }


    Vec4D qqk = m_momenta[0] + m_momenta[1] + m_photonSum;
    p_Pboost = new Poincare(qqk);
    p_boost  = new Poincare(m_oldmomenta[0] + m_oldmomenta[1]);

    p_rotate = new Poincare(m_beams[0], Vec4D(0., 0., 0., 1.));
    // p_rotate = new Poincare(Vec4D(125., 0., 0., 125.), Vec4D(0., 0., 0., 1.));
    for (size_t i = 0; i < 2; ++i)
    {
      Boost(m_momenta[i]);
      m_newmomenta[i]=m_momenta[i];
      Boost(m_ghost[i]);
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

  m_b1 = (Vec3D(m_momenta[0]).Abs() / m_momenta[0].E());
  m_b2 = (Vec3D(m_momenta[1]).Abs() / m_momenta[1].E());
  double logarg = (1+m_b1)*(1+m_b2);
  logarg /= (1-m_b1)*(1-m_b2);
  
  m_gamma  = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg)-2);
  m_gammap = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg));

  m_gamma  *= -m_alpi*m_QiQj;
  m_gammap *= -m_alpi*m_QiQj;

}

void Dipole::AddPhotonsToDipole(ATOOLS::Vec4D_Vector &Photons) {
  m_photonSum *= 0;
  if (m_dipolePhotons.size() != 0) {
    // msg_Error() << "Warning: Dipole still contains Photons, deleting old and adding new\n ";
    m_dipolePhotons.clear();
  }
  if (Photons.size() == 0) {
    DEBUG_FUNC("No Photons for this dipole" << this);
    return;
  }
  else {
    for (auto &k : Photons) {
      m_dipolePhotons.push_back(k);
      m_photonSum += k;
    }
  }
  DEBUG_FUNC("Photons added to this dipole " << this << "\n " << m_dipolePhotons);
}

ATOOLS::Vec4D Dipole::Sum() {
  ATOOLS::Vec4D sum;
  for (auto m : m_momenta) sum += m;
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
  if (Type() == dipoletype::initial) {
    double p1p2 = m_bornmomenta[0]*m_bornmomenta[1];
    double a = k*m_bornmomenta[0]/p1p2;
    double b = k*m_bornmomenta[1]/p1p2;
    double V = 1+m_gamma/2.;
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-b));
  }
  else if (Type() == dipoletype::final) {
    double p1p2 = m_bornmomenta[0]*m_bornmomenta[1];
    double ap = k*m_bornmomenta[0]/p1p2;
    double bp = k*m_bornmomenta[1]/p1p2;
    double V = 1+m_gamma/2.;
    double a = ap/(1.+ap+bp);
    double b = bp/(1.+ap+bp);
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-b));
  }
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
}

bool Dipole::IsDecayAllowed(){
  if(m_flavs[0].IsNeutrino() || m_flavs[1].IsNeutrino()){
    int diff = abs(m_flavs[0].Kfcode() -m_flavs[1].Kfcode());
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
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}


double Dipole::Eikonal(Vec4D k) {
  Vec4D p1 = m_bornmomenta[0];
  Vec4D p2 = m_bornmomenta[1];
  return m_QiQj*m_thetaij*m_alp / (4 * M_PI * M_PI) * (p1 / (p1 * k) - p2 / (p2 * k)).Abs2();
}



std::ostream& YFS::operator<<(std::ostream &out, const Dipole &Dip) {
  out << " Dipole Type is "<<Dip.m_type
      << "\n Dipole components are "
      << Dip.m_names[0] << " " << Dip.m_names[1] << std::endl;
  for (int i = 0; i < 2; ++i)
  {
    out << "Mass of " << Dip.m_names[i] << " = " << Dip.m_masses[i] << std::endl
        << "Charge of " << Dip.m_names[i] << " = " << Dip.m_charges[i] << std::endl;
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

