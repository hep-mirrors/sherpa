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


// for IFI corrections
double delf = 0;
double deli = 0;
int order = 0;
double SqLam(double x,double y,double z)
{
  return abs(x*x+y*y+z*z-2.*x*y-2.*x*z-2.*y*z);
  // double arg(sqr(s-s1-s2)-4.*s1*s2);
  // if (arg>0.) return sqrt(arg)/s;
  // return 0.;
}


Dipole::Dipole(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, 
              ATOOLS::Vec4D_Vector const &born, dipoletype::code ty, const double alpha):
  m_type(ty), m_alp(alpha)
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
  m_alpi = m_alp/M_PI;
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
    m_beams.push_back(v);
    m_ghost.push_back(v);
  }
  for (auto &v : born) m_bornmomenta.push_back(v);
  m_eikmomentum = m_bornmomenta;
  if (ty == dipoletype::code::initial) {
    if(fl[0].IsAnti()) m_thetai = -1;
    else m_thetai = 1;
    if(fl[1].IsAnti()) m_thetaj = 1;
      else m_thetaj = -1;
    if(IsEqual(m_Qi,m_Qj)){
      m_thetai = m_thetaj = -1;
    }
    // m_thetai = m_thetaj = -1;
    // for (int i = 0; i < 2; ++i) m_beams.push_back(m_bornmomenta[i]);
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
    double lamCM = 0.5*sqrt(SqLam(Q.Abs2(),m1*m1,m2*m2)/Q.Abs2());
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
  }
  else if (Type() == dipoletype::final) {
    if (m_dipolePhotons.size() == 0) return;
    if (m_dipolePhotons.size() != m_Nphotons){
      msg_Error()<<"Wrong Photon multiplicity in Boost \n";
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
    m_ranPhi = ran->Get()*2.*M_PI;
    // sqr(1.+2.*t/s)
    // m_eikmomentum = m_momenta;
    double s = (m_bornmomenta[0]+m_bornmomenta[1]).Abs2();
    double t = (m_bornmomenta[0]-m_bornmomenta[1]).Abs2();
    // m_ranTheta = acos(1.+2.*t/s);
    m_ranTheta = acos(1.-2.*ran->Get());
    // m_ranTheta = m_beams[0].Theta();
    Vec4D qqk = m_momenta[0] + m_momenta[1] + m_photonSum;
    p_Pboost = new Poincare(qqk);
    p_boost  = new Poincare(m_bornmomenta[0] + m_bornmomenta[1]);

    // p_rotate = new Poincare(m_bornmomenta[0], Vec4D(0., 0.,  0., 1.));
    p_rotatey = new Poincare(m_bornmomenta[0], Vec4D(0., 0., 1., 0.));
    p_rotatex = new Poincare(m_bornmomenta[0], Vec4D(0., 1., 0., 0.));
    for (size_t i = 0; i < 2; ++i)
    {
      p_Pboost->Boost(m_momenta[i]);
      if(i==0) p_rotate = new Poincare(m_bornmomenta[0],Vec4D(0., 0.,  0., 1.));
      Boost(m_momenta[i]);
      m_newmomenta[i]=m_momenta[i];
      if(m_ghost.size()!=0){
        Boost(m_ghost[i]);
      }
    }
    m_eikmomentum = m_momenta;
    m_photonSum*=0.;
    // m_dipolePhotonsEEX.clear();
    for (auto &k : m_dipolePhotons) {
      // Boost(k);
      p_Pboost->Boost(k);
      Boost(k);
      m_dipolePhotonsEEX.push_back(k);
      m_photonSum+=k;
    }
    if (p_rotate) delete p_rotate;
    if (p_rotatex) delete p_rotatey;
    if (p_rotatey) delete p_rotatex;
    if (p_Pboost) delete p_Pboost;
    if (p_boost) delete p_boost;
  }
}

void Dipole::Boost(ATOOLS::Vec4D &p) {
  // p_Pboost->Boost(p);
  p_rotate->Rotate(p);
  p_boost->BoostBack(p);
}

void Dipole::RandomRotate(Vec4D &p){
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta)*t1[2] - sin(m_ranTheta)*t1[3];
  p[3] = sin(m_ranTheta)*t1[2] + cos(m_ranTheta)*t1[3];
  t1 = p;
  p[1] = cos(m_ranPhi)*t1[1]-sin(m_ranPhi)*t1[2];
  p[2] = sin(m_ranPhi)*t1[1]+cos(m_ranPhi)*t1[2];
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
  poin.Boost(m_eikmomentum[0]);
  poin.Boost(m_eikmomentum[1]);
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


void Dipole::CalculateGamma(){
  m_b1 = (Vec3D(m_eikmomentum[0]).Abs() / m_eikmomentum[0].E());
  m_b2 = (Vec3D(m_eikmomentum[1]).Abs() / m_eikmomentum[1].E());
  double logarg = (1+m_b1)*(1+m_b2);
  logarg /= (1-m_b1)*(1-m_b2);
  m_gamma  = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg)-2);
  m_gammap = (1.+m_b1*m_b2)/(m_b1+m_b2)*(log(logarg));

  m_gamma  *= m_alpi*abs(ChargeNorm());
  m_gammap *= m_alpi*abs(ChargeNorm());
  if(Type()==dipoletype::final)   delf = 0.5*m_gamma;
  if(Type()==dipoletype::initial) deli = 0.5*m_gamma;
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

double Dipole::EEX(const int betaorder){
  double real=0;
  if(m_dipolePhotonsEEX.size()==0) return real;
  CalculateGamma();
  m_betaorder = betaorder;
  if(betaorder >= 1 && Type()!=dipoletype::ifi) {
    for(auto &k: m_dipolePhotonsEEX){
     real += Beta1(k)/Eikonal(k);
    }
  }
  if(betaorder >= 2 ) {
    for (int j = 1; j < m_dipolePhotonsEEX.size(); j++) {
      for (int i = 0; i < j; i++) {
        // m_betaorder-=1;
        Vec4D k1 = m_dipolePhotonsEEX[j];
        Vec4D k2 = m_dipolePhotonsEEX[i];
        real += Beta2(k1,k2)/Eikonal(k1)/Eikonal(k2);
        m_betaorder = betaorder;
      }
    }
  }
  if(betaorder >= 3){
    for (int j = 1; j < m_dipolePhotonsEEX.size(); j++) {
      for (int i = 0; i < j; i++) {
        for (int k = 0; k < i; k++) {
          Vec4D k1 = m_dipolePhotonsEEX[j];
          Vec4D k2 = m_dipolePhotonsEEX[i];
          Vec4D k3 = m_dipolePhotonsEEX[k];
          double eik1 = Eikonal(k1);
          double eik2 = Eikonal(k2);
          double eik3 = Eikonal(k3);
          real += Beta3(k1,k2,k3)/eik1/eik2/eik3;
          m_betaorder = betaorder;
        }
      }
    }
  }
  if(IsNan(real)){
    msg_Error()<<"YFS EEX is NaN at order "<<betaorder<<std::endl;
  }
  return real;//+virt;
}

double Dipole::Beta1(const Vec4D &k){
  double b1=0;
  if(Type()==dipoletype::initial) {
  //   // beta11
    if(m_betaorder==2) b1 = Hard(k)*(1+delf)-Eikonal(k)*(1+deli)*(1+delf);
    // beta12
    else if(m_betaorder==3) b1 = Hard(k)-Eikonal(k)*(1+delf+0.5*delf*delf)*(1+deli+0.5*deli*deli);
    else b1 = (Hard(k)-Eikonal(k))*(1+delf);
  }
  else if(Type()==dipoletype::final) {
    // if(m_betaorder==2) b1 = Hard(k)*(1+0.5*deli-0.25*deli)*(1+delf)-Eikonal(k);
    if(m_betaorder==2) b1 = Hard(k)*(1+deli)-Eikonal(k)*(1.+deli)*(1+delf);
    else if(m_betaorder==3) b1 = Hard(k)*(1+deli+0.5*deli*deli)-Eikonal(k)*(1+delf+0.5*delf*delf)*(1+deli+0.5*deli*deli);
    else b1 = Hard(k)-Eikonal(k);
  }
  else{
    b1 = Hard(k)-Eikonal(k);
  }
  return b1;
}

double Dipole::Beta2(const Vec4D &k1, const Vec4D &k2){
  double eik1 = Eikonal(k1);
  double eik2 = Eikonal(k2);
  double delta=1, hard;
  if(m_betaorder==3) delta = (1+delf)*(1+deli);
  hard = Hard(k1,k2);
  m_betaorder-=1;// Reduce order to calculate beta1(n-1)
  hard+= -eik1*(Beta1(k2))
         -eik2*(Beta1(k1))
         -eik1*eik2;
  m_betaorder+=1;
  return hard*delta;
}

double Dipole::Beta3(const Vec4D &k1, const Vec4D &k2, const Vec4D &k3){
  double eik1 = Eikonal(k1);
  double eik2 = Eikonal(k2);
  double eik3 = Eikonal(k3);
  double del = 0;
  double b3 = 0;
  if(Type()!=dipoletype::initial) return 0;
  if(Type()==dipoletype::initial) del = deli;
  else del = delf;
  double hard = Hard(k1,k2,k3);
  m_betaorder=-1;// Reduce order to calculate beta1(n-1)
  hard += -eik1*Beta2(k3,k2)
          -eik2*Beta2(k3,k1)
          -eik3*Beta2(k1,k2)
          -eik1*eik2*eik3;
  m_betaorder=-2;
  hard += -eik2*eik3*Beta1(k1)
          -eik1*eik3*Beta1(k2)
          -eik1*eik2*Beta1(k3);
  return hard;
}

double Dipole::VirtualEEX(const int betaorder){
  double virt{0};
  // For ISR+FSR virtuals are taken in for ISRxFSR not ISR+FSR
  if(betaorder==1)  virt =  0.5*m_gamma;
  else if(betaorder==2) virt = 0.5*m_gamma + 0.125*m_gamma*m_gamma;
  else if(betaorder==3) virt = 0.5*m_gamma + 0.125*m_gamma*m_gamma + pow(m_gamma,3)/48;
  return virt;
}

double Dipole::Hard(const Vec4D &k, int i){
  double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];
  double a = k*m_eikmomentum[0]/p1p2;
  double b = k*m_eikmomentum[1]/p1p2;
  double ap = a/(1.+a+b);
  double bp = b/(1.+a+b);
  double delta = 0;
  if (Type() == dipoletype::initial) {
    double z = (1-a)*(1-b);
    if(m_betaorder>=2){
      delta += 0.5*m_gamma
              +m_alpi*(log(a)*log(1-b)+log(b)*log(1-a)
                      +DiLog(a) + DiLog(b)
                      -0.5*sqr(log(1-a))-0.5*sqr(log(1-b))
                      +1.5*log(1-a)+1.5*log(1-b)
                      +0.5*a*(1-a)/(1+sqr(1-a))
                      +0.5*b*(1-b)/(1+sqr(1-b)));
    }
    if(m_betaorder>=3){
      delta += 0.125*sqr(m_gamma)*(1-log(z))
             +sqr(m_gamma)/24 *sqr(log(z));
    } 
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-b))*(1+delta);
  }
  else if (Type() == dipoletype::final) {
    double z = (1-ap)*(1-bp);
    if(m_betaorder>=2){
      delta += 0.5*m_gamma+0.25*m_gamma*log(z);
    }
    return 0.5*Eikonal(k)*(sqr(1-ap)+sqr(1-bp))*(1+delta);
  }
  else if (Type() == dipoletype::ifi) {
    return 0.5*Eikonal(k)*(sqr(1-a)+sqr(1-bp));
  }
  return 0;
}

double Dipole::Hard(const Vec4D &k1, const Vec4D &k2){
  double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];
  
  double a1 = k1*m_eikmomentum[0]/p1p2;
  double a2 = k2*m_eikmomentum[0]/p1p2;
  
  double b1 = k1*m_eikmomentum[1]/p1p2;
  double b2 = k2*m_eikmomentum[1]/p1p2;
  
  double eta1 = a1/(1+a1+b1);
  double eta2 = a2/(1+a2+b2);

  double zeta1 = b1/(1+a1+b1);
  double zeta2 = b2/(1+a2+b2);

  double etap1 = eta1/(1+eta2);
  double etap2 = eta2/(1+eta1);;

  double zetap1 = zeta1/(1+zeta2);
  double zetap2 = zeta2/(1+zeta1);

  double ap1 = a1/(1.-a2);
  double bp1 = b1/(1.-b2);

  double ap2 = a2/(1.-a1);
  double bp2 = b2/(1.-b1);
  
  double v1 = a1+b1;
  double v2 = a2+b2;
  double hard,delta{1};
  if (Type() == dipoletype::initial) {
    if(v1 > v2){
      hard = xi(a1,ap2,bp2) + xi(b1,ap2,bp2);
    } 
    else{
      hard = xi(a2,ap1,bp1) + xi(b2,ap1,bp1);
    }
    // if(m_betaorder==3){
    //   delta = 1+delf;
    // }
    return Eikonal(k1)*Eikonal(k2)*hard*delta;
  }
  else if (Type() == dipoletype::final) {
    if(v1 > v2){
      hard = xi(eta1,etap2,zetap2) + xi(zeta1,etap2,zetap2);
    }
    else {
      hard = xi(eta2,etap1,zetap1) + xi(zeta2,etap1,zetap1);
    }
    // if(m_betaorder==3){
    //   delta = 1+deli;
    // }
    return Eikonal(k1)*Eikonal(k2)*hard*delta;
  }
  else if(Type()== dipoletype::initial){
    if(v1 > v2){
      hard = xi(a1,ap2,bp2) + xi(zeta1,etap2,zetap2);
    } 
    else{
      hard = xi(a2,ap1,bp1) + xi(zeta2,etap1,zetap1);
    }
  }
  return 0;
}

double Dipole::Hard(const Vec4D &k1, const Vec4D &k2, const Vec4D &k3){
  double p1p2 = m_eikmomentum[0]*m_eikmomentum[1];
  
  double a1 = k1*m_eikmomentum[0]/p1p2;
  double a2 = k2*m_eikmomentum[0]/p1p2;
  double a3 = k3*m_eikmomentum[0]/p1p2;
  
  double b1 = k1*m_eikmomentum[1]/p1p2;
  double b2 = k2*m_eikmomentum[1]/p1p2;
  double b3 = k3*m_eikmomentum[1]/p1p2;
  
  double eta1 = a1/(1+a1+b1);
  double eta2 = a2/(1+a2+b2);
  double eta3 = a3/(1+a3+b3);

  double zeta1 = b1/(1+a1+b1);
  double zeta2 = b2/(1+a2+b2);
  double zeta3 = b3/(1+a3+b3);

  double etap1 = eta1/(1+eta2);
  double zetap1 = zeta1/(1+zeta2);

  double etap2 = eta2/(1+eta1);;
  double zetap2 = zeta2/(1+zeta1);

  double etap3  = eta3/(1+eta1+eta3);
  double zetap3 = zeta3/(1+zeta1+zeta2);

  double ap1 = a1/(1.-a2);
  double bp1 = b1/(1.-b2);

  double ap2 = a2/(1.-a1);
  double bp2 = b2/(1.-b1);

  double ap3 = a3/(1-a1-a2);
  double bp3 = b3/(1-b1-b2);
  
  double v1 = a1+b1;
  double v2 = a2+b2;
  double hard;
  if (Type() == dipoletype::initial) {
    if(v1 > v2){
      hard = xi(a1,ap2,bp2,ap3,bp3) + xi(b1,ap2,bp2,ap3,bp3);
    } 
    else{
      hard = xi(a2,ap1,bp1,ap3,bp3) + xi(b2,ap1,bp1,ap3,bp3);
    }
    return Eikonal(k1)*Eikonal(k2)*Eikonal(k3)*hard;
  }
  else if (Type() == dipoletype::final) {
    // return 0; // not implemented as << 0
    if(v1 > v2){
      hard = xi(eta1,etap2,zetap2,etap3,zetap3) + xi(zeta1,etap2,zetap2,etap3,zetap3);
    }
    else{
      hard = xi(eta2,etap1,zetap1,etap3,zetap3) + xi(zeta2,etap1,zetap1,etap3,zetap3);
    }
    return Eikonal(k1)*Eikonal(k2)*Eikonal(k3)*hard;
  }
  return 0;
}


double Dipole::xi(const double &alp, const double &beta, const double &gamma){
  return 0.25*sqr(1.-alp)*(sqr(1.-beta)+sqr(1.-gamma));
}

double Dipole::xi(const double &alp, const double &a1, const double &b1, const double &a2, const double &b2){
  return 0.125*sqr(1.-alp)*(sqr(1.-a1)+sqr(1.-b1))*(sqr(1.-a2)+sqr(1.-b2));
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
  m_eikmomentum.clear();
  m_beams.clear();
  m_ghost.clear();
  m_dipolePhotons.clear();
  m_dipolePhotonsEEX.clear();
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


double Dipole::Eikonal(const Vec4D &k,const Vec4D &p1,const Vec4D &p2) {
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

