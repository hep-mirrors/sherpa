#include "PHOTONS++/Tools/Dress_Blob_Base.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "PHOTONS++/Main/Photons.H"
#include "PHOTONS++/Tools/Generate_One_Photon.H"
#include "PHOTONS++/Tools/Weight_Dipole.H"
#include "PHOTONS++/Tools/Weight_Jacobian.H"
#include "PHOTONS++/Tools/Weight_YFS.H"
#include "PHOTONS++/Tools/Weight_Higher_Order_Corrections.H"
#include "PHOTONS++/Tools/Generate_One_Photon.H"
#include "PHOTONS++/PhaseSpace/Avarage_Photon_Number.H"
#include "ATOOLS/Org/Run_Parameter.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Dress_Blob_Base::Dress_Blob_Base() :
  m_photonsadded(false), m_success(false), m_soft(false),
  p_newinitialstate(NULL), m_accu(PHOTONS::Photons::s_accu),
  m_genweight(1.), m_genmaxweight(1.), m_nbar(0.), m_n(0.),
  m_omegaMin(PHOTONS::Photons::s_ircutoff)
{
  if (PHOTONS::Photons::s_mode == 1) m_soft = true;
}

Dress_Blob_Base::~Dress_Blob_Base() {
  if (m_dtype == Dipole_Type::ff)  delete p_newinitialstate;
}

void Dress_Blob_Base::BuildNewParticleVectorVector() {
  if ((m_dtype == Dipole_Type::ff) && (p_newinitialstate != NULL)) delete p_newinitialstate;
  m_pvv_new.clear();
  Particle_Vector cinit, ninit;
  if (m_dtype == Dipole_Type::ff) {
    p_newinitialstate = new Particle(*m_neutralinparticles[0]);
    p_newinitialstate->SetMomentum(m_P+m_PN+m_K);
    ninit.push_back(p_newinitialstate);
  }
  else if (m_dtype == Dipole_Type::fi) {
    p_newinitialstate = m_newdipole[0];
    cinit.push_back(p_newinitialstate);
  }
  m_pvv_new.push_back(cinit);
  m_pvv_new.push_back(ninit);
  m_pvv_new.push_back(m_newdipole);
  m_pvv_new.push_back(m_newspectator);
  m_pvv_new.push_back(m_softphotons);

#ifdef PHOTONS_DEBUG
  for (unsigned int i=0; i<m_pvv_new.size(); i++)
    for(unsigned int j=0; j<m_pvv_new[i].size(); j++)
      msg_Info()<<i<<j<<" "<<*m_pvv_new[i][j]<<endl;
#endif
}

double Dress_Blob_Base::CalculateBeta(const Vec4D& p) {
  return Vec3D(p).Abs()/p[0];
}

Vec4D Dress_Blob_Base::CalculateMomentumSum(const Particle_Vector& partvec) {
  Vec4D sum = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<partvec.size(); i++) sum+=partvec[i]->Momentum();
  return sum;
}

void Dress_Blob_Base::CalculateWeights() {
  BuildNewParticleVectorVector();

  Weight_Dipole           dip(m_olddipole,m_newdipole,m_softphotons,m_dtype);
  Weight_YFS              yfs(m_newdipole,m_olddipole,m_omegaMin,m_nbar);
  Weight_Jacobian_Mapping jacobM(m_newdipole,m_newspectator,m_olddipole,
                                 m_oldspectator,m_K,m_M,m_u,m_dtype);
  Weight_Jacobian_Lorentz jacobL(m_newdipole,m_newspectator,m_olddipole,
                                 m_oldspectator,m_softphotons,m_dtype);

  double wdipole    = dip.Get();
  double mwdipole   = dip.GetMax();
  double wyfs       = yfs.Get();
  double mwyfs      = yfs.GetMax();
  double wjacobianM = jacobM.Get();
  double mwjacobianM= jacobM.GetMax();
  double wjacobianL = jacobL.Get();
  double mwjacobianL= jacobL.GetMax();
  double whigher    = 1.;
  double mwhigher   = 1.;

  if (!m_soft) {
    Weight_Higher_Order_Corrections c(m_pvv,m_pvv_new,m_dtype);
    whigher    = c.Get();
    mwhigher   = c.GetMax();
  }

  m_genweight       = wdipole * wjacobianM * wjacobianL * whigher * wyfs;
  m_genmaxweight    = mwdipole * mwjacobianM * mwjacobianL * mwhigher * mwyfs;

#ifdef PHOTONS_DEBUG
  char s1w[10],s2w[10],s3w[10],s4w[10],s5w[10];
  sprintf(s1w,"%f",wdipole);
  sprintf(s2w,"%f",wjacobianM);
  sprintf(s3w,"%f",wjacobianL);
  sprintf(s4w,"%f",whigher);
  sprintf(s5w,"%f",wyfs);
  msg_Info()<<"weights:    "<<s1w<<" "<<s2w<<" "<<s3w<<" "<<s4w<<" "<<s5w<<" : "
                            <<m_genweight<<endl;
  char s1m[10],s2m[10],s3m[10],s4m[10],s5m[10];
  sprintf(s1m,"%f",mwdipole);
  sprintf(s2m,"%f",mwjacobianM);
  sprintf(s3m,"%f",mwjacobianL);
  sprintf(s4m,"%f",mwhigher);
  sprintf(s5m,"%f",mwyfs);
  msg_Info()<<"maxweights: "<<s1m<<" "<<s2m<<" "<<s3m<<" "<<s4m<<" "<<s5m
            <<" : "<<m_genmaxweight<<endl;
#endif
  if (Photons::s_strict && m_genweight > m_genmaxweight) {
    msg_Tracking()<<"weight: "<<m_genweight<<" > maxweight: "<<m_genmaxweight
                   <<"  ... event will be rejected. Retrying ... "<<endl;
    for (unsigned int i=0; i<m_softphotons.size(); i++)
      msg_Debugging()<<*m_softphotons[i]<<endl;
    m_genweight = 0.;
  }
}

void Dress_Blob_Base::CheckAvaragePhotonNumberForNumericalErrors() {
  // numerical values slightly < 0. do not lead to fail
  if (m_nbar<0.) {
    m_nbar    = 0.;
    m_success = true;
  }
  // maybe give error if too much below 0., such that physical error
}


void Dress_Blob_Base::DeleteAll(Particle_Vector& pv) {
  while (pv.size() != 0) {
    delete pv[pv.size()-1];
    pv.pop_back();
  }
}

void Dress_Blob_Base::DeleteAllPhotons() {
  DeleteAll(m_softphotons);
}

void Dress_Blob_Base::DetermineU() {
  double M2 = m_M*m_M;
  std::vector<double> mC2;
  std::vector<double> mN2;
  std::vector<Vec3D> q;
  for (unsigned int i=0; i<m_mC.size(); i++)
    mC2.push_back(sqr(m_mC[i]));
  for (unsigned int i=0; i<m_mN.size(); i++)
    mN2.push_back(sqr(m_mN[i]));
  // momenta in Q-CMS
  for (unsigned int i=0; i<m_olddipole.size(); i++)
    q.push_back(Vec3D(m_olddipole[i]->Momentum()));
  for (unsigned int i=0; i<m_oldspectator.size(); i++)
    q.push_back(Vec3D(m_oldspectator[i]->Momentum()));

  double c = Func(M2,mC2,mN2,q,1.);

  if (abs(c) < 1E-12)    m_u = 1.;
  else {
    double i = 0.;
    while ((c*Func(M2,mC2,mN2,q,1.-i) > 0.) && (i<=1))    i = i + 1E-2;
/*    if (i>1)  m_u = -1;
    else */if (abs(Func(M2,mC2,mN2,q,1.-i)) < 1E-14)    m_u = 1.-i;
    else {
      i = i - 1E-2;
      double j = 0.;
      while (c*Func(M2,mC2,mN2,q,1.-i-j) > 0.)    j = j + 1E-4;
      if (abs(Func(M2,mC2,mN2,q,1.-i-j)) < 1E-14)    m_u = 1.-i-j;
      else {
        j = j - 1E-4;
        double k = 0.;
        while (c*Func(M2,mC2,mN2,q,1.-i-j-k) > 0.)    k = k + 1E-6;
        if (abs(Func(M2,mC2,mN2,q,1.-i-j-k)) < 1E-14)    m_u = 1.-i-j-k;
        else {
          k = k - 1E-6;
          double l = 0.;
          while (c*Func(M2,mC2,mN2,q,1.-i-j-k-l) > 0.)    l = l + 1E-8;
          if (abs(Func(M2,mC2,mN2,q,1.-i-j-k-l)) < 1E-14)    m_u = 1.-i-j-k-l;
          else {
            l = l - 1E-8;
            double m = 0.;
            while (c*Func(M2,mC2,mN2,q,1.-i-j-k-l-m) > 0.)   m = m + 1E-10;
            if (abs(Func(M2,mC2,mN2,q,1.-i-j-k-l-m)) < 1E-14)    m_u = 1.-i-j-k-l-m;
            else {
              m = m - 1E-10;
              double n = 0.;
              while (c*Func(M2,mC2,mN2,q,1.-i-j-k-l-m-n) > 0.)    n = n + 1E-12;
              if (abs(Func(M2,mC2,mN2,q,1.-i-j-k-l-m-n)) < 1E-14)    m_u = 1.-i-j-k-l-m-n;
              else {
                n = n - 1E-12;
                double o = 0.;
                while (c*Func(M2,mC2,mN2,q,1.-i-j-k-l-m-n-o) > 0.)    o = o + 1E-14;
                if (abs(Func(M2,mC2,mN2,q,1.-i-j-k-l-m-n-o)) < 1E-14)    m_u = 1.-i-j-k-l-m-n-o;
                else {
                  o = o - 1E-14;
                  double p = 0.;
                  while ((c*Func(M2,mC2,mN2,q,1.-i-j-k-l-m-n-o-p) > 0.) && (p < 1E-14))   p = p + 1E-16;
                  if (p < 1E-14)  m_u = 1.-i-j-k-l-m-n-o-p+0.5E-16;
                  else  m_u = -1.;
                }
              }
            }
          }
        }
      }
    }
  }
}

std::vector<double> Dress_Blob_Base::GenerateNumberAndEnergies() {
//  std::vector<double> R, k;
//  R.push_back(0);
//  for (unsigned int i=0; R.back()<m_nbar; ++i) {
//    R.push_back(R.back()-log(ran->Get()));
//    k.push_back(m_omegaMax*pow(m_omegaMin/m_omegaMax,R.back()/m_nbar));
//    if (k.size()>Photons::s_nmax) break;
//  }
//  k.pop_back();
//  m_n = k.size();
  m_n = -1;
  double expnbar = exp(-m_nbar);
  double prod = 1.;
  while (true) {
    m_n++;
    prod = prod*ran->Get();
    if (prod <= expnbar) break;
  }
  m_n = Min(m_n,Photons::s_nmax);
  std::vector<double> k;
  for (unsigned int i=0; i<m_n; i++) {
    k.push_back(m_omegaMin*pow(m_omegaMax/m_omegaMin,ran->Get()));
  }
  return k;
}

void Dress_Blob_Base::GeneratePhotons(const double& beta1, const double& beta2) {
  m_softphotons.clear();
  std::vector<double> k;
  k = GenerateNumberAndEnergies();
  for (unsigned int i=0; i<m_n; i++) {
    Generate_One_Photon gamma(beta1,beta2,k[i],m_dtype);
    m_softphotons.push_back(gamma.GetPhoton());
  }
}

void Dress_Blob_Base::GeneratePhotons(const IdPairNbarVector& nbars) {
  m_softphotons.clear();
  std::vector<double> k;
  k = GenerateNumberAndEnergies();
  for (unsigned int i=0; i<m_n; i++) {
    Generate_One_Photon gamma(m_olddipole,nbars,k[i],m_dtype);
    m_softphotons.push_back(gamma.GetPhoton());
  }
}

double Dress_Blob_Base::KallenFunction(const double& x, const double& y,
                                       const double& z) {
  return ( x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z );
}

std::string Dress_Blob_Base::ProcessName()
{
  std::string str(" ");
  for (size_t i(0);i<m_chargedinparticles.size();++i)
    str+=m_chargedinparticles[i]->Flav().IDName()+" ";
  for (size_t i(0);i<m_neutralinparticles.size();++i)
    str+=m_neutralinparticles[i]->Flav().IDName()+" ";
  str+=" -> ";
  for (size_t i(0);i<m_chargedoutparticles.size();++i)
    str+=m_chargedoutparticles[i]->Flav().IDName()+" ";
  for (size_t i(0);i<m_neutraloutparticles.size();++i)
    str+=m_neutraloutparticles[i]->Flav().IDName()+" ";
  return str;
}

