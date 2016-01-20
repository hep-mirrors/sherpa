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

#include <iomanip>

#ifdef PHOTONS_DEBUG
#include "ATOOLS/Math/Histogram_2D.H"
#endif

using namespace PHOTONS;
using namespace ATOOLS;
using namespace std;

Dress_Blob_Base::Dress_Blob_Base() :
  m_photonsadded(false), m_success(false),
  m_soft((PHOTONS::Photons::s_mode==yfsmode::soft)?true:false),
  p_newinitialstate(NULL), m_accu(PHOTONS::Photons::s_accu),
  m_genweight(1.), m_genmaxweight(1.), m_nbar(0.), m_n(0.),
  m_omegaMax(PHOTONS::Photons::s_uvcutoff),
  m_omegaMin(PHOTONS::Photons::s_ircutoff)
{}

Dress_Blob_Base::~Dress_Blob_Base()
{
  if (m_dtype == Dipole_Type::ff)  delete p_newinitialstate;
}

void Dress_Blob_Base::BuildNewParticleVectorVector()
{
  if ((m_dtype == Dipole_Type::ff) && (p_newinitialstate != NULL))
    delete p_newinitialstate;
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
  msg_Debugging()<<"new particles:\n";
  for (unsigned int i=0; i<m_pvv_new.size(); i++)
    for(unsigned int j=0; j<m_pvv_new[i].size(); j++)
      msg_Debugging()<<i<<j<<" "<<*m_pvv_new[i][j]<<endl;
#endif
}

double Dress_Blob_Base::CalculateBeta(const Vec4D& p)
{
  return Vec3D(p).Abs()/p[0];
}

Vec4D Dress_Blob_Base::CalculateMomentumSum(const Particle_Vector& partvec)
{
  Vec4D sum = Vec4D(0.,0.,0.,0.);
  for (unsigned int i=0; i<partvec.size(); i++) sum+=partvec[i]->Momentum();
  return sum;
}

void Dress_Blob_Base::CalculateWeights()
{
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
  if (IsNan(wdipole) || IsNan(mwdipole) ||
      IsNan(wjacobianM) || IsNan(mwjacobianM) ||
      IsNan(wjacobianL) || IsNan(mwjacobianL) ||
      IsNan(whigher) || IsNan(mwhigher) ||
      IsNan(wyfs) || IsNan(mwyfs) ||
      IsNan(m_genweight) || IsNan(m_genmaxweight) || IsNan((double)m_n)) {
    msg_Error()<<METHOD<<"(): Encountered nan weight\n"
               <<"weights:    "
               <<std::setw(12)<<wdipole
               <<std::setw(12)<<wjacobianM
               <<std::setw(12)<<wjacobianL
               <<std::setw(12)<<whigher
               <<std::setw(12)<<wyfs<<" : "
               <<std::setw(12)<<m_genweight<<std::endl
               <<"maxweights: "
               <<std::setw(12)<<mwdipole
               <<std::setw(12)<<mwjacobianM
               <<std::setw(12)<<mwjacobianL
               <<std::setw(12)<<mwhigher
               <<std::setw(12)<<mwyfs<<" : "
               <<std::setw(12)<<m_genmaxweight<<std::endl;
    wdipole=wjacobianM=wjacobianL=whigher=wyfs=0.;
    mwdipole=mwjacobianM=mwjacobianL=mwhigher=mwyfs=1.;
    m_genweight=0.; m_genmaxweight=1.;
  }
  Photons::s_histo_dipole.Insert(wdipole/mwdipole,m_n);
  Photons::s_histo_jacobianM.Insert(wjacobianM/mwjacobianM,m_n);
  Photons::s_histo_jacobianL.Insert(wjacobianL/mwjacobianL,m_n);
  Photons::s_histo_higher.Insert(whigher/mwhigher,m_n);
  Photons::s_histo_yfs.Insert(wyfs/mwyfs,m_n);
  Photons::s_histo_total.Insert(m_genweight/m_genmaxweight,m_n);
  Vec4D sump(0.,0.,0.,0.);
  for (size_t i(0);i<m_softphotons.size();++i)
   sump+=m_softphotons[i]->Momentum();
  Poincare dec(m_newdipole[0]->Momentum());
  dec.Boost(sump);
  Vec4D test(m_newdipole[0]->Momentum()); dec.Boost(test);
  double t(sump[0]);
  Photons::s_histo_t_dipole.Insert(wdipole/mwdipole,t);
  Photons::s_histo_t_jacobianM.Insert(wjacobianM/mwjacobianM,t);
  Photons::s_histo_t_jacobianL.Insert(wjacobianL/mwjacobianL,t);
  Photons::s_histo_t_higher.Insert(whigher/mwhigher,t);
  Photons::s_histo_t_yfs.Insert(wyfs/mwyfs,t);
  Photons::s_histo_t_total.Insert(m_genweight/m_genmaxweight,t);

  msg_Debugging()<<"weights:    "<<std::setw(12)<<wdipole
                                 <<std::setw(12)<<wjacobianM
                                 <<std::setw(12)<<wjacobianL
                                 <<std::setw(12)<<whigher
                                 <<std::setw(12)<<wyfs<<" : "
                                 <<std::setw(12)<<m_genweight<<endl;
  msg_Debugging()<<"maxweights: "<<std::setw(12)<<mwdipole
                                 <<std::setw(12)<<mwjacobianM
                                 <<std::setw(12)<<mwjacobianL
                                 <<std::setw(12)<<mwhigher
                                 <<std::setw(12)<<mwyfs<<" : "
                                 <<std::setw(12)<<m_genmaxweight<<endl;
  msg_Debugging()<<"max. weight factor: "<<Photons::s_increasemaxweight<<endl;
#endif
  if (Photons::s_strict &&
      m_genweight > m_genmaxweight*Photons::s_increasemaxweight) {
    msg_Tracking()<<"weight: "<<m_genweight<<" > maxweight: "
                  <<m_genmaxweight*Photons::s_increasemaxweight
                  <<"  ... event will be rejected. Retrying ... "<<endl;
    for (unsigned int i=0; i<m_softphotons.size(); i++)
      msg_Debugging()<<*m_softphotons[i]<<endl;
    m_genweight = 0.;
  }
}

void Dress_Blob_Base::CheckAvaragePhotonNumberForNumericalErrors()
{
  // numerical values slightly < 0. do not lead to fail
  if (m_nbar<0.) {
    m_nbar    = 0.;
    m_success = true;
  }
  // maybe give error if too much below 0., such that physical error
}


void Dress_Blob_Base::DeleteAll(Particle_Vector& pv)
{
  while (pv.size() != 0) { delete pv[pv.size()-1]; pv.pop_back(); }
}

void Dress_Blob_Base::DetermineU()
{
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
    if (abs(Func(M2,mC2,mN2,q,1.-i)) < 1E-14)    m_u = 1.-i;
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
  msg_Debugging()<<"u:    "<<m_u<<std::endl;
}

std::vector<double> Dress_Blob_Base::GenerateNumberAndEnergies()
{
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
  m_n = Min(Max(m_n,Photons::s_nmin),Photons::s_nmax);
  std::vector<double> k;
  for (unsigned int i=0; i<m_n; i++) {
    k.push_back(m_omegaMin*pow(m_omegaMax/m_omegaMin,ran->Get()));
  }
  if (m_n && msg_LevelIsDebugging()) {
    msg_Debugging()<<"k_i = ["<<m_omegaMin<<" .. "<<m_omegaMax<<"]\n";
    for (size_t i(0);i<k.size();++i) {
      msg_Debugging()<<"  "<<std::setw(2)<<i<<": "<<std::setw(8)<<k[i]<<std::endl;
    }
  }
  return k;
}

void Dress_Blob_Base::GeneratePhotons(const double& beta1, const double& beta2)
{
  m_softphotons.clear();
  std::vector<double> k;
  k = GenerateNumberAndEnergies();
  for (unsigned int i=0; i<m_n; i++) {
    Generate_One_Photon gamma(beta1,beta2,k[i],m_dtype);
    m_softphotons.push_back(gamma.GetPhoton());
  }
}

void Dress_Blob_Base::GeneratePhotons(const IdPairNbarVector& nbars)
{
  m_softphotons.clear();
  std::vector<double> k;
  k = GenerateNumberAndEnergies();
  for (unsigned int i=0; i<m_n; i++) {
    Generate_One_Photon gamma(m_olddipole,nbars,k[i],m_dtype);
    m_softphotons.push_back(gamma.GetPhoton());
  }
}

double Dress_Blob_Base::KallenFunction(const double& x, const double& y,
                                       const double& z)
{
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

