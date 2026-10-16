#include "HADRON_RESCATTERING/Main/Deuteron_Production.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Histogram.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/Vector.H"
#include <algorithm>
#include <cmath>
#include <map>

using namespace HADRON_RESCATTERING;
using namespace ATOOLS;

Deuteron_Coalescence::Deuteron_Coalescence(double b_coal,
                                           double p_coal,
                                           double inv_sigma0,
                                           CoalescenceModel model)
  : m_b_coal(b_coal),
    m_p_coal(p_coal),
    m_inv_sigma0(inv_sigma0),
    m_model(model)
{}

std::vector<Deuteron_Channel> Deuteron_Coalescence::Channels() {
  return { Deuteron_Channel::pn_dgamma,   Deuteron_Channel::pn_dpi0,
           Deuteron_Channel::pp_dpiplus,  Deuteron_Channel::nn_dpiminus,
           Deuteron_Channel::pn_dpi0pi0,  Deuteron_Channel::pn_dpipluspiminus,
           Deuteron_Channel::pp_dpipluspi0, Deuteron_Channel::nn_dpiminuspi0 };
}

std::string Deuteron_Coalescence::Name(const Deuteron_Channel & channel) {
  switch (channel) {
  case Deuteron_Channel::pn_dgamma:         return "pn_dgamma";
  case Deuteron_Channel::pn_dpi0:           return "pn_dpi0";
  case Deuteron_Channel::pp_dpiplus:        return "pp_dpiplus";
  case Deuteron_Channel::nn_dpiminus:       return "nn_dpiminus";
  case Deuteron_Channel::pn_dpi0pi0:        return "pn_dpi0pi0";
  case Deuteron_Channel::pn_dpipluspiminus: return "pn_dpipluspiminus";
  case Deuteron_Channel::pp_dpipluspi0:     return "pp_dpipluspi0";
  case Deuteron_Channel::nn_dpiminuspi0:    return "nn_dpiminuspi0";
  default: break;
  }
  return "none";
}

std::vector<Flavour> Deuteron_Coalescence::Products(const Deuteron_Channel & channel,
                                                    const bool & anti) {
  Flavour d(kf_deuterium), gam(kf_photon), pi0(kf_pi);
  Flavour pip(kf_pi_plus), pim(Flavour(kf_pi_plus).Bar());
  std::vector<Flavour> fls;
  switch (channel) {
  case Deuteron_Channel::pn_dgamma:         fls = {d, gam};      break;
  case Deuteron_Channel::pn_dpi0:           fls = {d, pi0};      break;
  case Deuteron_Channel::pp_dpiplus:        fls = {d, pip};      break;
  case Deuteron_Channel::nn_dpiminus:       fls = {d, pim};      break;
  case Deuteron_Channel::pn_dpi0pi0:        fls = {d, pi0, pi0}; break;
  case Deuteron_Channel::pn_dpipluspiminus: fls = {d, pip, pim}; break;
  case Deuteron_Channel::pp_dpipluspi0:     fls = {d, pip, pi0}; break;
  case Deuteron_Channel::nn_dpiminuspi0:    fls = {d, pim, pi0}; break;
  default: break;
  }
  // only conjugate charged states, photon and pi0 are self-conjugate
  if (anti) for (auto & fl : fls) if (fl.IntCharge()!=0) fl = fl.Bar();
  return fls;
}

int Deuteron_Coalescence::PairType(Particle * A,Particle * B) {
  const Flavour & flA = A->Flav(), & flB = B->Flav();
  if (!flA.IsNucleon() || !flB.IsNucleon() || flA.IsAnti()!=flB.IsAnti()) return 0;
  return flA.IsAnti() ? -1 : 1;
}

bool Deuteron_Coalescence::Matches(const Deuteron_Channel & channel,
                                   Particle * A,Particle * B) {
  if (PairType(A,B)==0) return false;
  size_t np = (A->Flav().Kfcode()==kf_p_plus) + (B->Flav().Kfcode()==kf_p_plus);
  switch (channel) {
  case Deuteron_Channel::pn_dgamma:
  case Deuteron_Channel::pn_dpi0:
  case Deuteron_Channel::pn_dpi0pi0:
  case Deuteron_Channel::pn_dpipluspiminus: return np==1;
  case Deuteron_Channel::pp_dpiplus:
  case Deuteron_Channel::pp_dpipluspi0:     return np==2;
  case Deuteron_Channel::nn_dpiminus:
  case Deuteron_Channel::nn_dpiminuspi0:    return np==0;
  default: break;
  }
  return false;
}

double Deuteron_Coalescence::Threshold(const Deuteron_Channel & channel) const {
  double mass = 0.;
  for (const auto & fl : Products(channel,false)) mass += fl.HadMass();
  return mass;
}

double Deuteron_Coalescence::K(const Vec4D & pA,const Vec4D & pB) {
  Poincare boost(pA+pB);
  Vec4D pAcm = pA;
  boost.Boost(pAcm);
  return 2.*Vec3D(pAcm).Abs();
}

double Deuteron_Coalescence::ImpactParameter2(Particle * A,Particle * B) {
  // positions are in mm; in the pair frame A and B move along -+n, so the
  // separation transverse to n does not depend on the production times
  Poincare boost(A->Momentum()+B->Momentum());
  Vec4D xA = A->Position(), xB = B->Position(), pA = A->Momentum();
  boost.Boost(xA); boost.Boost(xB); boost.Boost(pA);
  Vec3D dx = Vec3D(xA)-Vec3D(xB), n = Vec3D(pA);
  if (n.Abs()>0.) n = n/n.Abs();
  double dl = dx*n;
  return std::max(0.,dx.Sqr()-dl*dl)*1.e24;   // mm^2 -> fm^2
}

double Deuteron_Coalescence::Sigma(const Deuteron_Channel & channel,
                                   Particle * A,Particle * B) {
  // The fits are not all zero below threshold (e.g. the d pi pi ones),
  // so enforce the kinematic thresholds here, as in the paper.
  if ((A->Momentum()+B->Momentum()).Mass()<=Threshold(channel)) return 0.;
  switch (channel) {
  case Deuteron_Channel::pn_dgamma:         return m_bb.pnToDGamma(A,B);
  case Deuteron_Channel::pn_dpi0:           return m_bb.pnToDPiZero(A,B);
  case Deuteron_Channel::pp_dpiplus:        return m_bb.ppToDPiPlus(A,B);
  case Deuteron_Channel::nn_dpiminus:       return m_bb.nnToDPiMinus(A,B);
  case Deuteron_Channel::pn_dpi0pi0:        return m_bb.pnToPiZeroPiZero(A,B);
  case Deuteron_Channel::pn_dpipluspiminus: return m_bb.pnToPiPlusPiMinus(A,B);
  case Deuteron_Channel::pp_dpipluspi0:     return m_bb.ppToPiPlusPiZero(A,B);
  case Deuteron_Channel::nn_dpiminuspi0:    return m_bb.nnToPiMinusPiZero(A,B);
  default: break;
  }
  return 0.;
}

std::vector<Channel_Probability>
Deuteron_Coalescence::Probabilities(Particle * A,Particle * B,const double & b2_fm2) {
  std::vector<Channel_Probability> probs;
  if (PairType(A,B)==0) return probs;
  if (m_model==CoalescenceModel::CrossSectionModel) {
    for (const auto & ch : Channels()) {
      if (!Matches(ch,A,B)) continue;
      double sigma = Sigma(ch,A,B);
      if (sigma<=0.) continue;
      // sigma in microbarn, 1/sigma_0 in 1/barn.
      probs.push_back({ch, sigma, std::min(1.,sigma*1.e-6*m_inv_sigma0)});
    }
    return probs;
  }
  // Coalescence-type models: only p n -> d gamma conserves charge.
  if (!Matches(Deuteron_Channel::pn_dgamma,A,B)) return probs;
  // p_coal is a cut on k = |p_p - p_n| in the pair frame, the usual
  // convention (e.g. arXiv:1504.07242, Eq. (1)); b_coal <= 0: no b condition
  double k = K(A->Momentum(),B->Momentum()), P = 0.;
  bool   space = m_b_coal>0.;
  if (m_model==CoalescenceModel::HardDisk) {
    P = ((!space || b2_fm2<m_b_coal*m_b_coal) && k<m_p_coal) ? 1. : 0.;
  }
  else {
    double C_b = (m_b_coal*m_b_coal)/(0.75*M_PI);
    double C_p = (m_p_coal*m_p_coal)/(0.75*M_PI);
    P = 0.75*(space ? std::exp(-b2_fm2/C_b) : 1.)*std::exp(-(k*k)/C_p);
  }
  if (P>0.) probs.push_back({Deuteron_Channel::pn_dgamma, 0., P});
  return probs;
}

Deuteron_Channel Deuteron_Coalescence::SelectChannel(Particle * A,Particle * B,
                                                     const double & b2_fm2) {
  std::vector<Channel_Probability> success;
  double sigsum = 0.;
  for (const auto & cp : Probabilities(A,B,b2_fm2)) {
    if (ran->Get()<cp.prob) { success.push_back(cp); sigsum += cp.sigma; }
  }
  if (success.empty())   return Deuteron_Channel::none;
  if (success.size()==1 || sigsum<=0.) return success[0].channel;
  double r = ran->Get()*sigsum;
  for (const auto & cp : success) {
    if ((r -= cp.sigma)<=0.) return cp.channel;
  }
  return success.back().channel;
}

void Deuteron_Coalescence::Isotropic(const double & p,Vec3D & dir) const {
  double cosT = 1.-2.*ran->Get(), sinT = sqrt(std::max(0.,1.-cosT*cosT));
  double phi  = 2.*M_PI*ran->Get();
  dir = p*Vec3D(sinT*cos(phi),sinT*sin(phi),cosT);
}

void Deuteron_Coalescence::TwoBody(const Vec4D & P,const double & m1,const double & m2,
                                   Vec4D & p1,Vec4D & p2) const {
  // isotropic two-body decay in the rest frame of P, boosted back with P
  double M2 = P.Abs2(), M = sqrt(M2);
  double lambda = (M2-sqr(m1+m2))*(M2-sqr(m1-m2));
  double p = sqrt(std::max(0.,lambda))/(2.*M);
  Vec3D dir;
  Isotropic(p,dir);
  p1 = Vec4D(sqrt(p*p+m1*m1), dir);
  p2 = Vec4D(sqrt(p*p+m2*m2),-dir);
  Poincare boost(P);
  boost.BoostBack(p1);
  boost.BoostBack(p2);
}

bool Deuteron_Coalescence::FinalState(const Deuteron_Channel & channel,
                                      const Vec4D & pA,const Vec4D & pB,
                                      std::vector<Vec4D> & moms) const {
  std::vector<Flavour> fls = Products(channel,false);
  Vec4D  P     = pA+pB;
  double s     = P.Abs2(), sqrtS = sqrt(std::max(0.,s));
  moms.assign(fls.size(),Vec4D());
  if (fls.empty() || sqrtS<=Threshold(channel)) return false;
  if (fls.size()==2) {
    TwoBody(P,fls[0].HadMass(),fls[1].HadMass(),moms[0],moms[1]);
    return true;
  }
  // Three-body final state d pi pi, Sec. III D: flat phase space, i.e. a
  // uniform point (m_pipi^2, m_dpi^2) inside the Dalitz region.  The
  // deuteron energy is fixed by the invariant mass of the other two,
  //    E_d = (s + m_d^2 - m_pipi^2)/(2 sqrt(s)),
  // (Eq. (15) of the paper is printed with m_dpi^2 here, which is a typo).
  // For given m_pipi, flat phase space means an isotropic pi pi decay in
  // their rest frame, so only m_pipi^2 has to be taken from the Dalitz point.
  double m1 = fls[0].HadMass(), m2 = fls[1].HadMass(), m3 = fls[2].HadMass();
  double m23min = sqr(m2+m3), m23max = sqr(sqrtS-m1);
  double m12min = sqr(m1+m2), m12max = sqr(sqrtS-m3);
  double m23sq  = -1.;
  for (size_t trial=0;trial<10000;trial++) {
    double m23 = m23min+ran->Get()*(m23max-m23min);
    double m12 = m12min+ran->Get()*(m12max-m12min);
    double M23 = sqrt(m23);
    double E2  = (m23-m3*m3+m2*m2)/(2.*M23), E1 = (s-m23-m1*m1)/(2.*M23);
    double q2  = sqrt(std::max(0.,E2*E2-m2*m2)), q1 = sqrt(std::max(0.,E1*E1-m1*m1));
    double lo  = sqr(E1+E2)-sqr(q1+q2), hi = sqr(E1+E2)-sqr(q1-q2);
    if (m12>=lo && m12<=hi) { m23sq = m23; break; }
  }
  if (m23sq<0.) return false;
  double E_d = (s+m1*m1-m23sq)/(2.*sqrtS);
  double p_d = sqrt(std::max(0.,E_d*E_d-m1*m1));
  Vec3D dir;
  Isotropic(p_d,dir);
  Vec4D pd_cm(E_d,dir), ppipi_cm(sqrtS-E_d,-dir);
  Vec4D p2_cm, p3_cm;
  TwoBody(ppipi_cm,m2,m3,p2_cm,p3_cm);
  Poincare boost(P);
  boost.BoostBack(pd_cm);
  boost.BoostBack(p2_cm);
  boost.BoostBack(p3_cm);
  moms[0] = pd_cm; moms[1] = p2_cm; moms[2] = p3_cm;
  return true;
}

Blob * Deuteron_Coalescence::BuildBlob(const Deuteron_Channel & channel,
                                       Particle * A,Particle * B) {
  std::vector<Vec4D> moms;
  if (!FinalState(channel,A->Momentum(),B->Momentum(),moms)) {
    msg_Out()<<METHOD<<": "<<Name(channel)<<" closed for sqrt(s) = "
             <<(A->Momentum()+B->Momentum()).Mass()<<", skipping.\n";
    return nullptr;
  }
  std::vector<Flavour> fls = Products(channel,PairType(A,B)<0);
  Vec4D pos = 0.5*(A->Position()+B->Position());
  Blob * blob = new Blob();
  blob->SetType(btp::Soft_Collision);
  blob->SetTypeSpec("Deuteron_Formation: "+Name(channel));
  blob->SetStatus(blob_status::needs_hadrondecays);
  blob->SetPosition(pos);
  blob->SetId();
  for (Particle * part : {A,B}) {
    part->SetStatus(part_status::decayed);
    blob->AddToInParticles(part);
  }
  for (size_t i=0;i<fls.size();i++) {
    Particle * out = new Particle(-1,fls[i],moms[i],'F');
    out->SetNumber();
    out->SetPosition(pos);
    blob->AddToOutParticles(out);
  }
  return blob;
}

Blob * Deuteron_Coalescence::operator()(Particle * A,Particle * B,
                                        const double & b2_fm2) {
  Deuteron_Channel channel = SelectChannel(A,B,b2_fm2);
  if (channel==Deuteron_Channel::none) return nullptr;
  return BuildBlob(channel,A,B);
}

//////////////////////////////////////////////////////////////////////////////
// Tests: probabilities, kinematics/conservation, and channel selection.
//////////////////////////////////////////////////////////////////////////////

namespace {
  Flavour NucleonA(const Deuteron_Channel & ch) {
    std::string name = Deuteron_Coalescence::Name(ch);
    return Flavour(name[0]=='p' ? kf_p_plus : kf_n);
  }
  Flavour NucleonB(const Deuteron_Channel & ch) {
    std::string name = Deuteron_Coalescence::Name(ch);
    return Flavour(name[1]=='p' ? kf_p_plus : kf_n);
  }
  // Back-to-back pair with c.m. momentum difference k = |pA - pB| along z
  // (the paper's k), then boosted with velocity beta along a fixed direction
  // (beta=0: pair at rest).
  void MakePair(const Flavour & flA,const Flavour & flB,const double & k,
                const double & beta,Vec4D & pA,Vec4D & pB) {
    double q = k/2.;
    pA = Vec4D(sqrt(q*q+sqr(flA.HadMass())),0.,0., q);
    pB = Vec4D(sqrt(q*q+sqr(flB.HadMass())),0.,0.,-q);
    if (beta<=0.) return;
    double M = (pA+pB).Mass(), gamma = 1./sqrt(1.-beta*beta);
    Vec3D  n = Vec3D(0.3,-0.5,0.8)/Vec3D(0.3,-0.5,0.8).Abs();
    Vec4D  Q(gamma*M,gamma*beta*M*n);
    Poincare boost(Q);
    boost.BoostBack(pA);
    boost.BoostBack(pB);
  }
}

void Deuteron_Coalescence::Tests(const std::string & dir) {
  MakeDir(dir,true);
  msg_Out()<<"=== "<<METHOD<<": writing test output to "<<dir<<"/\n";

  // (1) sigma(k) and P(k) for each channel, with the right nucleon flavours.
  size_t nbins = 400;
  double kmin  = 0.01, kmax = 5.0, dk = (kmax-kmin)/double(nbins);
  for (const auto & ch : Channels()) {
    Histogram hsig(0,kmin,kmax,nbins), hP(0,kmin,kmax,nbins);
    for (size_t i=0;i<nbins;i++) {
      double k = kmin+(i+0.5)*dk;
      Vec4D pA, pB;
      MakePair(NucleonA(ch),NucleonB(ch),k,0.,pA,pB);
      Particle A(-1,NucleonA(ch),pA), B(-1,NucleonB(ch),pB);
      double sigma = Sigma(ch,&A,&B);
      hsig.Insert(k,sigma);
      hP.Insert(k,std::min(1.,sigma*1.e-6*m_inv_sigma0));
    }
    hsig.Output(dir+"/test_sigma_"+Name(ch)+".dat");
    hP.Output(dir+"/test_P_"+Name(ch)+".dat");
  }

  // (2) Conservation and on-shell checks, nucleons and antinucleons, pair at
  //     rest and boosted.
  msg_Out()<<"--- kinematics: max |sum p_out - sum p_in|/sqrt(s), "
           <<"max |p^2-m^2| [GeV^2], charge violations\n";
  bool allgood = true;
  for (const auto & ch : Channels()) {
    double maxdp = 0., maxms = 0.;
    size_t nqviol = 0, nclosed = 0, ntried = 0;
    for (bool anti : {false,true}) {
      Flavour flA = NucleonA(ch), flB = NucleonB(ch);
      if (anti) { flA = flA.Bar(); flB = flB.Bar(); }
      std::vector<Flavour> fls = Products(ch,anti);
      int qin = flA.IntCharge()+flB.IntCharge(), qout = 0;
      for (const auto & fl : fls) qout += fl.IntCharge();
      for (double k : {0.4, 1.2, 2.0, 3.0, 5.0}) {
        for (double beta : {0., 0.9}) {
          Vec4D pA, pB;
          MakePair(flA,flB,k,beta,pA,pB);
          double sqrtS = (pA+pB).Mass();
          for (size_t n=0;n<2000;n++) {
            ntried++;
            std::vector<Vec4D> moms;
            if (!FinalState(ch,pA,pB,moms)) { nclosed++; continue; }
            Vec4D sum;
            for (size_t i=0;i<moms.size();i++) {
              sum   += moms[i];
              maxms  = std::max(maxms,std::abs(moms[i].Abs2()-sqr(fls[i].HadMass())));
            }
            Vec4D diff = sum-pA-pB;
            for (size_t mu=0;mu<4;mu++)
              maxdp = std::max(maxdp,std::abs(diff[mu])/sqrtS);
            if (qin!=qout) nqviol++;
          }
        }
      }
    }
    bool good = (maxdp<1.e-10 && maxms<1.e-8 && nqviol==0);
    allgood  &= good;
    msg_Out()<<"   "<<Name(ch)<<": "<<maxdp<<", "<<maxms<<", "<<nqviol
             <<"   (closed below threshold: "<<nclosed<<"/"<<ntried<<")"
             <<(good ? "   OK" : "   FAILED")<<"\n";
  }
  msg_Out()<<"--- kinematics "<<(allgood ? "all OK" : "FAILED")<<"\n";

  // (3) Deuteron c.m. momentum for d pi pi at sqrt(s) for k = 3 GeV, vs.
  //     flat three-body phase space, and cos(theta) of the deuteron
  //     (isotropy) for all channels.
  for (const auto & ch : Channels()) {
    std::vector<Flavour> fls = Products(ch,false);
    Vec4D pA, pB;
    MakePair(NucleonA(ch),NucleonB(ch),3.,0.,pA,pB);
    double s = (pA+pB).Abs2(), sqrtS = sqrt(s), m1 = fls[0].HadMass();
    double pmax = sqrt(sqr((s+sqr(m1)-sqr(Threshold(ch)-m1))/(2.*sqrtS))-sqr(m1));
    size_t nev = 200000, nb = 50;
    Histogram hpd(1,0.,1.05*pmax,nb), hcos(1,-1.,1.,20);
    for (size_t n=0;n<nev;n++) {
      std::vector<Vec4D> moms;
      if (!FinalState(ch,pA,pB,moms)) continue;
      Vec3D pd(moms[0]);
      hpd.Insert(pd.Abs(),1./(nev*hpd.BinSize()));
      hcos.Insert(pd[3]/pd.Abs(),1./(nev*hcos.BinSize()));
    }
    hcos.Output(dir+"/test_costheta_d_"+Name(ch)+".dat");
    if (fls.size()<3) continue;
    hpd.Output(dir+"/test_pd_cm_"+Name(ch)+".dat");
    // Expected: dN/dm23^2 ~ (m12^2_max - m12^2_min)(m23^2), with
    // m23^2 = s + m1^2 - 2 sqrt(s) E_d and |dm23^2/dp_d| = 2 sqrt(s) p_d/E_d.
    double m2 = fls[1].HadMass(), m3 = fls[2].HadMass();
    Histogram hexp(0,0.,1.05*pmax,nb);
    std::vector<double> vals(nb,0.);
    double norm = 0.;
    for (size_t i=0;i<nb;i++) {
      double pd  = (i+0.5)*hexp.BinSize(), Ed = sqrt(pd*pd+m1*m1);
      double m23 = s+m1*m1-2.*sqrtS*Ed;
      if (m23<=sqr(m2+m3)) continue;
      double M23 = sqrt(m23);
      double E2  = (m23-m3*m3+m2*m2)/(2.*M23), E1 = (s-m23-m1*m1)/(2.*M23);
      double q2  = sqrt(std::max(0.,E2*E2-m2*m2)), q1 = sqrt(std::max(0.,E1*E1-m1*m1));
      vals[i]    = 4.*q1*q2*(2.*sqrtS*pd/Ed);
      norm      += vals[i]*hexp.BinSize();
    }
    for (size_t i=0;i<nb;i++)
      if (norm>0.) hexp.Insert((i+0.5)*hexp.BinSize(),vals[i]/norm);
    hexp.Output(dir+"/test_pd_cm_"+Name(ch)+"_phasespace.dat");
  }

  // (4) Channel selection: observed fractions vs. the exact expectation of
  //     the paper's rule (one random number per channel, ties broken by
  //     relative cross section).  Use a large 1/sigma_0 for visible rates.
  msg_Out()<<"--- channel selection (1/sigma_0 = 1000/b): "
           <<"pair, k, channel, expected, observed, pull\n";
  Deuteron_Coalescence test(m_b_coal,m_p_coal,1000.,CoalescenceModel::CrossSectionModel);
  size_t ntrials = 200000;
  double maxpull = 0.;
  for (const auto & pair : std::vector<std::pair<kf_code,kf_code> >
         {{kf_p_plus,kf_n},{kf_p_plus,kf_p_plus},{kf_n,kf_n}}) {
    for (double k : {0.6, 1.1, 1.6}) {
      Vec4D pA, pB;
      Flavour flA(pair.first), flB(pair.second);
      MakePair(flA,flB,k,0.,pA,pB);
      Particle A(-1,flA,pA), B(-1,flB,pB);
      std::vector<Channel_Probability> probs = test.Probabilities(&A,&B);
      size_t nch = probs.size();
      std::vector<double> expected(nch,0.);
      for (size_t mask=1;mask<(size_t(1)<<nch);mask++) {
        double weight = 1., sigsum = 0.;
        for (size_t j=0;j<nch;j++) {
          bool in = mask&(size_t(1)<<j);
          weight *= in ? probs[j].prob : 1.-probs[j].prob;
          if (in) sigsum += probs[j].sigma;
        }
        for (size_t j=0;j<nch;j++)
          if (mask&(size_t(1)<<j)) expected[j] += weight*probs[j].sigma/sigsum;
      }
      std::map<Deuteron_Channel,size_t> counts;
      for (size_t n=0;n<ntrials;n++) counts[test.SelectChannel(&A,&B)]++;
      for (size_t j=0;j<nch;j++) {
        double obs  = double(counts[probs[j].channel])/double(ntrials);
        double err  = sqrt(expected[j]*(1.-expected[j])/double(ntrials));
        double pull = err>0. ? (obs-expected[j])/err : 0.;
        maxpull = std::max(maxpull,std::abs(pull));
        msg_Out()<<"   "<<flA<<" "<<flB<<", k = "<<k<<", "<<Name(probs[j].channel)
                 <<": "<<expected[j]<<", "<<obs<<", "<<pull<<"\n";
      }
    }
  }
  msg_Out()<<"--- channel selection: max |pull| = "<<maxpull
           <<(maxpull<5. ? "   OK" : "   SUSPICIOUS")<<"\n";

  // (5) coalescence models: impact parameter and probabilities against
  //     hand-calculated values
  {
    const double fm = 1.e-12;                  // 1 fm in mm
    Flavour flp(kf_p_plus), fln(kf_n);
    Vec4D pA, pB;
    MakePair(flp,fln,0.1,0.,pA,pB);            // k = 0.1 GeV along z, at rest
    Particle A(-1,flp,pA), B(-1,fln,pB);
    // B 2 fm away transverse, 5 fm along z and produced 3 fm/c later: b = 2 fm
    A.SetPosition(Vec4D(0.,0.,0.,0.));
    B.SetPosition(Vec4D(3.*fm,2.*fm,0.,5.*fm));
    double b2 = ImpactParameter2(&A,&B);
    // the same configuration boosted with beta = 0.9
    double M = (pA+pB).Mass(), g = 1./std::sqrt(1.-0.81);
    Vec3D  n = Vec3D(0.3,-0.5,0.8)/Vec3D(0.3,-0.5,0.8).Abs();
    Poincare lab(Vec4D(g*M,0.9*g*M*n));
    Vec4D pAb = pA, pBb = pB, xA = A.Position(), xB = B.Position();
    lab.BoostBack(pAb); lab.BoostBack(pBb); lab.BoostBack(xA); lab.BoostBack(xB);
    Particle Ab(-1,flp,pAb), Bb(-1,fln,pBb);
    Ab.SetPosition(xA); Bb.SetPosition(xB);
    double b2boost = ImpactParameter2(&Ab,&Bb);
    // HardDisk: p_coal = 0.2 GeV, b_coal = 3 fm; Exponential with the same
    Deuteron_Coalescence hd(3.,0.2,1.,CoalescenceModel::HardDisk);
    Deuteron_Coalescence ex(3.,0.2,1.,CoalescenceModel::Exponential);
    Deuteron_Coalescence hdk(0.,0.2,1.,CoalescenceModel::HardDisk);   // no b cut
    double Pin  = hd.Probabilities(&A,&B,4.)[0].prob;                // k=0.1<0.2, b=2<3
    double Pout = hd.Probabilities(&A,&B,16.).empty() ? 0. : 1.;     // b=4>3
    double Pnob = hdk.Probabilities(&A,&B,1.e6)[0].prob;             // no b cut
    Vec4D qA, qB; MakePair(flp,fln,0.3,0.,qA,qB);
    Particle C(-1,flp,qA), D(-1,fln,qB);
    double Pk   = hd.Probabilities(&C,&D,0.).empty() ? 0. : 1.;      // k=0.3>0.2
    double Pex  = ex.Probabilities(&A,&B,4.)[0].prob;
    double Pexp = 0.75*std::exp(-4./(9./(0.75*M_PI)))*std::exp(-0.01/(0.04/(0.75*M_PI)));
    bool good = std::abs(b2-4.)<1.e-6 && std::abs(b2boost-4.)<1.e-6 &&
                Pin==1. && Pout==0. && Pnob==1. && Pk==0. && std::abs(Pex/Pexp-1.)<1.e-9;
    msg_Out()<<"--- coalescence: b = "<<std::sqrt(b2)<<" fm (2), boosted "<<std::sqrt(b2boost)
             <<" fm; HardDisk in/out(b)/no-b-cut/out(k) = "<<Pin<<"/"<<Pout<<"/"<<Pnob<<"/"<<Pk
             <<" (1/0/1/0); Exponential = "<<Pex<<" ("<<Pexp<<")"
             <<(good ? "   OK" : "   FAILED")<<"\n";
  }
}
