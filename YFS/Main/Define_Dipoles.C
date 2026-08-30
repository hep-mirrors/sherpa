#include "YFS/Main/Define_Dipoles.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"
#include <algorithm>
#include <iomanip>
#include <limits>
#include <set>
#ifdef USING__LOOPTOOLS
  #include "clooptools.h"
#endif

using namespace YFS;
using namespace ATOOLS;
using namespace std;


Define_Dipoles::Define_Dipoles() {
  m_in = 2; // This is fine in YFS. It will not work for any other inital state multiplicity
  m_softphotonSum = {0., 0., 0., 0.};
  p_yfsFormFact = new YFS::YFS_Form_Factor();
}

Define_Dipoles::~Define_Dipoles() {
  if (m_pole_tried > 0)
    msg_Info() << "YFS WW pole expansion: " << m_pole_built << " / "
               << m_pole_tried << " events inside the leading-pole window ("
               << (100.*m_pole_built)/m_pole_tried << " %); the rest used the "
               << "flat four-fermion dipoles.\n";
  if (m_recoil_tried > 0)
    msg_Info() << "YFS WW production recoil: " << m_recoil_massfail << " / "
               << m_recoil_tried << " rejected on the W-mass check.\n";
  if(p_yfsFormFact) delete p_yfsFormFact;
}


void Define_Dipoles::MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  m_N_born_Gamma=1;
  for(auto f: fl) if(f.IsPhoton()) m_N_born_Gamma+=1;
  if(!HasISR()) return;
  m_test_dip.clear();
  m_flav_label.clear();
  for(size_t i(0); i<fl.size(); ++i) m_flav_label[fl[i]] = i;
  m_softphotonSum *= 0;
  m_out = fl.size() - m_in;
  m_bornmomenta = born;
  m_set.BuildInitial(fl, mom, born, m_alpha);
  // BuildInitial leaves no dipole if the initial state carries no charge.
  if (m_set.HasII()) { m_g = m_set.II().m_gamma; m_gp = m_set.II().m_gammap; }
}


void Define_Dipoles::MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born) {
  if (m_mode == yfsmode::fsr) return;
  if (!HasFSR()) return;
  if (mom.size() != fl.size())
    THROW(fatal_error, "Incorrect dipole size in YFS for initial-final dipoles");
  m_set.BuildIF(fl, mom, born, m_alpha);
}

void Define_Dipoles::MakeDipolesFF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  // Never built FF dipoles - only the charged/neutral bookkeeping, which is
  // reported by operator<< and used nowhere else. MakeDipoles builds them.
  Dipole_FF(fl, mom);
}

void Define_Dipoles::MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born ) {
  if ((mom.size() != fl.size()) || (born.size() != fl.size()))
    THROW(fatal_error, "Incorrect dipole size in YFS for final-state dipoles");
  m_test_dip.clear();
  m_flav_label.clear();
  for(size_t i(0); i<fl.size(); ++i) m_flav_label[fl[i]] = i;
  m_bornmomenta = born;
  m_out = fl.size() - m_in;
  if (!HasFSR()) return;
  m_set.BuildFinal(fl, mom, born, m_alpha,
                   [this](const YFS::Dipole &D) {
                     return ResonanceWidthDistance(const_cast<YFS::Dipole&>(D));
                   });
  Dipole_FF(fl, mom);
}

bool Define_Dipoles::MakeDipolesPole(ATOOLS::Flavour_Vector const &fl,
                                     ATOOLS::Vec4D_Vector const &mom,
                                     ATOOLS::Vec4D_Vector const &born) {
  if (m_wwscheme != wwscheme::pole) return false;
  if (!HasFSR()) return false;
  const bool built(m_set.BuildPole(fl, mom, born, m_alpha, m_wwonshell != 0,
                                   m_resonace_max, m_wwpoleemission != 0));
  ++m_pole_tried;
  if (built) ++m_pole_built;
  m_pole_active = built;
  return built;
}


bool Define_Dipoles::ApplyPoleRecoil(ATOOLS::Vec4D_Vector &plab) {
  if (!m_pole_active) return false;
  ++m_recoil_tried;
  YFS::DipoleView ff(m_set.FF());
  if (ff.size() != 1) return false;
  Dipole &D(ff[0]);
  const DipoleSet::WWLegs &w(m_set.WW());

  // BuildPole adds the pair as (W-, W+), so leg 0 carries w.lm/w.nm and leg 1
  // carries w.lp/w.np.
  const Vec4D oldw[2] = {D.GetMomenta(0),    D.GetMomenta(1)};
  const Vec4D neww[2] = {D.GetNewMomenta(0), D.GetNewMomenta(1)};
  const std::size_t dau[2][2] = {{w.lm, w.nm}, {w.lp, w.np}};

  for (int i(0); i < 2; ++i) {
    const double m2o(oldw[i].Abs2()), m2n(neww[i].Abs2());
    if (m2o <= 0. || m2n <= 0.) return false;
    if (!IsEqual(sqrt(m2o), sqrt(m2n), 1e-6)) {
      ++m_recoil_massfail;
      if (m_recoil_massfail < 4)
        msg_Error() << METHOD << ": W mass not preserved by the dipole recoil, "
                    << sqrt(m2o) << " -> " << sqrt(m2n)
                    << " GeV; rejecting the event.\n";
      return false;
    }
    Poincare toRest(oldw[i]), fromNew(neww[i]);
    for (int j(0); j < 2; ++j) {
      Vec4D p(plab[dau[i][j]]);
      toRest.Boost(p);        // into the pre-emission W rest frame
      fromNew.BoostBack(p);   // out along the post-emission W
      plab[dau[i][j]] = p;
    }
  }
  return true;
}











void Define_Dipoles::Get4Mom(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector mom) {
  Vec4D_Vector P;
  for(size_t i = 2; i < fl.size(); ++i)
  {
    if (fl[i].IntCharge()!=0) {
      m_test_dip[fl[i]] = mom[i];
      P.push_back(mom[i]);
      if (P.size() == 2) break;
    }
  }
  if (P.size() != 2) {
    PRINT_VAR(P.size());
    THROW(fatal_error, "Wrong size dipole");
  }
}





void Define_Dipoles::Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  CleanOutParticles();
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_FF Momenta");
  }
  for(size_t i = 0; i < fl.size(); ++i)
  {
    if (fl[i].IsQED()) {
      m_chargedoutparticles.push_back(mom[i]);
      m_massOutC.push_back(mom[i].Mass());

    }
    else {
      m_neutraloutparticles.push_back(mom[i]);
      m_massOutN.push_back(mom[i].Mass());
    }
  }
}







double Define_Dipoles::CalculateRealSub(const Vec4D &k) {
  double sub(0);
  // if(FixedOrder()!=fixed_order::full) return sub;
  Vec4D eik{0.,0.,0.,0.};
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
       Vec4D p = D.GetMomenta(i);
      eik += D.m_Q[i]*p/(p*k);
    }
  }
  for (auto &D : m_set.FF()) {
     if(!D.IsResonance()) continue;
    for(size_t i = 0; i < D.GetBornMomenta().size(); ++i)
    {
      Vec4D p = D.GetMomenta(i);
      eik += -D.m_Q[i]*p/(p*k);
    }
  }
  sub = -m_alpha / (4 * M_PI * M_PI)*eik*eik;
  return sub/(m_N_born_Gamma!=0?m_N_born_Gamma:1.0);
}

double Define_Dipoles::CalculateRealSubIF(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_set.IF()){
    if(m_massless_sub) sub += D.EikonalMassless(k, D.GetMomenta(0), D.GetMomenta(1));
    else sub +=  D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  }
  return sub;
}


double Define_Dipoles::CalculateVirtualSub() {
  double sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  if(m_dim_reg==1) return CalculateVirtualSubEps();
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  for (auto &D : m_set.FF()) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }

  for (auto &D : m_set.IF()){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  return sub;
}

double Define_Dipoles::CalculateVirtualSubEps() {
  DivArrD sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  for (auto &D : m_set.FF()) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }

  for (auto &D : m_set.IF()){
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateVVSubEps() {
  DivArrD sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  for (auto &D : m_set.FF()) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }

  for (auto &D : m_set.IF()){
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  m_vvSub=0.5*sub*sub;
  return (0.5*sub*sub).Finite();
}

double Define_Dipoles::CalculateRealVirtualSubEps(const Vec4D &k) {
  // Virtual YFS B-hat in dim-reg at the reduced (photon-removed) kinematics.
  // Must use the same BVV function as CalculateVirtualSubEps: the RV finite
  // remainder is defined by subtracting 2*alpha*B_fin off the one-loop real
  // emission ME, with the identical finite convention as the Born virtual —
  // otherwise beta_1^1 does not vanish in the soft limit.
  DivArrD sub(0);
  if(m_tchannel>=2) return CalculateVirtualSubTchannel();
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_set.FF()) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }

  for (auto &D : m_set.IF()){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }
  m_virtSub=sub;
  return sub.Finite();
}



double Define_Dipoles::FormFactorSum(){
  double form = 0;

  for(auto &D: m_set.ByType(dipoletype::initial)){
    form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D, sqrt(m_s)/2);
  }
  // if(!m_hidephotons){
      for(auto &D: m_set.FF()){
        form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D, sqrt(D.Sprime())/2); 
        form-= CoulombSubtraction(D);
      }
    // }
  form += FormFactorSumIF();
  form += FormFactorSumDecay();
  return form;
}

// The decay-stage term of the pole expansion: one (W, charged lepton) dipole
// per W, each evaluated at half that W's own invariant mass. Empty in the flat
// scheme, where there is no decay stage.
//
// The cutoff is sqrt(D.Sprime())/2 -- the maximum photon energy in the decaying
// W's rest frame -- which is what makes this reproduce eq (12) of
// hep-ph/0302065 rather than an arbitrary fraction of it. The production stage
// above uses sqrt(s)/2 for the same reason: each stage's soft integral runs to
// its own kinematic limit, and the two stages are separately charge conserving,
// so each is separately free of the photon mass.
double Define_Dipoles::CoulombSubtraction(YFS::Dipole &D){
  // eq (9) of hep-ph/9606429. The YFS virtual B for a pair of massive charged
  // particles contains the Coulomb singularity: as their relative velocity
  // beta -> 0 it diverges like pi/(4 beta), so exp(2a ReB) is unbounded near
  // threshold. That regime is not something the exponentiation can describe --
  // it is resummed separately, here by YFS::Coulomb using the width-regulated
  // Fadin-Khoze-Martin-Stirling form -- so the singular piece is removed from
  // the exponent to avoid both the divergence and the double counting:
  //
  //   B'(beta) = B(beta) - (pi/4beta) theta(beta_t - beta)
  //
  // The exponent this sum builds is 2a ReB + 2a Btilde, so the amount to
  // remove is 2a * pi/(4 beta) = a pi / (2 beta).
  //
  // beta_t = 0.382 (YFS: WW_BETAT) is the paper's transition velocity, fixed by
  // requiring the O(a^2) Coulomb term (pi a/beta)^2/12 to stay below 0.03%
  // above it. Above beta_t the Coulomb series is an ordinary part of the
  // perturbative expansion and needs no special treatment.
  //
  // Only for a pair of massive charged vector bosons, i.e. only in the pole
  // expansion: in the flat scheme the final-state pair is light fermions, this
  // regime never arises, and subtracting would remove something real.
  if (m_wwscheme != wwscheme::pole) return 0.;
  if (!D.GetFlav(0).IsVector() || !D.GetFlav(1).IsVector()) return 0.;
  const double sp(D.Sprime());
  if (sp <= 0.) return 0.;
  const double m1s(sqr(D.GetMass(0))), m2s(sqr(D.GetMass(1)));
  const double arg(1. - 2.*(m1s + m2s)/sp + sqr((m1s - m2s)/sp));
  if (arg <= 0.) return 0.;
  const double beta(sqrt(arg));
  if (beta >= m_betatWW || beta <= 0.) return 0.;
  return m_alpha*M_PI/(2.*beta);
}


double Define_Dipoles::FormFactorSumDecay(){
  double form = 0;
  for(auto &D: m_set.Decay()){
    // Leg 0 is the resonance, leg 1 its charged daughter -- BuildPole adds
    // them in that order.
    const double M(D.GetMass(0)), ml(D.GetMass(1));
    // Soft cutoff: the largest photon energy this decay can produce in the
    // resonance rest frame, (M^2 - ml^2)/2M. The whole soft integral then sits
    // in the exponent, matching what the production stage does with sqrt(s)/2.
    if (M <= ml) continue;
    const double ks((M*M - ml*ml)/(2.*M));
    form += fabs(D.m_QiQj)*p_yfsFormFact->BVR_decay(M, ml, ks);
  }
  return form;
}

void Define_Dipoles::DumpDipoles(){
  static bool done_pole(false), done_flat(false);
  bool &done(m_pole_active ? done_pole : done_flat);
  if (done) return;
  done = true;

  auto line = [this](YFS::Dipole &D, const std::string &type, double omega) {
    double y(0.);
    if (type == "DEC") {
      const double M(D.GetMass(0)), ml(D.GetMass(1));
      if (M > ml) y = fabs(D.m_QiQj)*p_yfsFormFact->BVR_decay(M, ml, omega);
    } else if (type == "IF") {
      y = m_ifisub == 1
              ? D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega()) : 0.;
    } else {
      y = D.ChargeNorm()*p_yfsFormFact->BVR_full(D, omega);
    }
    msg_Info() << "  " << std::setw(3) << type
               << "  (" << D.GetFlav(0) << "," << D.GetFlav(1) << ")"
               << "  legs=(" << D.Left() << "," << D.Right() << ")"
               << "  QiQj=" << D.m_QiQj
               << "  ChargeNorm=" << D.ChargeNorm()
               << "  decay=" << (D.IsDecayAllowed() ? "yes" : "no ")
               << "  radiates=" << (type == "FF" ? (D.IsResonance() ? "yes" : "no ") : "n/a")
               << "  omega=" << omega << " GeV"
               << "  Y=" << y
               << std::endl;
  };

  msg_Info() << "YFS dipole set (YFS: Dump_Dipoles), " << m_set.size()
             << " dipoles"
             << (m_wwscheme == wwscheme::pole
                     ? (m_pole_active ? " [pole applied]" : " [fell back to flat]")
                     : "")
             << ":" << std::endl;
  msg_Info() << "  scheme: " << m_wwscheme
             << (m_wwscheme == wwscheme::pole && m_wwonshell ? " (W's on shell)" : "")
             << std::endl;
  for (auto &D : m_set.ByType(dipoletype::initial)) line(D, "II", sqrt(m_s)/2.);
  for (auto &D : m_set.FF())                        line(D, "FF", sqrt(D.Sprime())/2.);
  for (auto &D : m_set.Decay()) {
    const double M(D.GetMass(0)), ml(D.GetMass(1));
    line(D, "DEC", M > ml ? (M*M - ml*ml)/(2.*M) : 0.);
    msg_Info() << "        (M=" << M << " GeV, m_l=" << ml << " GeV)"
               << std::endl;
  }
  // IF dipoles only contribute to the form factor when IFI_Sub is on, so say
  // so rather than printing a cutoff that is not being used.
  for (auto &D : m_set.IF())
    line(D, "IF", m_ifisub == 1 ? IFIOmega() : 0.);
  if (m_set.FF().empty())
    msg_Info() << "  (no final-state dipoles: fewer than two charged out-legs)"
               << std::endl;
  msg_Info() << "  sum Y = " << FormFactorSum()
             << "   form factor = " << FormFactor() << std::endl;
}


double Define_Dipoles::FormFactorSumIF(){
  if(m_ifisub!=1) return 0.;
  double form = 0;
  // IFForFac = Btilda + t-channel virtual, i.e. KKMC's TForFac, and
  // ChargeNorm() = -QiQj*thetaij reproduces KKMC's +/- pattern across its
  // four TForFac calls (KKceex.cxx:315). +ChargeNorm is also what
  // TFormFactor() and every CalculateVirtualSub*() use on these same
  // dipoles, so at NLO the exponentiated IF form factor and the IF virtual
  // subtracted from the one-loop ME now carry the same sign.
  //
  // omega is sqrt(s)/2 -- one cutoff shared by all four pairs, as in KKMC,
  // and the same maximal choice the II term above makes. It was
  // sqrt(D.Sprime())/2, which for an initial-final pair is
  // (s/2)(1 -+ beta cos theta) and so is ANGLE-DEPENDENT: a soft cutoff odd
  // in cos(theta) manufactures A_FB by construction. That, not the
  // interference, was supplying essentially all of Sherpa's asymmetry.
  // Measured standalone at 0.7 GeV, |cos| < 0.55 (YFS/Tools/IFI_Budget.C):
  //
  //   assembly                              A_FB(mu+)
  //   -ChargeNorm, omega = sqrt(Sprime)/2    -0.01614   <- was live here
  //   +ChargeNorm, omega = sqrt(Sprime)/2    +0.01614
  //   +ChargeNorm, omega = sqrt(s)/2         -0.00056   <- this line
  //   KKMC CEEX2 / BabaYaga measured         -0.0103 / -0.0110
  //   Sherpa measured, before this change    -0.01605
  //
  // The -0.01614 reproduces the measured -0.01605, i.e. the whole of the old
  // asymmetry came from this one line, and it was a cutoff artifact whose
  // sign had been flipped to make it land near KKMC.
  //
  // A_FB is exactly linear in log(omega) and crosses KKMC's -0.0103 at
  // omega ~ 56 MeV, which is not a scale in the problem: no common cutoff is
  // a prediction. The remaining gap to KKMC is the real interference, which
  // KKMC gets from summing over photon partitions in the CEEX amplitude and
  // which Sherpa can only get on the emission side -- see RealIFWeight() and
  // the IFI_Real switch.
  for(auto &D: m_set.IF()){
    form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega());
  }
  return form;
}


double Define_Dipoles::IFIOmega() const {
  // Soft cutoff for the IF form factor. It must be ONE scale for all four
  // pairs, and it must be the boundary above which the interference is carried
  // by explicit photons instead.
  //
  // With IFI_Real off nothing carries it, so the only consistent choice is the
  // kinematic maximum: the exponent then holds the whole soft integral and the
  // result is cutoff-independent by construction, matching what the II term
  // does with sqrt(s)/2.
  //
  // With IFI_Real on, RealIFWeight() supplies the region above the generation
  // cutoff, so the exponent must stop there or the two double-count.
  // Honoured whenever set, independently of IFI_Real. The two used to be
  // welded together, so every run changed the exponent AND the real
  // restoration at once and no measurement could separate them. Now:
  //
  //   IFI_Real  IFI_Omega     what is being tested
  //   0         unset         baseline: whole soft integral in the exponent
  //   0         <small>       exponent lowered, nothing restoring it
  //   1         sqrt(s)/2     restoration only, exponent left inclusive
  //   1         unset         both (exponent at m_Emin + restoration)
  //
  // Only the last is a physics configuration; the middle two are for bisecting
  // which half moves sigma.
  if (m_ifiomega > 0.) return m_ifiomega;
  return sqrt(m_s)/2.;
}


double Define_Dipoles::RealIFWeight(const ATOOLS::Vec4D_Vector &photons) {
  // The interference the IF form factor no longer carries once IFIOmega() has
  // been lowered to the generation cutoff, supplied per generated photon as the
  // ratio of the COHERENT radiation function to the crude one it was generated
  // from:
  //
  //     W_IF = prod_i  CalculateRealSub(k_i) / CalculateRealSubEEX(k_i)
  //
  // Written that way the cancellation against beta_1 is exact, not approximate.
  // For one photon the weight Sherpa builds is
  //
  //     subb * [ W_IF*Born + tot ],   tot = (R*flux - subloc*Born)/subb
  //
  // with subb = CalculateRealSubEEX the crude and subloc = CalculateRealSub the
  // coherent -alpha/(4pi^2) J^2 (NLO_Base.C:524,568). Substituting
  // W_IF = subloc/subb:
  //
  //     subb*(subloc/subb)*Born + R*flux - subloc*Born  =  R*flux
  //
  // identically - the exact real matrix element, for any photon, hard or soft.
  // No soft approximation enters, which is why nothing here may be clamped: any
  // r != subloc/subb leaves a residue (subloc - subb*r)*Born, largest precisely
  // where a clamp would fire (hard, wide-angle photons, where |S_IF| rivals the
  // diagonal). That residue is the m_mumu / E_gamma distortion.
  //
  // This is also why the ratio is taken as subloc/subb rather than the
  // algebraically equal 1 + S_IF/S_crude. The three functions do NOT use
  // matching conventions: CalculateRealSub takes GetMomenta and only resonant
  // FF dipoles, CalculateRealSubEEX takes GetBornMomenta over all of them, and
  // CalculateRealSubIF takes GetMomenta. So S_crude + S_IF = S_coh holds on
  // paper but not in code, and only the explicit ratio cancels beta_1 exactly.
  //
  // The ISR and FSR generation cutoffs differ, but FSR::YFS_FORM()'s Piatek
  // term m_DelYFS puts the FSR bookkeeping back onto m_Emin, which is what
  // IFIOmega() defaults to - the same single-Emin arrangement KKMC uses for
  // Yisr, Yfsr and Yint. Still scan IFI_Omega and check the answer is flat
  // before believing a number: that is the observable statement that the
  // exponent and these photons are cancelling each other.
  if (!m_ifireal || m_set.IF().empty()) return 1.;
  const double omega = IFIOmega();
  double w(1.);
  for (const auto &k : photons) {
    if (IsZero(k.E())) continue;
    // Only photons ABOVE the cutoff. Everything below omega is already held by
    // the IF form factor, so reweighting it here counts that region twice.
    //
    // ISR photons cannot trip this - their generation cutoff is
    // (sqrt(s)/2)*m_isrcut = IR_CUTOFF/2, exactly IFIOmega()'s default - but
    // FSR photons can: FSR_CUT defaults to 1e-2*IR_CUTOFF, so the FSR
    // generation reaches about two decades lower.
    if (k.E() <= omega) continue;
    const double crude = CalculateRealSubEEX(k);
    const double coh   = CalculateRealSub(k);
    if (IsZero(crude) || IsBad(crude) || IsBad(coh)) continue;
    double r = coh/crude;
    // Off by default (IFI_RClip <= 0). The exact cancellation above says no
    // clamp is justified: r is allowed to be negative or large, because the
    // radiation function including interference genuinely is, and beta_1
    // cancels it term by term. Kept only as a diagnostic - if a run needs it to
    // stay stable, something upstream is wrong and the count below says so.
    if (m_ifi_rclip > 0.) {
      if (r < m_ifi_rclip)         { r = m_ifi_rclip;    ++m_ifi_clipped; }
      else if (r > 1./m_ifi_rclip) { r = 1./m_ifi_rclip; ++m_ifi_clipped; }
    }
    w *= r;
  }
  if (IsBad(w)) return 1.;
  return w;
}

double Define_Dipoles::FormFactor(){
  double form = FormFactorSum();
  if(FixedOrder()==fixed_order::nlo){
    return 1.+form;
  }
  return exp(form);
}


double Define_Dipoles::TFormFactor(){
  double form = 0;
  for(auto &D: m_set.ByType(dipoletype::initial)){
    form+= D.ChargeNorm()*p_yfsFormFact->R1(D);
  }
    for(auto &D: m_set.FF()){
      form += D.ChargeNorm()*p_yfsFormFact->R1(D);
  }
  if(m_ifisub==1){
    // IF dipoles use IFForFac here too, NOT R1. An initial-final pair is
    // t-channel-like by construction, so its form factor does not depend on the
    // TChannel setting -- that flag is about how the II/FF (s-channel) dipoles
    // are treated. Routing IF through R1 made the interference term change
    // whenever TChannel changed, which is why the IF treatment appeared entangled
    // with a switch that has nothing to do with it. Both FormFactorSum() and
    // TFormFactor() now give the IF dipoles the same object.
    for(auto &D: m_set.IF()){
      form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega());
    }
  }
  if(FixedOrder()==fixed_order::nlo){
    return 1.+form;
  }
  return exp(form);
}

double Define_Dipoles::CalculateVirtualSubTchannel(){
   // YFSij = 2.d0*B0ij - B0ii - B0jj
   //   .         + 4.d0 * mi2 * C0singular(mi2,phmass)
   //   .         + 4.d0 * mj2 * C0singular(mj2,phmass)
   //   .         + 8.d0*pi_pj * C0ij
  if(m_dim_reg) return CalculateVirtualSubTchannelEps();
  double sub(0);
  // Vec4D_Vector pvirt;
  // std::vector<double> z,th;
  // pvirt.push_back((*m_set.II()).GetNewMomenta(0));
  // pvirt.push_back((*m_set.II()).GetNewMomenta(1));
  // pvirt.push_back(m_set.FF()[0].GetBornMomenta(0));
  // pvirt.push_back(m_set.FF()[0].GetBornMomenta(1));
  // th.push_back(1);
  // th.push_back(1);
  // th.push_back(-1);
  // th.push_back(-1);
  // z.push_back((*m_set.II()).m_Qi);
  // z.push_back((*m_set.II()).m_Qj);
  // z.push_back(m_set.FF()[0].m_Qi);
  // z.push_back(m_set.FF()[0].m_Qj);
  // // double m1 = pvirt[0].Mass();
  // // double m2 = pvirt[1].Mass();
  // // double m3 = pvirt[2].Mass();
  // // double m4 = pvirt[3].Mass();
  // for(size_t i = 0; i < pvirt.size(); ++i)
  // {
  //   for(size_t j=i; j<pvirt.size(); ++j ){
  //     double etaij = z[i]*z[j]*th[i]*th[j];
  //     Complex YFSij = 0.;
  //     double mi = pvirt[i].Mass();
  //     double mj = pvirt[j].Mass();
  //     double s = (pvirt[i]-pvirt[j]).Abs2();

  //     Complex bii = B0(0,mi*mi,mi*mi);
  //     Complex bjj = B0(0,mj*mj,mj*mj);
  //     Complex bij = B0(0,mi*mi,mj*mj);
  //     Complex cij = C0(mi*mi,(th[i]*pvirt[i]+th[j]*pvirt[j]).Abs2(),mj*mj,
  //                      0.,mi*mi,mj*mj);
  //     Complex cii = C0(mi*mi, 0., mi*mi, 0.0, mi*mi, mi*mi);
  //     // YFSij = 8*pvirt[i]*pvirt[j]*cij;
  //     if(i==j){
  //       YFSij = bii-4.*mi*mi*0.5*log(m_photonMass*m_photonMass/mi/mi);
  //     }
  //     else{
  //       YFSij = 2.*pvirt[i]*pvirt[j]*cij+0.5*bij;
  //       }
  //     if(IsBad(YFSij)){
  //       msg_Error()<<"YFS Virtual Sub is NaN"<<endl
  //                  <<"bii = "<<bii<<endl
  //                  <<"bij = "<<bij<<endl
  //                  <<"cii = "<<cii<<endl
  //                  <<"cij = "<<cij<<endl;
  //     }
  //     sub+=m_alpi*etaij*YFSij.real();
  // //     // PRINT_VAR(etaij*YFSij);
  //   }
  // }
  // PRINT_VAR(count);
  for (auto &D : m_set.ByType(dipoletype::initial)){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_set.FF()){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  for (auto &D : m_set.IF()){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtT(D,sqrt(m_s) / 2.);
  }
  // clearcache();
  return sub;
}

double Define_Dipoles::CalculateVirtualSubTchannelEps() {
  DivArrD sub(0);
  DivArrD massph(0,-1.,0,0,0,0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
    #endif
  }
  for (auto &D : m_set.FF()) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
  }

  for (auto &D : m_set.IF()){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtTEps(D,sqrt(m_s) / 2.);
    #ifdef USING__LOOPTOOLS
      //ii term
      Vec4D p1 = D.GetBornMomenta(0);
      Vec4D p2 = D.GetBornMomenta(1);
      double m1 = D.GetMass(0);
      double m2 = D.GetMass(1);
      double irloop = p_yfsFormFact->p_virt->IRscale();
      double epsloop = p_yfsFormFact->p_virt->Eps_Scheme_Factor({p1,p1});
      DivArrD c0 = (-massph-log(4.*M_PI*sqr(irloop)/m1/m2/epsloop));
    #endif
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateRealVirtualSub(const Vec4D & k) {
  double sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetNewMomenta(0), D.GetNewMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_set.FF()) {
    sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D.GetOldMomenta(0), D.GetOldMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);

  }

  for (auto &D : m_set.IF()){
    // change to + for IFI terms
    // Note Born momenta are redifined
    // for IFI terms.
    sub += D.m_QiQj*p_yfsFormFact->BVV_full(D.GetBornMomenta(0), D.GetBornMomenta(1), m_photonMass, sqrt(m_s) / 2., 3);
  }
  return sub;
}


double Define_Dipoles::CalculateEEX(){
  double eex=0;
  for (auto &D: m_set.ByType(dipoletype::initial)){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_set.FF()){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  for (auto &D: m_set.IF()){
    D.SetRealOnly(m_real_only);
    eex += D.EEX(m_betaorder);
  }
  return eex;
}

double Define_Dipoles::CalculateEEXVirtual(){
  return CalculateEEXVirtual(m_betaorder);
}

// Same, at an EXPLICIT order rather than the runcard's BETA. Differencing two
// orders isolates one term of the EEX virtual series
// (Dipole::VirtualEEX: 0.5*gamma at order 1, +0.125*gamma^2 at order 2), which
// is how the approximate double-virtual is built in
// YFS_Handler::GenerateWeight when no exact VV provider exists.
double Define_Dipoles::CalculateEEXVirtual(int betaorder){
  double vint{1.}, vfin{1};
  for (auto &D: m_set.ByType(dipoletype::initial)){
    vint*=1+D.VirtualEEX(betaorder);
  }
  for (auto &D: m_set.FF()){
    vfin*=1+D.VirtualEEX(betaorder);
  }
  return vint*vfin;
}

double Define_Dipoles::EEXRealVirtual(const Vec4D &k){
  double eex = 0;
  for(auto &D: m_set.ByType(dipoletype::initial)){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  for(auto &D: m_set.FF()){
    D.m_betaorder = 2;
    eex += D.Beta1(k)/D.Eikonal(k);
  }
  // for(auto &D: m_set.IF()){
  //   D.m_betaorder = 2;
  //   eex += D.Beta1(k)/D.Eikonal(k);
  // }
  if(IsNan(eex)) return 0;
  return eex;
}

double Define_Dipoles::CalculateRealSubEEX(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  for (auto &D : m_set.FF()) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  // for (auto &D : m_set.IF()) {
  //   sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  // }

  return sub;
}


void Define_Dipoles::CleanInParticles() {
  m_chargedinparticles.clear();
  m_neutralinparticles.clear();
  m_massInC.clear();
  m_massInN.clear();
}

void Define_Dipoles::CleanOutParticles() {
  m_chargedoutparticles.clear();
  m_neutraloutparticles.clear();
  m_massOutC.clear();
  m_massOutN.clear();
}

void Define_Dipoles::CleanUp() {
}

double Define_Dipoles::CalculateFlux(const Vec4D &k){
  if(!HasISR()) return 1;
  double sq, sx;
  double flux = 1;
  dipoletype::code fluxtype;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  if(HasISR()&&HasFSR()){
    fluxtype = dipoletype::initial;
  }
  else if(!HasFSR()){
    fluxtype = dipoletype::initial;
  }
  else if(!HasISR()){
    fluxtype = dipoletype::final;
  }
  else{
    msg_Error()<<"Unknown dipole type in "<<METHOD<<std::endl;
  }
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_set.ByType(dipoletype::initial)) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetBornMomenta(0)+D.GetBornMomenta(1);
      sq = (QX).Abs2(); 
      sx = (QX-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_set.FF()) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_set.FF().size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, dipoletype::code &fluxtype){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  if(m_noflux==1) return 1;
  if(fluxtype==dipoletype::initial){
    for (auto &D : m_set.ByType(dipoletype::initial)) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2(); 
      sx = (Q-k).Abs2();
      flux = (sx/sq);
      return flux;
    }

  }
  if(fluxtype==dipoletype::final){
    flux=0;
    for (auto &D : m_set.FF()) {
      Q  = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux += (sq/sx);
      // flux = Propagator(sx)/Propagator(sq);
    }
    return flux/m_set.FF().size();
  }
  return flux;
}

double Define_Dipoles::CalculateFlux(const Vec4D &k, const Vec4D &kk){
  double sq, sx;
  double flux = 1;
  Vec4D Q,QX;
  dipoletype::code fluxtype1, fluxtype2;
  if(m_noflux==1) return 1;
  fluxtype1 = WhichResonant(k);
  fluxtype2 = WhichResonant(kk);
  if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::initial){
    for (auto &D : m_set.ByType(dipoletype::initial)) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k-kk).Abs2();
      flux = sx/sq;
      return flux;
    }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::final){
    for (auto &D : m_set.FF()) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k+kk).Abs2();
      flux = sq/sx;
    }
    return flux;
  }
  else if(fluxtype1==dipoletype::initial && fluxtype2==dipoletype::final){
    for (auto &D : m_set.ByType(dipoletype::initial)) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-k).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_set.FF()) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+kk).Abs2();
      flux*= sq/sx;
   }
  }
  else if(fluxtype1==dipoletype::final && fluxtype2==dipoletype::initial){
    for (auto &D : m_set.ByType(dipoletype::initial)) {
      QX = D.GetNewMomenta(0)+D.GetNewMomenta(1);
      Q =  D.GetMomenta(0)+D.GetMomenta(1);
      sq = Q.Abs2();
      sx = (Q-kk).Abs2();
      flux = sx/sq;
    }
    for (auto &D : m_set.FF()) {
      Q = D.GetBornMomenta(0)+D.GetBornMomenta(1);
      QX = D.GetMomenta(0)+D.GetMomenta(1);
      sq = (Q).Abs2();
      sx = (Q+k).Abs2();
      flux*= sq/sx;
    }
  }
  else{
    msg_Error()<<"Unknown flux type in "<<METHOD<<" with \nfluxtype1 = "<<fluxtype1<<"\n and fluxtype2 = "<< fluxtype2 <<std::endl;
  }
  return flux;
}


double Define_Dipoles::Propagator(const double &s, int width){
  Flavour fl;
  Complex Prop = Complex(0.,0.);///Complex(s,0.0);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      fl = v->in[0];
      if(IsZero(fl.Mass())) continue;
      Prop += Complex(1.,0.)/Complex(s-sqr(fl.Mass()),fl.Width()*fl.Mass());
    }
  }
  return ((Prop*conj(Prop)).real());
}

double Define_Dipoles::ResonanceWidthDistance(YFS::Dipole &D) {
  // |m_ij - M| / Gamma, minimised over the resonances of the process. Returns a
  // large sentinel if the process has no usable resonance, so that a pair stays
  // selectable and the ordering simply falls back to the order of the pairs.
  const double mass_d((D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass());
  double mdist(std::numeric_limits<double>::max());
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if (IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width())) continue;
      mdist = std::min(mdist, std::abs(mass_d - v->in[0].Mass()) / v->in[0].Width());
    }
  }
  return mdist;
}



void Define_Dipoles::IsResonant(YFS::Dipole &D) {
double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      // if(D.IsDecayAllowed()){
      //   D.SetResonance(true);
      //   continue;
      // }
      if(D.m_QiQj==1 || !D.IsDecayAllowed()){
        D.SetResonance(false);
        continue;
        }   
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<m_resonace_max) {
        D.SetResonance(true);
        return;
      }
      else D.SetResonance(false);
    }
    D.SetResonance(false);
  }
  if(D.IsDecayAllowed())  D.SetResonance(true);
  else  D.SetResonance(false);
}

bool Define_Dipoles::CheckResonant(YFS::Dipole &D) {
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist;
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
      if(mdist<5) {
        return true;
      }
    }
  }
  return false;
}

bool Define_Dipoles::IsResonant() {
  for(auto &D: m_set.FF()){
    double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
    double mdist;
    for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
      for (auto *v : it->second) {
        mdist = abs(mass_d - v->in[0].Mass()) / v->in[0].Width();
        if(mdist<5) {
          return true;
        }
      }
    }
  }
  return false;
}

bool Define_Dipoles::CheckResonant(){
  bool isres = false;
  for(auto &D: m_set.ByType(dipoletype::initial)){
    if(CheckResonant(D)) isres=true;
  }
  for(auto &D: m_set.FF()){
    if(CheckResonant(D)) isres=true;
  }
  return isres;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D, const Vec4D &k){
  double mass_d = (D.GetBornMomenta(0) + D.GetBornMomenta(1)-k).Mass();
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_d- mass_i);
      if(mcheck < mdist) {
        mdist = mcheck;
        D.SetResonaceFlavour(v->in[0]);
      }
    }
  }
  return mdist;
}

double Define_Dipoles::ResonantDist(YFS::Dipole &D){
  double mass_i = (D.GetBornMomenta(0) + D.GetBornMomenta(1)).Mass();
  double mdist(100000000);
  double mcheck(100000000);
  for (auto it = m_proc_restab_map.begin(); it != m_proc_restab_map.end(); ++it) {
    for (auto *v : it->second) {
      if(IsZero(v->in[0].Mass()) || IsZero(v->in[0].Width()) ) continue;
      mcheck = abs(mass_i - v->in[0].Mass());
      if(mcheck < mdist) mdist = mcheck;
    }
  }
  return mdist;
}

dipoletype::code Define_Dipoles::WhichResonant(const Vec4D &k){
  if(!HasFSR()) return dipoletype::initial;
  if(!HasISR()) return dipoletype::final;
  double mdistisr(10000), mdistfsr(100000),mdistifi(100000);
  double mindis(10000);
  dipoletype::code min(dipoletype::initial);
  for(auto &D: m_set.ByType(dipoletype::initial)){
    mdistisr = ResonantDist(D,k);
    mindis = mdistisr;
    min = dipoletype::initial;
  }
  for(auto &D: m_set.FF()){
    mdistfsr = ResonantDist(D,k);  
    if(mdistfsr < mdistisr){
      min = dipoletype::final;
    }
  }
  return min;
}

void Define_Dipoles::generate_pairings(std::vector<std::vector<int>>& pairings, std::vector<int>& curr_pairing, std::vector<int>& available_nums) {
  if (available_nums.empty()) {
    pairings.push_back(curr_pairing);
    return;
  }
  int curr_num = available_nums[0];
  available_nums.erase(available_nums.begin());
  for(size_t i = 0; i < available_nums.size(); i++) {
    int next_num = available_nums[i];
    available_nums.erase(available_nums.begin() + i);
    curr_pairing.push_back(curr_num);
    curr_pairing.push_back(next_num);
    generate_pairings(pairings, curr_pairing, available_nums);
    curr_pairing.pop_back();
    curr_pairing.pop_back();
    available_nums.insert(available_nums.begin() + i, next_num);
  }
  available_nums.insert(available_nums.begin(), curr_num);
}

std::ostream& Define_Dipoles::operator<<(std::ostream &out) {
  out << "N_in = " << m_in << "\n m_out = " << m_out <<
      "Number of Charged incoming particles = " << m_chargedinparticles.size() << std::endl <<
      "Number of Charged outgoing particles = " << m_chargedoutparticles.size() << std::endl <<
      "Number of Neutral incoming particles = " << m_neutralinparticles.size() << std::endl <<
      "Number of Neutral outgoing particles = " << m_neutraloutparticles.size() << std::endl;
  return out;
}

std::vector<ParticleInfo> Define_Dipoles::ExtractChargedParticles(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta,
    size_t start_index,
    size_t end_index,
    bool is_initial_state) const {
    
    std::vector<ParticleInfo> charged_particles;
    
    for (size_t i = start_index; i < end_index; ++i) {
        if (flavors[i].IntCharge() != 0) {
            charged_particles.emplace_back(
                flavors[i], momenta[i], born_momenta[i], i, is_initial_state);
        }
    }
    
    return charged_particles;
}

std::vector<ParticleInfo> Define_Dipoles::ExtractInitialStateCharged(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta) const {
    
    constexpr size_t INITIAL_STATE_PARTICLES = 2;
    return ExtractChargedParticles(flavors, momenta, born_momenta, 
                                   0, INITIAL_STATE_PARTICLES, true);
}

std::vector<ParticleInfo> Define_Dipoles::ExtractFinalStateCharged(
    const ATOOLS::Flavour_Vector& flavors,
    const ATOOLS::Vec4D_Vector& momenta,
    const ATOOLS::Vec4D_Vector& born_momenta) const {
    
    constexpr size_t INITIAL_STATE_PARTICLES = 2;
    return ExtractChargedParticles(flavors, momenta, born_momenta, 
                                   INITIAL_STATE_PARTICLES, flavors.size(), false);
}

Dipole Define_Dipoles::CreateDipole(
    const ParticleInfo& particle1,
    const ParticleInfo& particle2,
    dipoletype::code type) const {
    
    ATOOLS::Flavour_Vector flavors = {particle1.flavor, particle2.flavor};
    ATOOLS::Vec4D_Vector momenta = {particle1.momentum, particle2.momentum};
    ATOOLS::Vec4D_Vector born_momenta = {particle1.born_momentum, particle2.born_momentum};
    
    Dipole dipole(flavors, momenta, born_momenta, type, m_alpha);
    dipole.SetFlavLab(particle1.index, particle2.index);

    return dipole;
}
