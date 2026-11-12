#include <sstream>
#include "YFS/Main/Define_Dipoles.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
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
  p_yfsFormFact = std::make_unique<YFS::YFS_Form_Factor>();
}

void Define_Dipoles::SetS(const double s) {
  YFS_Base::SetS(s);
  p_yfsFormFact->SetS(s);
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
}


void Define_Dipoles::MakeDipolesII(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born) {
  m_N_born_Gamma=1;
  for(auto f: fl) if(f.IsPhoton()) m_N_born_Gamma+=1;
  if(!HasISR()) return;
  if (m_flav_label_fl != fl) {
    m_flav_label.clear();
    for(size_t i(0); i<fl.size(); ++i) m_flav_label[fl[i]] = i;
    m_flav_label_fl = fl;
  }
  m_bornmomenta = born;
  m_set.BuildInitial(fl, mom, born, m_alpha);
  // BuildInitial leaves no dipole if the initial state carries no charge.
  if (m_set.HasII()) { m_g = m_set.II().m_gamma; m_gp = m_set.II().m_gammap; }
  // the ISR crude's region: photons kept above IR_CUTOFF/2 in P (ISR.C)
  if (m_set.HasII())
    m_set.II().SetGenRegion(SoftRegion::Lab(0.5*sqrt(m_s)*m_isrcut, "P"));
}


void Define_Dipoles::MakeDipolesIF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const mom, ATOOLS::Vec4D_Vector const born) {
  if (m_mode == yfsmode::fsr) return;
  if (!HasFSR()) return;
  if (mom.size() != fl.size())
    THROW(fatal_error, "Incorrect dipole size in YFS for initial-final dipoles");
  m_set.BuildIF(fl, mom, born, m_alpha);
}

void Define_Dipoles::MakeDipoles(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom, ATOOLS::Vec4D_Vector const &born ) {
  if ((mom.size() != fl.size()) || (born.size() != fl.size()))
    THROW(fatal_error, "Incorrect dipole size in YFS for final-state dipoles");
  if (m_flav_label_fl != fl) {
    m_flav_label.clear();
    for(size_t i(0); i<fl.size(); ++i) m_flav_label[fl[i]] = i;
    m_flav_label_fl = fl;
  }
  m_bornmomenta = born;
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
  const bool built(m_set.BuildPole(fl, mom, born, m_alpha, m_wwonshell,
                                   m_resonace_max, m_wwpoleemission));
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

  // BuildPole adds the pair as (W-, W+): leg 0 carries w.lm/w.nm, leg 1 w.lp/w.np.

  const Vec4D oldw[2] = {D.GetOldMomenta(0), D.GetOldMomenta(1)};
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



void Define_Dipoles::Dipole_FF(ATOOLS::Flavour_Vector const &fl, ATOOLS::Vec4D_Vector const &mom) {
  if (fl.size() != mom.size()) {
    THROW(fatal_error, "Inconsistent flavour vector for Dipole_FF Momenta");
  }
}

double Define_Dipoles::CalculateRealSub(const Vec4D &k) {
  double sub(0);
  Vec4D eik{0.,0.,0.,0.}, eikI{0.,0.,0.,0.}, eikF{0.,0.,0.,0.};
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    for(size_t i = 0; i < 2; ++i)
    {
       Vec4D p = D.GetMomenta(i);
      eik  += D.m_Q[i]*p/(p*k);
      eikI += D.m_Q[i]*p/(p*k);
    }
  }
  for (auto &D : m_set.FF()) {
     if(!D.IsResonance()) continue;
    for(size_t i = 0; i < 2; ++i)
    {
      Vec4D p = D.GetMomenta(i);
      eik  += -D.m_Q[i]*p/(p*k);
      eikF += -D.m_Q[i]*p/(p*k);
    }
  }
  if (m_ifi) sub = -m_alpha / (4 * M_PI * M_PI)*eik*eik;
  else sub = -m_alpha / (4 * M_PI * M_PI)*(eikI*eikI + eikF*eikF);
  /*
    Identical photons in the Born (YFS: REAL_BORN_PHOTON_SYM). The eikonal is
    returned UNDIVIDED: the n-extra-photon real matrix elements carry the
    symmetry factor 1/(N+n)! against the Born's 1/N!, and NLO_Base multiplies
    each of them by (N+n)!/N! (BornPhotonSym), so every subtraction here is
    the plain soft factor.
  */
  static const bool bornphsym(ATOOLS::Settings::GetMainSettings()["YFS"]
                              ["REAL_BORN_PHOTON_SYM"].SetDefault(true).Get<bool>());
  if (bornphsym) return sub;
  return sub/(m_N_born_Gamma!=0?m_N_born_Gamma:1.0);
}

double Define_Dipoles::CalculateRealSubIF(const Vec4D &k) {
  double sub(0);
  if (!m_ifi) return sub;   // YFS: IFI false: no initial-final eikonal
  for (auto &D : m_set.IF()){
    if(m_massless_sub) sub += D.EikonalMassless(k, D.GetMomenta(0), D.GetMomenta(1));
    else sub +=  D.Eikonal(k, D.GetMomenta(0), D.GetMomenta(1));
  }
  return sub;
}


double Define_Dipoles::CalculateVirtualSub() {
  double sub(0);
  if(m_dim_reg==1) return CalculateVirtualSubEps();
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  for (auto &D : m_set.FF()) {
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full(D, m_photonMass, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }

  for (auto &D : m_set.IF()){
    sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  }
  // Pole-scheme decay stage: its exponent is in the YFS weight, so its virtual B
  // is subtracted like every other dipole's.
  if (DecayVirtualSubtraction())
    for (auto &D : m_set.Decay())
      sub += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  return sub;
}

/*
  YFS: IF_VIRTUAL_TCHANNEL (default true): the dim.-reg. virtual subtraction of an
  initial-final dipole uses its t-channel B (YFS_Form_Factor::BVirtT_eps), the B the
  exponent carries (IFForFac). false restores BVV_full_eps, the s-channel form, whose
  finite part differs by a C-odd amount (see BVirtT_eps).
*/
DivArrD Define_Dipoles::IFVirtualBEps(YFS::Dipole &D) {
  static const bool tch(ATOOLS::Settings::GetMainSettings()["YFS"]["IF_VIRTUAL_TCHANNEL"]
                        .SetDefault(true).Get<bool>());
  return tch ? p_yfsFormFact->BVirtT_eps(D) : p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
}

void Define_Dipoles::VirtualSubParts(double &ii, double &ff, double &ifi) {
  ii = ff = ifi = 0.;
  if (m_dim_reg == 1) {   // CalculateVirtualSubEps, finite parts
    for (auto &D : m_set.ByType(dipoletype::initial))
      ii += (D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3)).Finite();
    for (auto &D : m_set.FF())
      if (!D.IsFinite())
        ff += (D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3)).Finite();
    for (auto &D : m_set.IF())
      if (!D.IsFinite())
        ifi += (D.ChargeNorm()*IFVirtualBEps(D)).Finite();
    return;
  }
  for (auto &D : m_set.ByType(dipoletype::initial))
    ii += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  for (auto &D : m_set.FF())
    ff += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
  for (auto &D : m_set.IF())
    ifi += D.ChargeNorm()*p_yfsFormFact->BVirtGeneral(D, sqrt(m_s) / 2.);
}

/*
  YFS: WW_DECAY_VIRTUAL_SUB (default 0). Whether the pole expansion's decay
  dipoles enter the virtual subtraction. They should: the exponent carries
  FormFactorSumDecay, so the loop's decay-stage collinear logs have to be
  subtracted with a B of the same object. But the once-per-run
  CheckDecayFormFactor shows the closed-form exponent 
*/
bool Define_Dipoles::DecayVirtualSubtraction() const {
  static const bool on(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["WW_DECAY_VIRTUAL_SUB"].Get<bool>());
  return on;
}

/*
  YFS: RV_MODE 1 with RR_MODE 1 (NLO_Base::CalculateRealVirtualRemainder):
  sum over the same dipoles and charge factors as CalculateVirtualSub of the
  real-soft integral B-tilde(K) (YFS_Form_Factor::BVR_full, mode 0), with the
  photon mass of the subtraction's regularisation: m_photonMass, or with
  Dim_Reg the finite-part convention of BVV_full_eps, lambda^2 = 4 pi mu^2
  / eps_scheme, so that B(virtual, finite) + B-tilde is the lambda-free Y.
*/
double Define_Dipoles::RealSoftSum(double K) {
  auto one = [&](YFS::Dipole &D) {
    Vec4D p1, p2;
    D.LegsBeforeRadiation(p1, p2);
    double lam(m_photonMass);
    if (m_dim_reg == 1) {
      PHASIC::Virtual_ME2_Base *v(p_yfsFormFact->p_virt);
      if (v == nullptr) return std::numeric_limits<double>::quiet_NaN();
      const double eps(v->Eps_Scheme_Factor({p1, p2}));
      lam = sqrt(4.*M_PI*sqr(v->IRscale())/eps);
    }
    return p_yfsFormFact->BVR_full(p1*p2, p1.E(), p2.E(), D.GetMass(0), D.GetMass(1),
                                   K, lam, 0); };
  double sum(0.);
  for (auto &D : m_set.ByType(dipoletype::initial)) sum += D.ChargeNorm()*one(D);
  for (auto &D : m_set.FF()) {
    if (m_dim_reg == 1 && D.IsFinite()) continue;
    if (m_mode == yfsmode::fsr) sum += -D.m_QiQj*one(D);
    else sum += D.ChargeNorm()*one(D);
  }
  for (auto &D : m_set.IF()) {
    if (m_dim_reg == 1 && D.IsFinite()) continue;
    sum += D.ChargeNorm()*one(D);
  }
  return sum;
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
    sub += D.ChargeNorm()*IFVirtualBEps(D);
    if(IsBad(sub.Finite())) {
      msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
      // THROW(fatal_error, "YFS Subtraction fails");
    }
  }
  /*
    The pole expansion's decay stage (W incoming, charged daughter outgoing).
    Its exponent FormFactorSumDecay is in the YFS weight, so its virtual B
    must be subtracted from the loop like the production dipoles' - without
    it the decay-stage collinear logs log(M_W/m_l) sit in exp(Y) AND in the
    provider's virtual
  */
  for (auto &D : m_set.Decay()) {
    CheckDecayFormFactor(D);
    if(!DecayVirtualSubtraction() || D.IsFinite()) continue;
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    if(IsBad(sub.Finite())) msg_Error()<<"YFS subtraction is Nan For dipole:"<<D<<std::endl;
  }
  m_virtSub=sub;
  return sub.Finite();
}

double Define_Dipoles::CalculateRealVirtualSubEps(const Vec4D &k) {
  // Virtual YFS B-hat in dim-reg at the reduced kinematics. Uses the same BVV
  // function as CalculateVirtualSubEps so beta_1^1 vanishes in the soft limit.
  DivArrD sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }
  for (auto &D : m_set.FF()) {
    if(D.IsFinite()) continue;
    if(m_mode==yfsmode::fsr) sub += -D.m_QiQj*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
    else sub += D.ChargeNorm()*p_yfsFormFact->BVV_full_eps(D, sqrt(m_s) / 2., 3);
  }

  for (auto &D : m_set.IF()){
    // Born momenta are redefined for IFI terms.
    if(D.IsFinite()) continue;
    sub += D.ChargeNorm()*IFVirtualBEps(D);
  }
  m_virtSub=sub;
  return sub.Finite();
}



double Define_Dipoles::FormFactorSumFinalLegs(){
  double form = 0;
  for(auto &D: m_set.FF()){
    form+= D.ChargeNorm()*p_yfsFormFact->BVR_full(D, sqrt(D.Sprime())/2);
    form-= CoulombSubtraction(D);
  }
  form += FormFactorSumIF();
  form += FormFactorSumDecay();
  return form;
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

// Decay-stage term of the pole expansion (one (W, lepton) dipole per W, empty in
// the flat scheme). The cutoff is the maximum photon energy in the decaying W's
// rest frame, sqrt(D.Sprime())/2, reproducing eq (12) of hep-ph/0302065.
double Define_Dipoles::CoulombSubtraction(YFS::Dipole &D){
  // eq (9) of hep-ph/9606429. The virtual B diverges as pi/(4 beta) at threshold;
  // that Coulomb piece is resummed separately (YFS::Coulomb), so it is removed
  // from the exponent: B'(beta) = B(beta) - (pi/4beta) theta(beta_t - beta).
  // The exponent is 2a ReB + 2a Btilde, so the amount removed is a pi/(2 beta).
  // beta_t = 0.382 (YFS: WW_BETAT). Only for massive charged vector pairs, i.e.
  // the pole expansion.
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


void Define_Dipoles::CheckDecayFormFactor(YFS::Dipole &D) {
  if (m_decaycheckdone || m_dim_reg != 1) return;
  m_decaycheckdone = true;
  // Legs in the W rest frame, where Kmax = (M^2 - m_l^2)/2M is defined.
  const double M(D.GetMass(0)), ml(D.GetMass(1));
  if (M <= ml) return;
  const double ks((M*M - ml*ml)/(2.*M));
  Vec4D_Vector legs{D.GetBornMomenta(0), D.GetBornMomenta(1)};
  Poincare toW(legs[0]);
  for (Vec4D &p : legs) toW.Boost(p);
  YFS::Dipole Dw(D.GetFlavors(), legs, legs, D.Type(), m_alpha);
  Dw.SetMass(0, M); Dw.SetMass(1, ml);
  const DivArrD bv(p_yfsFormFact->BVV_full_eps(Dw, ks, 3));
  const DivArrD br(p_yfsFormFact->BVR_full_eps(Dw, ks, 0));
  const double ymach(D.ChargeNorm()*(bv + br).Finite());
  const double yclosed(fabs(D.m_QiQj)*p_yfsFormFact->BVR_decay(M, ml, ks));
  msg_Info()<<"YFS pole scheme, decay-stage exponent check ("<<D.GetFlav(0)
            <<" -> "<<D.GetFlav(1)<<", M = "<<M<<" GeV): closed form Y = "
            <<yclosed<<", dim-reg virtual + real B = "<<ymach
            <<" (1/eps pole of the sum "<<(bv + br).GetIR()
            <<"); the two must agree for the O(alpha) subtraction to match "
            <<"the exponent."<<std::endl;
}

double Define_Dipoles::FormFactorSumDecay(){
  double form = 0;
  for(auto &D: m_set.Decay()){
    // Leg 0 is the resonance, leg 1 its charged daughter.
    const double M(D.GetMass(0)), ml(D.GetMass(1));
    // Soft cutoff: largest photon energy in the resonance rest frame, (M^2 - ml^2)/2M.
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
      y = m_ifisub
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
    // exponent region and generated region
    msg_Info() << "        exponent: "
               << (D.HasExpRegion() ? D.ExpRegion().Describe() : std::string("-"))
               << "   generated: "
               << (D.HasGenRegion() ? D.GenRegion().Describe() : std::string("-"))
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
  // IF dipoles contribute to the form factor only when IFI_Sub is on.
  for (auto &D : m_set.IF())
    line(D, "IF", m_ifisub ? IFIOmega() : 0.);
  if (m_set.FF().empty())
    msg_Info() << "  (no final-state dipoles: fewer than two charged out-legs)"
               << std::endl;
  msg_Info() << "  sum Y = " << FormFactorSum()
             << "   form factor = " << FormFactor() << std::endl;
}


double Define_Dipoles::FormFactorSumIF(){
  if (!m_ifisub) return 0.;
  double form = 0;

  for(auto &D: m_set.IF()){
    form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, IFIOmega());
  }
  return form;
}


/*
  The exponent the sampled photon multiplicities belong to. FormFactorSum()
  integrates each initial dipole's eikonal up to sqrt(s)/2 and each final one
  up to sqrt(s')/2 in its own rest frame, while every generated (visible)
  photon has energy above emin = IR_CUTOFF/2 in the beams' rest frame: the
  event weight exp(FormFactorSum()) times the sampled photons stands for
  exp(Y(emin)) prod S~(k) with Y(emin) = FormFactorSum() + this shift. Only the
  real B-tilde moves; the virtual is the same on both sides and drops out.
  Read-only on the dipoles (the low-level BVR_full, not the one that records
  the region on the dipole).
*/
double Define_Dipoles::FormFactorShiftToRegion(double emin){
  const SoftRegion gen(SoftRegion::Lab(emin, "P-gen"));
  auto real = [this](const YFS::Dipole &D, const SoftRegion &reg) {
    Vec4D p1, p2;
    D.LegsBeforeRadiation(p1, p2);
    return p_yfsFormFact->BVR_full(reg, p1, p2, p1.Mass(), p2.Mass(), m_photonMass, 0);
  };
  double shift(0.);
  for (auto &D : m_set.ByType(dipoletype::initial))
    shift += D.ChargeNorm()*(real(D, gen)
             - real(D, p_yfsFormFact->ExponentRegion(D, sqrt(m_s)/2)));
  for (auto &D : m_set.FF())
    shift += D.ChargeNorm()*(real(D, gen)
             - real(D, p_yfsFormFact->ExponentRegion(D, sqrt(D.Sprime())/2)));
  return shift;
}


double Define_Dipoles::FormFactorShiftToRegionFinal(double emin){
  const SoftRegion gen(SoftRegion::Lab(emin, "P-gen"));
  auto real = [this](const YFS::Dipole &D, const SoftRegion &reg) {
    Vec4D p1, p2;
    D.LegsBeforeRadiation(p1, p2);
    return p_yfsFormFact->BVR_full(reg, p1, p2, p1.Mass(), p2.Mass(), m_photonMass, 0);
  };
  double shift(0.);
  for (auto &D : m_set.FF())
    shift += D.ChargeNorm()*(real(D, gen)
             - real(D, p_yfsFormFact->ExponentRegion(D, sqrt(D.Sprime())/2)));
  return shift;
}


double Define_Dipoles::FormFactorSumIF(double omega){
  if (!m_ifisub) return 0.;
  double form = 0;
  for(auto &D: m_set.IF())
    form += D.ChargeNorm()*p_yfsFormFact->IFForFac(D, omega);
  return form;
}


double Define_Dipoles::IFIOmega() const {
  // Soft cutoff for the IF form factor: one scale for all four pairs, above which
  // the interference is carried by explicit photons.
  if (m_ifiomega > 0.) return m_ifiomega;
  return sqrt(m_s)/2.;
}


double Define_Dipoles::RealIFWeight(const ATOOLS::Vec4D_Vector &photons) {
  // Interference not carried by the IF form factor once IFIOmega() is lowered to
  // the generation cutoff: per generated photon, coherent over crude radiation function.
  
  if (!m_ifireal || m_set.IF().empty()) return 1.;
  const double omega = IFIOmega();
  double w(1.);
  for (const auto &k : photons) {
    if (IsZero(k.E())) continue;
    // Only photons above omega; below it the IF form factor already holds them.
    if (k.E() <= omega) continue;
    const double crude = CalculateRealSubEEX(k);
    const double coh   = CalculateRealSub(k);
    if (IsZero(crude) || IsBad(crude) || IsBad(coh)) continue;
    double r = coh/crude;
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
  return FormFactor(FormFactorSum());
}

double Define_Dipoles::FormFactor(double form){
  return exp(form);
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

double Define_Dipoles::CalculateRealSubEEX(const Vec4D &k) {
  double sub(0);
  for (auto &D : m_set.ByType(dipoletype::initial)) {
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  static const bool rad(ATOOLS::Settings::GetMainSettings()["YFS"]
                        ["REAL_CRUDE_RADIATING"].SetDefault(false).Get<bool>());
  for (auto &D : m_set.FF()) {
    if (rad && !D.IsResonance()) continue;
    sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  }
  // for (auto &D : m_set.IF()) {
  //   sub += D.Eikonal(k, D.GetBornMomenta(0), D.GetBornMomenta(1));
  // }

  return sub;
}


double Define_Dipoles::CalculateFlux(const Vec4D &k){
  double sq, sx;
  double flux = 1;
  dipoletype::code fluxtype;
  Vec4D Q,QX;
  if (m_noflux) return 1;
  /*
    MODE: FSR. The Jacobian of the (n+1)-body point against the crude,
    (Q - k)^2/Q^2 with Q the pre-emission final state, is the same whether
    ISR is on or off; only its bookkeeping differs (with ISR on it is read
    off the II dipole's reduced momenta below). Returning 1 here left the
    single-FSR-photon Born+real weight 1/(1 - x) above ACRAIC in FSR-only mode
  */
  if(!HasISR()){
    for (size_t i(2); i < m_bornmomenta.size(); ++i) Q += m_bornmomenta[i];
    sq = Q.Abs2();
    if (sq <= 0.) return 1;
    return (Q-k).Abs2()/sq;
  }
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
  if (m_noflux) return 1;
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
  if (m_noflux) return 1;
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
  // |m_ij - M| / Gamma, minimised over the process resonances; a large sentinel
  // if there is none, so the pair order is used.
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

