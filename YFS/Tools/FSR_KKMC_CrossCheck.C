// Standalone comparison harness: drives YFS::FSR/Dipole with fixed,
// hand-picked multi-photon configurations (bypassing NPhotons()/
// GeneratePhotonMomentum()'s RNG sampling) so the result can be compared
// number-for-number against KKMC::KKarFin::YFSfin on identical input. The
// matching KKMC-side driver lives at
// Test/SherpaCompare/kkmc_fsr_crosscheck.cxx in the KKMC repo. Currently
// injects 2 photons (see k0s/rns/phis below) - trivially extendable to more
// by adding entries to those three vectors on both sides.
//
// This does NOT modify FSR.C/Dipole.C - it only calls their already-public
// methods/members in the same sequence YFS_Handler::CalculateFSR does,
// substituting a fixed photon for the randomly-sampled one.
//
// Build/run: see build_FSR_KKMC_CrossCheck.sh in this directory. Requires a
// prior `ninja libYFS.dylib` (or full build) so the .dylib it links against
// is up to date.
//
// History: this harness found and confirmed two bugs (2026-08-04):
// 1. GenerateAngles()'s m_kkmcAngles==1 branch: m_fbarvec used
//    (1+beta1*beta2)/(beta1+beta2) instead of (1+beta1*beta2)/2,
//    overweighting every FSR photon's mass-weight denominator by a factor
//    of beta relative to KKMC. Fixed in FSR.C; the "rigorous check" block
//    below calls the real GenerateAngles() (not a hand-derived mirror)
//    specifically so this stays a live regression check.
// 2. YFS_FORM()'s m_BtiXcru/m_BtiQcru defaulted to BVR_full (exact Btilda)
//    via FSR_CRU=0, but KKMC's Piatek() always uses the crude/truncated
//    Btildc for this ghost-momentum role (no equivalent flag on their
//    side). Confirmed via an 8-9 sig-fig match of BtiXcru/BtiQcru/DelVol
//    against KKMC's real KKbvir::Btildc (see the matching KKMC-side
//    driver's Piatek()-equivalent block). Fixed by changing FSR_CRU's
//    default to 1 in FSR.C.
#include "ATOOLS/Org/Settings.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Terminator_Objects.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Library_Loader.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/KF_Table.H"
#include "PDF/Main/PDF_Base.H"
#include "YFS/Main/Dipole.H"
#include "YFS/Main/FSR.H"
#include "MODEL/Main/Model_Base.H"
#include <iostream>
#include <iomanip>

using namespace ATOOLS;
using namespace YFS;

// Minimal stand-in for a real MODEL::Model_Base: YFS_Base::RegisterSettings()
// queries s_model->ScalarConstant("alpha_QED") unconditionally on one branch
// even when USE_MODEL_ALPHA is off, so *some* model has to be present.
class FakeModel : public MODEL::Model_Base {
public:
  FakeModel(double alpha) : MODEL::Model_Base(true) {
    (*p_constants)["alpha_QED"] = alpha;
  }
  void InitVertices() override {}
  void ParticleInit() override {}
  bool ModelInit() override { return true; }
};

// One phase-space point, driven end to end in the same sequence
// YFS_Handler::CalculateFSR uses. Factored out of main() so the harness can run
// it twice: once for a Born dipole at rest (s' = s) and once for a dipole that
// has recoiled against ISR (s' < s, and moving in the CMS).
//
// The second point is the one that was missing. What it adds is NOT the
// CMS-vs-Q-frame split - the FSR photons already recoil the dipole, so the leg
// energies entering m_btilStar differ from the Eq1/Eq2 entering m_btil even at
// rest (measured 0.988/1.027 there, 0.947/1.071 with ISR on). What it adds is
// that m_dip_sp stops being the nominal m_s: 0.49 -> 0.343 at v = 0.3. Any
// place that reaches for m_s where it wants the dipole's own s', or the other
// way round, is invisible until those two differ, and until now they never did
// in this harness.
static int RunPoint(const char *label, const Vec4D &pA, const Vec4D &pB,
                    double mpi, double alpha) {
  std::cout << "\n===============================================================\n"
            << "=== " << label << "\n"
            << "===============================================================\n";
  std::cout << std::setprecision(10);
  std::cout << "leg A     = " << pA << "\n"
            << "leg B     = " << pB << "\n"
            << "dipole s' = " << (pA+pB).Abs2()
            << "    nominal s = " << sqr(rpa->gen.Ecms()) << "\n";
  Flavour flPip((kf_code)211, false);
  Flavour flPim((kf_code)211, true);

  Flavour_Vector fl{flPip, flPim};
  Vec4D_Vector   mom{pA, pB};
  Vec4D_Vector   born{pA, pB};

  Dipole dip(fl, mom, born, dipoletype::final, alpha);
  // Make sure the dipole uses the physical pion mass/charge regardless of
  // whether the KF table has kf=211 registered in this bare-bones setup.
  dip.m_masses  = {mpi, mpi};
  dip.m_charges = {1., -1.};
  dip.SetFlavLab(0, 1);
  dip.SetResonance(true);
  dip.SetBorn(1.0);
  dip.BoostToQFM(false); // Q already at rest -> no-op boost, but computes CalculateGamma()

  FSR fsr;
  fsr.Reset();
  fsr.SetV(0.0);
  if (!fsr.Initialize(dip)) { std::cerr << "FSR::Initialize failed\n"; return 1; }
  // Initialize() reads masses off the dipole via p_dipole->m_masses, but does
  // its own bookkeeping too - force-correct in case Flavour::Mass() silently
  // returned 0 for an unregistered kf in this minimal setup.
  fsr.m_mass = {mpi, mpi};

  // ---- rigorous weight check: call the REAL (rebuilt) GenerateAngles() ----
  // Unlike the hand-derived mirror below (needed for the kinematics test,
  // which requires a chosen, reproducible photon direction), this calls the
  // actual library code so a future FSR.C change can't silently go
  // unreflected here. ran is seeded deterministically (seed 1234 above), so
  // this is reproducible run to run, just not hand-picked.
  fsr.GenerateAngles();
  double real_costh = fsr.m_c, real_sinth = fsr.m_st, real_fbar = fsr.m_fbarvec[0];
  {
    double am2r = sqr(fsr.m_mass[0]+fsr.m_mass[1]) / fsr.m_dip_sp;
    // del1*del2 here, NOT 1-costh^2 - they only coincide when am2==0.
    double del1r = 1. - fsr.m_beta1*real_costh;
    double del2r = 1. + fsr.m_beta2*real_costh;
    double dist0_for_real_costh = 1./(del1r*del2r) * (1.-am2r/2.);
    std::cout << "--- rigorous check: actual GenerateAngles() output vs KKMC dist0 for that costh ---\n";
    std::cout << "real costh = " << real_costh << "  real sinth = " << real_sinth << "\n";
    std::cout << "real m_fbarvec[0] (from rebuilt FSR.C) = " << real_fbar << "\n";
    std::cout << "KKMC dist0 for this SAME costh          = " << dist0_for_real_costh << "\n";
    std::cout << "ratio (should be 1.0 after the fix)     = " << real_fbar/dist0_for_real_costh << "\n";
  }
  fsr.m_cos.clear(); fsr.m_sin.clear(); fsr.m_fbarvec.clear(); fsr.m_MassWls.clear();

  // ---- inject TWO fixed photons instead of sampling them ----
  // Each photon's costh/sinth is derived from its OWN crude-sampling random
  // draw "rn_i" using the SAME formula as FSR::GenerateAngles()'s
  // m_kkmcAngles==1 branch / KKMC's AngBre() (no 50/50 symmetrization swap),
  // so the matching KKMC driver can derive bit-identical costh/sinth/dist0
  // per photon from the same inputs. am2_crude/beta_crude/eps_crude are
  // dipole-level (shared across all photons in the event), matching
  // GenerateAngles() being called once per photon but m_beta1/m_beta2 being
  // fixed for the whole dipole.
  std::vector<double> k0s = {0.3, 0.25, 0.08};
  std::vector<double> rns = {0.02, 0.97, 0.5};
  std::vector<double> phis = {0.0, 3.0, 5.5};
  double am2_crude = sqr(fsr.m_mass[0]+fsr.m_mass[1]) / fsr.m_dip_sp;
  double beta_crude = sqrt(1.-am2_crude);
  double eps_crude  = am2_crude/(1.+beta_crude);

  std::vector<double> cosths, sinths, dist0_crudes, fbarvecs;
  Vec4D_Vector photons_in;
  Vec4D photonSum_in;
  for (size_t i = 0; i < k0s.size(); ++i) {
    double rn = rns[i], k0 = k0s[i], phi = phis[i];
    double del1_c = (2.-eps_crude)*pow(eps_crude/(2.-eps_crude), rn);
    double del2_c = 2.-del1_c;
    double costh  = (del2_c-del1_c)/(2.*beta_crude);
    double sinth  = sqrt(del1_c*del2_c - am2_crude*costh*costh);
    double dist0_crude = 1./(del1_c*del2_c)*(1.-am2_crude/2.); // == KKMC's dis0[i]
    // Matches GenerateAngles()'s m_kkmcAngles==1 branch: after sampling
    // costh/sinth from the crude density above, it recomputes del1/del2
    // with the EXACT per-particle betas for the weight denominator - for
    // equal masses m_beta1==m_beta2==beta_crude so this is a no-op here.
    double del1_s = 1. - fsr.m_beta1*costh;
    double del2_s = 1. + fsr.m_beta2*costh;
    // Fixed formula: was erroneously /(m_beta1+m_beta2), now /2.
    double fbarvec = 1./(del1_s*del2_s) * (1.+fsr.m_beta1*fsr.m_beta2)/2.;

    cosths.push_back(costh);
    sinths.push_back(sinth);
    dist0_crudes.push_back(dist0_crude);
    fbarvecs.push_back(fbarvec);
    Vec4D photon(k0, k0*sinth*cos(phi), k0*sinth*sin(phi), k0*costh);
    photons_in.push_back(photon);
    photonSum_in += photon;
  }

  fsr.m_n   = (int)k0s.size();
  fsr.m_c   = cosths.back();
  fsr.m_st  = sinths.back();
  fsr.m_phi = phis.back();
  fsr.m_cos = cosths;
  fsr.m_sin = sinths;
  fsr.m_fbarvec = fbarvecs;
  fsr.m_MassWls = std::vector<double>(k0s.size(), 1.0);
  fsr.m_k0  = k0s;
  fsr.m_photons        = photons_in;
  fsr.m_photonSum      = photonSum_in;
  fsr.m_photonspreboost = fsr.m_photons;

  // ---- replicate MakeFSR()'s post-sampling logic (see YFS/Main/FSR.C) ----
  fsr.m_dist1.clear(); fsr.m_dist2.clear(); fsr.m_del1.clear(); fsr.m_del2.clear();
  fsr.m_cut = 1; fsr.m_wt2 = 1.0; fsr.m_yy = 1.0; fsr.m_xfact = 1.0;
  fsr.m_sQ = fsr.m_dip_sp;
  double smin = sqr(fsr.m_mass[0]+fsr.m_mass[1]);
  if (fsr.m_photons.size()==0) {
    fsr.m_sprim = fsr.m_dip_sp;
    fsr.m_sQ = fsr.m_sprim;
  } else {
    if (fsr.m_photonSum.E() >= 1) { std::cerr << "photon too hard\n"; return 1; }
    fsr.RescalePhotons();
    fsr.m_sQ = fsr.m_dip_sp * fsr.m_yy;
    fsr.m_sX = fsr.m_sQ*(1.+fsr.m_photonSum[0]+0.25*fsr.m_photonSum*fsr.m_photonSum);
    if (fsr.m_sQ < smin || fsr.m_sX < smin) { std::cerr << "below threshold\n"; return 1; }
  }
  fsr.m_u = 1. - fsr.m_sprim/fsr.m_dip_sp;
  fsr.MakePair(sqrt(fsr.m_sprim), fsr.m_dipole[0], fsr.m_dipole[1]);
  fsr.m_px = fsr.m_dipole[0]+fsr.m_dipole[1]+fsr.m_photonSum;
  fsr.m_Q  = fsr.m_dipole[0]+fsr.m_dipole[1];

  // ---- ghost/crude momenta (FSR::MakeFSR(), see YFS/Main/FSR.C ~304-317) ----
  // r1,r2 use rescaled masses masc1,masc2 (matching sQ, not the physical
  // masses) - this is what the crude Btilde (m_BtiXcru/m_BtiQcru inside
  // YFS_FORM()) is meant to be evaluated against.
  double masc1 = fsr.m_mass[0] * sqrt(fsr.m_sQ / fsr.m_dip_sp);
  double masc2 = fsr.m_mass[1] * sqrt(fsr.m_sQ / fsr.m_dip_sp);
  double eta1_ghost, eta2_ghost;
  fsr.MakePair(sqrt(fsr.m_sprim), fsr.m_r1, fsr.m_r2, masc1, masc2, eta1_ghost, eta2_ghost);
  fsr.CalculateBetaBar();

  for (int i=0;i<2;++i) dip.SetMomentum(i, fsr.m_dipole[i]);
  dip.AddToGhosts(fsr.m_r1);
  dip.AddToGhosts(fsr.m_r2);

  // ---- mass-weight (FSR::F(), see YFS/Main/FSR.C) ----
  fsr.F();

  // ---- dipole recoil boost (same call YFS_Handler::CalculateFSR makes) ----
  dip.SetNPhoton((int)fsr.m_photons.size()); // silences Boost()'s "Wrong Photon multiplicity" check
  Vec4D_Vector photons = fsr.m_photons;
  dip.AddPhotonsToDipole(photons);
  dip.Boost();

  // ---- sanity check: does GetGhost(0)/(1) actually return the r1/r2 ghosts
  // just added? Looked like a real bug at first read of AddToGhosts()
  // (it only ever push_back()s, and the constructor's mom loop seeds
  // m_ghost with 2 entries) - but Dipole's constructor has a later
  // `if (ty==final) m_ghost.clear();` that empties it again for FF dipoles
  // specifically, so AddToGhosts()'s two calls land at indices 0/1 as
  // intended. Confirmed correct below (kept as a standing regression check
  // rather than removed, since it's exactly the kind of thing that's easy
  // to break again). ----
  std::cout << "--- ghost/crude-momentum indexing sanity check ---\n";
  std::cout << "fsr.m_r1 (ghost, as computed, pre-boost) = " << fsr.m_r1 << "\n";
  std::cout << "fsr.m_r2 (ghost, as computed, pre-boost) = " << fsr.m_r2 << "\n";
  std::cout << "dip.GetGhost(0) (what YFS_FORM() actually reads) = " << dip.GetGhost(0) << "\n";
  std::cout << "dip.GetGhost(1) (what YFS_FORM() actually reads) = " << dip.GetGhost(1) << "\n";
  std::cout << "dip.GetNewMomenta(0) (physical dipole leg 0, for comparison) = " << dip.GetNewMomenta(0) << "\n";

  // ---- form factor / hide-photon reweighting (FSR::YFS_FORM(), see
  // YFS/Main/FSR.C) - this is the piece that mirrors KKMC's Piatek(). ----
  bool yfsform_ok = fsr.YFS_FORM();
  std::cout << "--- YFS_FORM() (Piatek-equivalent), ok=" << yfsform_ok << " ---\n";
  std::cout << std::setprecision(10);
  std::cout << "m_bvrA      = " << fsr.m_bvrA << "\n";
  std::cout << "m_BtiXcru   = " << fsr.m_BtiXcru << "\n";
  std::cout << "m_BtiQcru   = " << fsr.m_BtiQcru << "\n";
  std::cout << "m_volmc     = " << fsr.m_volmc << "\n";
  std::cout << "m_btilStar  = " << fsr.m_btilStar << "\n";
  std::cout << "m_btil      = " << fsr.m_btil << "\n";
  std::cout << "m_DelYFS    = " << fsr.m_DelYFS << "\n";
  std::cout << "m_delvol    = " << fsr.m_delvol << "\n";
  std::cout << "m_hideW     = " << fsr.m_hideW << "\n";
  std::cout << "m_YFS_IR    = " << fsr.m_YFS_IR << "\n";
  std::cout << "m_Emin      = " << fsr.m_Emin << "  m_EminQ = " << fsr.m_EminQ << "\n";
  std::cout << "--- EminQ breakdown ---\n";
  std::cout << "fsr.m_fsrcut = " << fsr.m_fsrcut << "\n";
  std::cout << "fsr.m_sQ     = " << fsr.m_sQ << "\n";
  std::cout << "fsr.m_Q      = " << fsr.m_Q << "\n";
  std::cout << "fsr.m_photonSum (post-YFS_FORM recompute) = " << fsr.m_photonSum << "\n";
  {
    double Eqq_check = 0.5*sqrt(fsr.m_sQ);
    double QdotK = fsr.m_Q*fsr.m_photonSum;
    double Delta_check = fsr.m_fsrcut*(1.+2.*QdotK/fsr.m_sQ);
    std::cout << "Eqq_check = " << Eqq_check << "  Q.photonSum = " << QdotK
               << "  Delta_check = " << Delta_check
               << "  EminQ_check = " << Eqq_check*Delta_check << "\n";
  }

  // ---- partial-variant test: does X or Q alone need to switch cru/full,
  // rather than both together (FSR_CRU)? YFS_IR_local/m_DelYFS/m_volmc don't
  // depend on the X/Q formula choice (verified algebraically: the A4/A4sng
  // mass-singular terms are identical for BtiXcru and BtiQcru since both use
  // the same r1/r2, so they cancel in delvol regardless of which formula is
  // used - only the "-1" in BVR_full's log term survives the subtraction).
  // Back out YFS_IR_local from the already-computed, verified m_hideW rather
  // than re-deriving the formula by hand a second time. ----
  {
    double YFS_IR_local = log(fsr.m_hideW) - fsr.m_DelYFS - fsr.m_volmc + fsr.m_delvol;
    double r1r2 = fsr.m_r1 * fsr.m_r2;
    double BtiXcru_full = fsr.p_fsrFormFact->BVR_full(r1r2, fsr.m_r1[0], fsr.m_r2[0], fsr.m_r1.Mass(), fsr.m_r2.Mass(), fsr.m_Emin, fsr.m_photonMass, 0);
    double BtiXcru_cru  = fsr.p_fsrFormFact->BVR_cru (r1r2, fsr.m_r1[0], fsr.m_r2[0], fsr.m_r1.Mass(), fsr.m_r2.Mass(), fsr.m_Emin);
    double BtiQcru_full = fsr.p_fsrFormFact->BVR_full(r1r2, fsr.m_r1[0], fsr.m_r2[0], fsr.m_r1.Mass(), fsr.m_r2.Mass(), fsr.m_EminQ, fsr.m_photonMass, 0);
    double BtiQcru_cru  = fsr.p_fsrFormFact->BVR_cru (r1r2, fsr.m_r1[0], fsr.m_r2[0], fsr.m_r1.Mass(), fsr.m_r2.Mass(), fsr.m_EminQ);
    // KKMC's Piatek() uses EQQ=0.5*sqrt(svarQ) for BOTH particles' energy in
    // the Q-scale crude evaluation (Btildc(...,EQQ,EQQ,...,EminQ,...)), not
    // the ghosts' own post-boost energies - a second, separate discrepancy
    // from the full/cru formula choice tested above.
    double Eqq_val = 0.5*sqrt(fsr.m_sQ);
    double BtiQcru_cru_Eqq  = fsr.p_fsrFormFact->BVR_cru (r1r2, Eqq_val, Eqq_val, fsr.m_r1.Mass(), fsr.m_r2.Mass(), fsr.m_EminQ);
    std::cout << "--- partial-variant test (same YFS_IR/DelYFS/volmc, varying only the X/Q formula choice) ---\n";
    struct Variant { const char* name; double delvol; };
    Variant variants[5] = {
      {"X=full,Q=full (current default, FSR_CRU=0)", BtiXcru_full - BtiQcru_full},
      {"X=cru, Q=cru  (FSR_CRU=1)                  ", BtiXcru_cru  - BtiQcru_cru},
      {"X=full,Q=cru                               ", BtiXcru_full - BtiQcru_cru},
      {"X=cru, Q=full                               ", BtiXcru_cru  - BtiQcru_full},
      {"X=cru, Q=cru with EQQ,EQQ (matches KKMC)   ", BtiXcru_cru  - BtiQcru_cru_Eqq},
    };
    for (auto &v : variants) {
      double hideW = exp(YFS_IR_local + fsr.m_DelYFS + fsr.m_volmc - v.delvol);
      std::cout << v.name << ": delvol=" << v.delvol
                << "  hideW=" << hideW
                << "  ratio_to_current_default=" << hideW/fsr.m_hideW << "\n";
    }
  }

  std::cout << std::setprecision(10);
  std::cout << "dip_sp        = " << fsr.m_dip_sp << "  sqrt = " << sqrt(fsr.m_dip_sp) << "\n";
  std::cout << "sprim         = " << fsr.m_sprim   << "  sqrt = " << sqrt(fsr.m_sprim) << "\n";
  std::cout << "yy            = " << fsr.m_yy << "\n";
  std::cout << "wt2 (dilatation jacobian, FSR::RescalePhotons()) = " << fsr.m_wt2 << "\n";
  std::cout << "dipole[0] (pre-boost, XFM frame)  = " << fsr.m_dipole[0] << "\n";
  std::cout << "dipole[1] (pre-boost, XFM frame)  = " << fsr.m_dipole[1] << "\n";
  std::cout << "newmomenta[0] (final, lab frame)  = " << dip.GetNewMomenta(0) << "\n";
  std::cout << "newmomenta[1] (final, lab frame)  = " << dip.GetNewMomenta(1) << "\n";
  Vec4D_Vector out_photons = dip.GetPhotons();
  Vec4D total = dip.GetNewMomenta(0)+dip.GetNewMomenta(1);
  for (size_t i=0;i<out_photons.size();++i) {
    std::cout << "photon["<<i<<"]    (final, lab frame)  = " << out_photons[i] << "\n";
    total += out_photons[i];
  }
  std::cout << "total out = " << total << "\n";
  std::cout << "born  sum = " << (pA+pB) << "\n";

  std::cout << "--- weight (FSR::F(), KKMC_ANG-style crude, " << k0s.size() << " photons) ---\n";
  double massW_total = 1.0;
  for (size_t i=0;i<k0s.size();++i) {
    std::cout << "photon["<<i<<"]: rn="<<rns[i]<<" k0="<<k0s[i]<<" costh="<<cosths[i]<<" sinth="<<sinths[i]<<"\n";
    std::cout << "  dist0_crude (== KKMC dis0[i]) = " << dist0_crudes[i] << "\n";
    std::cout << "  fbarvec (Sherpa's own m_fbar) = " << fbarvecs[i] << "\n";
    std::cout << "  m_dist1[i]  (Sherpa's m_f, exact) = " << fsr.m_dist1[i] << "\n";
    std::cout << "  m_dist2[i]  (Sherpa's m_fbar used) = " << fsr.m_dist2[i] << "\n";
    std::cout << "  m_MassWls[i] (Sherpa's f/fbar) = " << fsr.m_MassWls[i] << "\n";
    massW_total *= fsr.m_MassWls[i];
  }
  std::cout << "massW_total (product over photons) = " << massW_total << "\n";
  std::cout << "eta1=" << fsr.m_eta1 << " eta2=" << fsr.m_eta2 << "\n";
  std::cout << "--- combined FSR::Weight()-equivalent total (wgt::full: massW*hideW*wt2) ---\n";
  std::cout << "total = " << massW_total * fsr.m_hideW * fsr.m_wt2 << "\n";


  // ---- the quantities ISR actually moves ----
  // Run both points and diff this block. m_Emin must NOT follow s' - it is the
  // single nominal-CMS scale the Piatek term translates everything onto, the
  // same role KKMC's m_Emin plays for Yisr, Yfsr and Yint alike - while
  // m_EminQ must. The leg-energy ratio is already off 1 at rest (the FSR
  // photons see to that) and simply grows.
  {
    const double sq  = fsr.m_sQ;
    const double Eqq = 0.5*sqrt(sq);
    const double Eq1 = (sq + sqr(fsr.m_mass[0]) - sqr(fsr.m_mass[1]))/(2.*sqrt(sq));
    const double Eq2 = (sq + sqr(fsr.m_mass[1]) - sqr(fsr.m_mass[0]))/(2.*sqrt(sq));
    const double nominal = sqrt(sqr(rpa->gen.Ecms()));
    std::cout << "--- ISR-sensitive quantities (diff this block between the two points) ---\n";
    std::cout << "leg energies entering m_btilStar (CMS)     = "
              << dip.GetNewMomenta(0)[0] << ", " << dip.GetNewMomenta(1)[0] << "\n";
    std::cout << "leg energies entering m_btil     (Q frame) = "
              << Eq1 << ", " << Eq2 << "   (Eqq = " << Eqq << ")\n";
    std::cout << "  ratio CMS/Q per leg                      = "
              << dip.GetNewMomenta(0)[0]/Eq1 << ", "
              << dip.GetNewMomenta(1)[0]/Eq2 << "\n";
    // m_Emin is the single scale the Piatek term translates onto, so it must be
    // the nominal-CMS soft cutoff and must NOT follow s'. Kept as a standing
    // check against KKMC's own definition, CMSene/2 * vvmin (KKarFin.cxx:56),
    // rather than against a value handed over from the Sherpa side - which is
    // the only way the factor sqrt(s) this once carried could have been caught.
    std::cout << "m_Emin  Sherpa  (FSR.C:89)                 = " << fsr.m_Emin << "\n";
    std::cout << "m_Emin  KKMC    CMSene/2 * vvmin           = "
              << 0.5*nominal*fsr.m_isrcut << "\n";
    std::cout << "  ratio Sherpa/KKMC (must be 1)            = "
              << fsr.m_Emin/(0.5*nominal*fsr.m_isrcut) << "\n";
    std::cout << "m_EminQ (tracks s', as it should)          = " << fsr.m_EminQ << "\n";
  }

  return 0;
}

int main(int argc, char* argv[]) {
#ifdef USING__MPI
  MPI_Init(&argc, &argv);
#endif
  // Same minimal singleton bring-up SHERPA::Sherpa's constructor does
  // (SHERPA/Main/Sherpa.C), skipping everything not needed by FSR/Dipole.
  ATOOLS::mpi = new My_MPI();
  ATOOLS::exh = new Terminator_Object_Handler();
  ATOOLS::msg = new Message();
  ATOOLS::rpa = new Run_Parameter();
  // No MODEL::Model_Base is loaded in this standalone harness, so
  // YFS_Base::RegisterSettings() must not try to pull alpha_QED from
  // s_model. Override it the same way a runcard/CLI arg would, since
  // SetDefault() can't override RegisterDefaults()'s own SetDefault(1).
  // KKMC_ANG:1 selects FSR::GenerateAngles()'s KKMC-matching branch, which is
  // what the production runcard (Sherpa.Pion.NLO.New.yaml) actually uses.
  // IR_CUTOFF matches the production runcard (Sherpa.Pion.NLO.New.yaml) -
  // the default (1e-6) gives a very different Emin/EminQ ratio, which
  // matters a lot for m_delvol's magnitude.
  // FSR_CRU is left at its (now-fixed) default of 1 deliberately - this
  // exercises the actual production default, not an override.
  char override_arg[] = "YFS: {USE_MODEL_ALPHA: 0, KKMC_ANG: 1, IR_CUTOFF: 7e-11}";
  char prog_name[] = "sherpa_fsr_driver";
  char* fake_argv[] = {prog_name, override_arg};
  Settings::InitializeMainSettings(2, fake_argv);
  ATOOLS::ran = new Random(1234);
  ATOOLS::s_loader = new Library_Loader();
  PDF::pdfdefs = new PDF::PDF_Defaults();

  try {

  // ---- fixed phase-space point: e+e- -> pi+ pi- at sqrt(s)=0.7 GeV ----
  double mpi    = 0.13957039;
  double sqrts  = 0.7;
  double pmag   = sqrt(sqr(sqrts/2.)-sqr(mpi));

  // YFS_Base::YFS_Base() sets m_s = sqr(rpa->gen.Ecms()) at construction
  // time (YFS_Base.C:32) - without this, m_s=0 propagates into m_isrcut
  // (divided by sqrt(m_s)) and m_Emin, silently NaN-ing YFS_FORM(). Must be
  // set before any YFS_Base subclass (Dipole, FSR) is constructed below.
  rpa->gen.SetEcms(sqrts);

  Vec4D pPip(sqrts/2., 0., 0.,  pmag);
  Vec4D pPim(sqrts/2., 0., 0., -pmag);

  // No Model is loaded in this standalone harness, so kf=211 isn't
  // registered in the global KF table yet - Flavour::Charge()/Mass() would
  // dereference a null Particle_Info* without this. icharge is in units of
  // e/3, so 3 -> charge +1.
  s_kftable[211] = new Particle_Info(211, mpi, 0., 3, 0, 0, "pi+", "\\pi^+");
  double alpha = 1./137.035999084;
  MODEL::s_model = new FakeModel(alpha);

  Flavour flPip((kf_code)211, false);
  Flavour flPim((kf_code)211, true);


  // 1. Dipole at rest, s' = s. What this harness has always run, and the
  //    configuration whose numbers are already matched against KKMC.
  if (RunPoint("no ISR : dipole at rest, s' = s", pPip, pPim, mpi, alpha))
    std::cerr << "no-ISR point failed\n";

  // 2. ISR on. An ISR photon of energy v*sqrt(s)/2 along +z has been emitted,
  //    so the pi+pi- system carries invariant mass s' = s(1-v) and recoils
  //    against it. This is what YFS_Handler hands to MakeDipoles(): the FF
  //    dipole is built from m_plab[2..], i.e. the final state AFTER the ISR
  //    recoil, while m_s stays nominal.
  {
    const double v  = 0.30;
    const double kz = v*sqrts/2.;
    const Vec4D  Q(sqrts - kz, 0., 0., -kz);       // pi pair recoiling against k
    const double sp = Q.Abs2();
    const double pstar = sqrt(sp/4. - sqr(mpi));
    const double ct = 0.6, st = sqrt(1.-ct*ct);
    Vec4D q1(sqrt(sp)/2.,  pstar*st, 0.,  pstar*ct);
    Vec4D q2(sqrt(sp)/2., -pstar*st, 0., -pstar*ct);
    Poincare toCMS(Q);
    toCMS.BoostBack(q1);
    toCMS.BoostBack(q2);
    if (RunPoint("ISR on : dipole recoiling, s' = (1-v) s, v = 0.30", q1, q2, mpi, alpha))
      std::cerr << "ISR-on point failed\n";
  }

  } catch (const ATOOLS::Exception& e) {
    std::cerr << "ATOOLS::Exception: " << e << std::endl;
    return 1;
  } catch (const std::exception& e) {
    std::cerr << "std::exception: " << e.what() << std::endl;
    return 1;
  }

  return 0;
}
