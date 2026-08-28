#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "YFS/NLO/NLO_Base.H"
#include "YFS/NLO/Virtual.H"
#include "YFS/NLO/VirtualVirtual.H"
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

using namespace YFS;
using namespace MODEL;
using namespace ATOOLS;
using namespace std;

double massmin = 2220;
double rcount = 1;
double sumw = 0;


// Lambda (Kaellen function) now lives once in YFS/Tools/Dipole.H.

template <typename Get>
static void HardestBetas(const YFS::Photon_Vector &ph, Get get,
                         double &h1, double &h2) {
  std::vector<const YFS::Photon *> s;
  s.reserve(ph.size());
  for (const YFS::Photon &k : ph) s.push_back(&k);
  std::sort(s.begin(), s.end(),
            [](const YFS::Photon *a, const YFS::Photon *b) {
              return a->E() > b->E();
            });
  h1 = s.empty() ? 0. : get(*s[0]);
  h2 = h1 + (s.size() < 2 ? 0. : get(*s[1]));
}

NLO_Base::NLO_Base() {
  p_yfsFormFact = new YFS::YFS_Form_Factor();
  p_nlodipoles = new YFS::Define_Dipoles();
  p_real = NULL;
  p_virt = NULL;
  p_realvirt = NULL;
  p_vv = NULL;
  m_evts = 0;
  m_recola_evts = 0;
  m_realtool = 0;
  m_realvirt = 0;
  m_looptool = 0;
  m_rrtool = 0;
  m_vvtool = 0;
  m_zeroRV = 0;
  m_zeroRR = 0;
  m_nonZeroRR = 0;
  m_zeroV = 0;
  m_nonZeroRV=0;
  m_real_hard1 = 0.;
  m_rv_hard1 = 0.;
  m_rr_hard2 = 0.;
  m_real_hard2 = 0.;
  m_rv_hard2 = 0.;
  m_zero_real_amp = 0;
  m_ceex_done = false;
  m_softRV = 0;
  m_softRR = 0;
  m_rvUnstable = 0;
  m_rvHiC = 0;
  m_rvBlowup = 0;
  m_rvBlowupRtree0 = 0;
  m_rvBlowupHiC = 0;
  m_rvBlowupSoft = 0;
  m_rvBlowupHardWide = 0;
  BookHistograms();
  if (m_check_poles == 1) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    m_histograms1d["SinglePoleCD"] = new Histogram(0, 0, 25, 25);
    m_histograms1d["SinglePoleVV"] = new Histogram(0, 0, 25, 25);
    m_histograms1d["DoublePoleVV"] = new Histogram(0, 0, 25, 25);
    m_histograms1d["OneLoopEpsLP"] = new Histogram(0, -1.5, -0.5, 50);
    m_histograms1d["OneLoopEpsYFS"] = new Histogram(0, -1.5, -0.5, 50);
    m_histograms1d["RealLoopEpsLP"] = new Histogram(0, -5, 0.0, 50);
    m_histograms1d["RealLoopEpsYFS"] = new Histogram(0, -5, 0.0, 50);
    m_histograms1d["relativediff"] = new Histogram(0, -20., -5.0, 50);
    m_histograms1d["RVSinglePoleCD"] = new Histogram(0, 0, 25, 25);
    m_histograms2d["REAL_SUB"] =
        new Histogram_2D(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
    m_histograms2d["REAL"] =
        new Histogram_2D(0, 0, sqrt(m_s) / 2., 200, 0, 2 * M_PI, 20);
  }
  if (m_rv_cancel_hist) {
    if (!ATOOLS::DirectoryExists(m_debugDIR_NLO))
      ATOOLS::MakeDir(m_debugDIR_NLO);
    // Plateau-scan diagnostics for the RV soft/cancellation cut. Each RV photon
    // that reaches the beta_1^1 (i.e. passes the coarse energy pre-filter) is
    // filled here BEFORE the cancellation guard, weighted by its contribution
    // tot. Summing the weight from the hard side down to a given bin gives
    // RV_total(cut): a plateau followed by jitter marks where physics ends and
    // roundoff begins - the safe placement for RV_CANCEL_EPS / RV_SOFT_CUT.
    // *_w = weighted by tot (the integrand); *_n = unweighted counts per bin.
    m_histograms1d["RV_tot_by_logC_w"] = new Histogram(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_logC_n"] = new Histogram(0, -16., 0., 80);
    m_histograms1d["RV_tot_by_Efrac_w"] = new Histogram(0, 0., 0.5, 100);
    m_histograms1d["RV_tot_by_Efrac_n"] = new Histogram(0, 0., 0.5, 100);
    // Matrix-element stability ratio log10(|rv| / subtraction-scale) for every
    // RV photon. ~0 when healthy (rv tracks its subtraction), large when the
    // one-loop provider is unstable. "_all" is the full sample; "_hardwide" is
    // restricted to energetic, well-separated photons (E/sqrt(s)>0.1, no charged
    // leg within ~26deg). A tail at large log10 in "_hardwide" means the
    // instability reaches hard wide-angle (cannot simply skip); if only "_all"
    // has the tail, the failures are purely soft/collinear (safe to skip).
    m_histograms1d["RV_MEstab_all"] = new Histogram(0, -2., 40., 84);
    m_histograms1d["RV_MEstab_hardwide"] = new Histogram(0, -2., 40., 84);
    // Same binning but weighted by the contribution tot, so cumulative-from-low
    // gives RV_total(cut) as a function of the RV_ME_MAX_RATIO threshold: the
    // physical RV is the plateau reached before the instability shoulder/spike
    // starts inflating the integral. Answers "how much does the surviving
    // ratio=1e2..1e4 shoulder inflate RV?" from a single run.
    m_histograms1d["RV_tot_by_MEstab_w"] = new Histogram(0, -2., 40., 84);
  }
}

NLO_Base::~NLO_Base() {
  WriteHistograms();
  // PRINT_VAR(massmin)
  if (p_yfsFormFact)
    delete p_yfsFormFact;
  if (p_nlodipoles)
    delete p_nlodipoles;
  // p_virt, p_real, p_realvirt, p_realreal and p_vv are NOT deleted here: they
  // are owned by the YFS_Process they were built for (see SetProviders), since
  // one NLO_Base is shared by every process in the run card.
  msg_Out()<<"Total zero V: "<<m_zeroV<<std::endl;
  msg_Out()<<"Total zero RV: "<<m_zeroRV<<std::endl;
  msg_Out()<<"Total zero RR: "<<m_zeroRR<<std::endl;
  msg_Out()<<"Total non-zero RR: "<<m_nonZeroRR<<std::endl;
  msg_Out()<<"Total non-zero RV: "<<m_nonZeroRV<<std::endl;
  // Job-wide totals for the production RV guards. These are plain ints printed
  // only on rank 0, so Allreduce them (collective: every rank runs this
  // destructor and reaches the call). Always on - independent of the
  // RV_CANCEL_HIST diagnostics below.
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    int gbuf[3] = {m_softRV, m_rvUnstable, m_softRR};
    mpi->Allreduce(gbuf, 3, MPI_INT, MPI_SUM);
    m_softRV = gbuf[0];
    m_rvUnstable = gbuf[1];
    m_softRR = gbuf[2];
  }
#endif
  msg_Out()<<"Total soft RV skipped: "<<m_softRV<<std::endl;
  msg_Out()<<"Total unstable-ME RV skipped (RV_ME_MAX_RATIO): "<<m_rvUnstable<<std::endl;
  msg_Out()<<"Total soft RR pairs skipped: "<<m_softRR<<std::endl;
  if (m_rv_cancel_hist) {
    // Sum the per-rank diagnostic counters across all ranks. Unlike the
    // histograms (Allreduced in MPISync), these are plain ints printed only on
    // rank 0, so without this a blow-up landing on a non-zero rank would be
    // invisible. Collective: every rank runs this destructor and m_rv_cancel_hist
    // is read identically everywhere, so all ranks reach the Allreduce together.
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      int buf[6] = {m_rvHiC,      m_rvBlowup,    m_rvBlowupRtree0,
                    m_rvBlowupHiC, m_rvBlowupSoft, m_rvBlowupHardWide};
      mpi->Allreduce(buf, 6, MPI_INT, MPI_SUM);
      m_rvHiC = buf[0];
      m_rvBlowup = buf[1];
      m_rvBlowupRtree0 = buf[2];
      m_rvBlowupHiC = buf[3];
      m_rvBlowupSoft = buf[4];
      m_rvBlowupHardWide = buf[5];
    }
#endif
    msg_Out()<<"RV photons with C>=1 (subtraction not cancelling): "<<m_rvHiC<<std::endl;
    msg_Out()<<"RV blow-ups |tot|>1e3*|Born|: "<<m_rvBlowup
             <<"  (of these: rtree==0: "<<m_rvBlowupRtree0
             <<", C>=1: "<<m_rvBlowupHiC
             <<", soft E/sqrt(s)<0.01: "<<m_rvBlowupSoft
             <<", HARD WIDE-ANGLE: "<<m_rvBlowupHardWide<<")"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> hard wide-angle fraction of blow-ups: "
               <<(100.*m_rvBlowupHardWide/m_rvBlowup)<<"%"
               <<" (if >0, instability is NOT confined to soft/collinear)"<<std::endl;
    if (m_rvBlowup>0)
      msg_Out()<<"  -> rtree==0 fraction of blow-ups: "
               <<(100.*m_rvBlowupRtree0/m_rvBlowup)<<"%"
               <<" (mechanism confirmed if ~100%)"<<std::endl;
  }
  msg_Out()<<"Total zero real amplitudes: "<<m_zero_real_amp<<std::endl;
  msg_Out()<<"Total events : "<<m_evts<<std::endl;
}

void NLO_Base::SetProviders(YFS::Virtual *virt, YFS::Real *real,
                            YFS::RealVirtual *realvirt, YFS::RealReal *realreal,
                            YFS::VirtualVirtual *vv) {
  // Non-owning: the providers belong to the YFS_Process they were built for.
  // The "has this correction" flags follow the pointers, so switching process
  // switches the available corrections consistently.
  p_virt     = virt;
  p_real     = real;
  p_realvirt = realvirt;
  p_realreal = realreal;
  p_vv       = vv;
  m_looptool = (p_virt     != NULL);
  m_realtool = (p_real     != NULL);
  m_realvirt = (p_realvirt != NULL);
  m_rrtool   = (p_realreal != NULL);
  m_vvtool   = (p_vv       != NULL);
}

void NLO_Base::Init(Flavour_Vector &flavs, Vec4D_Vector &plab,
                    Vec4D_Vector &born) {
  m_flavs = flavs;
  m_plab = plab;
  m_bornMomenta = born;
  m_borngamma.clear();
  for (int i = 0; i < flavs.size(); i++) {
    if (m_flavs[i].IsPhoton())
      m_borngamma.push_back(m_plab[i]);
  }
}

double NLO_Base::CalculateNLO() {
  m_failcut = false;
  double result{0.0};
  // PRINT_VAR(m_looptool);
  result += CalculateReal();
  result += CalculateVirtual();
  result += CalculateRealReal();
  result += CalculateRealVirtual();
  result += CalculateVV();
  if(result>1e10){
    msg_Error()<<"Perturbative corrections to large!"<<std::endl;
    return 0;
  }
  return result;
}

double NLO_Base::CalculateVirtual() {
  if (m_eex_virt) {
    // subtract born to avoid double counting
    // already present in eex!!
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  if (!m_looptool)
    return 0;
  double virt;
  double sub;
  p_dipoles->p_yfsFormFact->p_virt = p_virt->p_loop_me;
  CheckMassReg();
  if (!HasISR())
    virt = p_virt->Calc(m_plab, m_born);
  else
    virt = p_virt->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_virt->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  m_oneloop = (virt - sub * m_born / m_rescale_alpha);
  if (IsZero(virt)){
    m_zeroV++;
    return 0;
  }
  if (p_virt->p_loop_me->Mode() == 1)
    m_oneloop /= m_rescale_alpha;
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  if (m_check_poles == 1) {
    if (m_virt_sub == 0)
      sub = p_dipoles->CalculateVirtualSub();
    double p1 = p_virt->p_loop_me->ME_E1() * p_virt->m_factor;
    double yfspole = p_dipoles->Get_E1();
    int ncorrect = ::countMatchingDigits(p1, -yfspole);
    double reldiff = (p1 + yfspole) / p1;
    if (!IsEqual(p1, -yfspole, 1e-4)) {
      msg_Error() << "Poles do not cancel in YFS Virtuals" << std::endl
                  << "Correct digits =  " << ncorrect << std::endl
                  << "Relative diff =  " << reldiff
                  << std::endl
                  // <<"Process =  "<<p_virt->p_loop_me->Name()<<std::endl
                  << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                  << "Sherpa V eps^{-1} = " << yfspole << std::endl
                  << "Sherpa/One-Loop = " << yfspole / p1 << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Debugging() << "Poles cancel in YFS Virtuals to " << ncorrect
                      << " digits" << std::endl
                      << "Relative diff =  " << reldiff << std::endl;
      m_histograms1d["SinglePoleCD"]->Insert(ncorrect);
      m_histograms1d["OneLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
      m_histograms1d["OneLoopEpsLP"]->Insert(log10(fabs(p1)));
      m_histograms1d["relativediff"]->Insert(log10(fabs(reldiff)));
      msg_Debugging() << std::setprecision(32)
                      << "One-Loop Provider V eps^{-1}  = " << p1 << std::endl
                      << "Sherpa V eps^{-1}  = " << yfspole << std::endl;
    }
  }
  return m_oneloop;
}

double NLO_Base::CalculateReal() {
  if (m_coll_real)
    return p_dipoles->CalculateEEX() * m_born;
  if (!m_realtool)
    return 0;
  double real(0);
  m_real_hard1 = 0.;
  m_real_hard2 = 0.;
  m_ifi_prod = 1.;
  for (YFS::Photon &g : m_photons) {
    const Vec4D k(g.K());
    if (m_check_real_sub == 1 && (g.IsFSR() || !HasFSR())) {
      if (k.E() < 0.2 * sqrt(m_s))
        continue;
      CheckRealSub(k, 0);
    }
    double contrib;
    if (g.IsISR() && (m_isr_debug || m_fsr_debug)) {
      contrib = CalculateReal(k);
      double coll = p_dipoles->GetDipoleII().Beta1(k);
      coll /= p_dipoles->GetDipoleII().Eikonal(k);
      if (contrib != 0)
        m_histograms2d["REAL_COLL_RATIO"]->Insert(k.E(),
                                                  coll * m_born / contrib);
    } else {
      contrib = CalculateReal(k);
    }
    real += contrib;
    if (m_check_real_sub == 2)
      RecordSubScatter(k, contrib, g.IsISR() ? "realISR" : "realFSR", m_eikeex);
    g.m_beta10 = contrib;
  }
  HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta10(); },
               m_real_hard1, m_real_hard2);
  if (m_ifireal && !IsBad(m_ifi_prod)) real += m_born*(m_ifi_prod - 1.);
  return real;
}

double NLO_Base::CalculateReal(Vec4D k, bool raw) {
  double norm = 2. * pow(2 * M_PI, 3);
  Vec4D_Vector p(m_plab), pi(m_bornMomenta), pf(m_bornMomenta);
  dipoletype::code fluxtype;
  Vec4D kk = k;
  m_evts += 1;

  msg_Debugging() << METHOD << " raw=" << raw
                  << " k=" << k << " E=" << k.E() << " pt=" << k.PPerp() << "\n";

  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  fluxtype = p_nlodipoles->WhichResonant(k);

  // if (fsr_recoil) {
  // msg_Debugging() << METHOD << " FSR branch, fluxtype=" << fluxtype << "\n";
  // if (!HasFSR() && m_borngamma.size() == 0)
  // msg_Error() << "Wrong dipole type in " << METHOD << "\n";
  // YFS::DipoleView Dip_view1(p_nlodipoles->GetDipoleFF());
  // for (auto Dip = Dip_view1.begin(); Dip != Dip_view1.end(); ++Dip) {
  // double scalek = p_fsr->ScalePhoton(k);
  // Dip->SetPhotonScale(scalek);
  // Dip->AddPhotonToDipole(k);
  // if (!Dip->BoostNLO()) {
  // msg_Debugging() << METHOD << " BoostNLO failed, returning 0\n";
  // return 0;
  // }
  // int i(0);
  // for (auto f : Dip->m_flavs) {
  // p[p_nlodipoles->m_flav_label[f]] = Dip->GetNewMomenta(i++);
  // }
  // }
  // } else {
  //
  msg_Debugging() << METHOD << " fluxtype=" << fluxtype << "\n";
  MapMomenta(p, k);

  p.push_back(k);
  CheckMasses(p, 1);

  Vec4D_Vector pp = p;
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, pp, m_plab);

  double r = p_real->Calc_R(p) / norm;
  m_real = r;
  if (p_real->FailCut()) m_failcut = true;

  double flux;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == 2)
    flux = 0.5 * (p_nlodipoles->CalculateFlux(kk) + p_nlodipoles->CalculateFlux(k));
  else
    flux = p_dipoles->CalculateFlux(kk);

  double subloc = p_nlodipoles->CalculateRealSub(k);
  double subb   = p_dipoles->CalculateRealSubEEX(kk);
  m_eikeex = subb;
  m_subloc = subloc;

  msg_Debugging() << METHOD << " r=" << r << " flux=" << flux
                  << " (mode=" << m_flux_mode << ")"
                  << " subloc=" << subloc << " subb=" << subb
                  << " born=" << m_born << " alpha=" << m_rescale_alpha << "\n";

  if (!CheckMomentumConservation(p)) {
    msg_Debugging() << METHOD << " momentum conservation failed"
                    << " k.E=" << k.E() << " dip_mass=" << (p[2]+p[3]).Mass() << "\n";
    msg_Error() << "Momentum Conservation fails in " << METHOD << "\n";
    if (m_isr_debug || m_fsr_debug) {
      m_histograms1d["k_E"]->Insert(k.E());
      m_histograms1d["k_pt"]->Insert(k.PPerp());
      m_histograms1d["dip_mass"]->Insert((p[2] + p[3]).Mass());
    }
    return 0;
  }

  if ((p[2] + p[3]).Mass() < massmin)
    massmin = (p[2] + p[3]).Mass();

  if (m_isr_debug || m_fsr_debug) {
    m_histograms1d["k_E_pass"]->Insert(k.E());
    m_histograms1d["k_pt_pass"]->Insert(k.PPerp());
    m_histograms1d["dip_mass_pass"]->Insert((p[2] + p[3]).Mass());
  }

  if (IsZero(r)) {
    msg_Debugging() << METHOD << " r=0, returning 0\n";
    m_zero_real_amp++;
    return 0;
  }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << r << " flux=" << flux << "\n";
    msg_Error() << "Bad point for YFS Real\n"
                << "  Real ME : " << r << "\n"
                << "  Flux    : " << flux << "\n";
    return 0;
  }

  double tot, rcoll;
  if (m_submode == submode::local)
    tot = (r * flux - subloc * m_born / m_rescale_alpha) / subloc;
  else if (m_submode == submode::global)
    tot = (r * flux - subloc * m_born / m_rescale_alpha) / subb;
  else if (m_submode == submode::off)
    tot = (r * flux) / subb;
  else
    msg_Error() << METHOD << " unknown YFS subtraction mode " << m_submode << "\n";

  const bool ifi_above = (kk.E() > p_dipoles->IFIOmega());
  if (m_ifireal && ifi_above && m_submode == submode::global &&
      !IsZero(subb) && !IsBad(subloc) && !IsBad(subb)) {
    const double cru_gen = p_dipoles->CalculateRealSubEEX(kk);
    const double if_gen  = p_dipoles->CalculateRealSubIF(kk);
    const double ratio = (IsZero(cru_gen) || IsBad(cru_gen) || IsBad(if_gen))
                       ? 1. : 1. + if_gen/cru_gen;
    if (!IsBad(ratio)) {
      m_ifi_prod *= ratio;
      ++m_ifi_n; m_ifi_sum += ratio; m_ifi_sum2 += ratio*ratio;
      {
        const int ib = std::min(4, (int)(10.*kk.E()/sqrt(m_s)));
        if (ib >= 0) {
          ++m_ifi_x_n[ib];
          m_ifi_x_r[ib] += ratio;
          const double ex = subloc/(m_rescale_alpha*subb);
          m_ifi_x_e[ib] += IsBad(ex) ? ratio : ex;
        }
      }
      if (ratio < m_ifi_min) m_ifi_min = ratio;
      if (ratio > m_ifi_max) m_ifi_max = ratio;
    }
    msg_Debugging() << METHOD << " IFI_Real ratio=" << ratio << "\n";
  }

  msg_Debugging() << METHOD << " submode=" << m_submode
                  << " r*flux=" << r*flux
                  << " sub=" << subloc * m_born / m_rescale_alpha
                  << " tot=" << tot << "\n";

  if (m_isr_debug || m_fsr_debug) {
    double diff = ((r / subloc - m_born) - (rcoll / subb - m_born)) /
                  ((r / subloc - m_born) + (rcoll / subb - m_born));
    m_histograms1d["Real_diff"]->Insert(diff);
    if (m_isr_debug)
      m_histograms2d["Real_Flux"]->Insert(
          flux, sqrt(p_dipoles->GetDipoleII().Sprime()));
  }

  if (m_no_subtraction) {
    msg_Debugging() << METHOD << " no_subtraction: returning r/subloc=" << r/subloc << "\n";
    return r / subloc;
  }

  if (IsBad(tot)) {
    msg_Debugging() << METHOD << " tot is NaN/Inf"
                    << " r*flux=" << r*flux
                    << " subloc*born=" << subloc*m_born
                    << " subb=" << subb << "\n";
    msg_Error() << "NLO real is NaN\n"
                << "  R        : " << r << "\n"
                << "  Local  S : " << subloc * m_born << "\n"
                << "  Global S : " << subb << "\n";
  }

  if (m_isr_debug || m_fsr_debug) {
    m_histograms2d["IFI_EIKONAL"]->Insert(k.Y(), k.PPerp(),
                                          p_nlodipoles->CalculateRealSubIF(k));
    m_histograms2d["REAL_SUB"]->Insert((p[0] + p[1]).Mass(), k.E(), tot / m_born);
    m_histograms2d["REAL"]->Insert(k.E(), k.Theta(), r);
    m_histograms2d["REAL_SUB"]->Insert(k.E(), k.Theta(), tot);
  }

  sumw += tot;
  rcount += 1;
  double avg = sumw / rcount;
  if (rcount == 1000)
    m_ravg = avg;
  if (rcount > 1000) {
    double diff = fabs(1. - m_ravg / avg) * 100;
    if (diff > 10) {
      msg_Debugging() << METHOD << " large weight jump: " << diff << "%"
                      << " prev=" << m_ravg << " curr=" << avg
                      << " n=" << rcount << "\n";
      m_ravg = avg;
    }
  }

  if (raw) {
    double rawval = r * flux - subloc * m_born / m_rescale_alpha;
    msg_Debugging() << METHOD << " raw: returning " << rawval << "\n";
    return rawval;
  }

  msg_Debugging() << METHOD << " returning tot=" << tot << "\n";
  return tot;
}

double NLO_Base::CalculateRealVirtual() {
  if (!m_realvirt)
    return 0;
  m_rv_hard1 = 0.;
  m_rv_hard2 = 0.;
  if (m_rv_hard_photon==2) {
    // Hard photons are rare in the YFS spectrum, so RV_Hard_Photon=1 means
    // waiting many events for MostEnergeticPhoton() to clear the 0.2*sqrt(s)
    // CHECK_RV threshold. RV_Hard_Photon=2 exercises the RV formula/
    // subtraction on the same fixed photon every event instead.
    Vec4D k = FixedTestPhoton();
    if (m_check_rv) CheckRealVirtualSub(k);
    m_rv_hard1 = CalculateRealVirtual(k);
    m_rv_hard2 = m_rv_hard1;  // single photon: 2-photon sum == 1-photon
    return m_rv_hard1;
  }
  if (m_rv_hard_photon==1) {
    Vec4D k = MostEnergeticPhoton();
    if (k.E() == 0.) return 0;
    if (m_check_rv) {
      if (k.E() < 0.2 * sqrt(m_s)) return 0;
      CheckRealVirtualSub(k);
    }
    m_rv_hard1 = CalculateRealVirtual(k);
    m_rv_hard2 = m_rv_hard1;  // single photon: 2-photon sum == 1-photon
    return m_rv_hard1;
  }
  double realvirtual(0);
  m_rv_hard2 = 0.;
  for (YFS::Photon &g : m_photons) {
    const Vec4D k(g.K());
    if (m_check_rv) {
      if (k.E() < 0.2 * sqrt(m_s))
        continue;
      CheckRealVirtualSub(k);
    }
    double contrib = CalculateRealVirtual(k);
    realvirtual += contrib;
    g.m_beta11 = contrib;
  }
  HardestBetas(m_photons, [](const YFS::Photon &g) { return g.beta11(); },
               m_rv_hard1, m_rv_hard2);
  return realvirtual;
}

double NLO_Base::CalculateRealVirtual(Vec4D k) {
  if (!m_realvirt) return 0;
  if (m_rv_soft_cut > 0. && k.E() < m_rv_soft_cut * sqrt(m_s)) {
    m_softRV++;
    msg_Debugging() << METHOD << " skipping soft photon: E=" << k.E()
                    << " < RV_SOFT_CUT*sqrt(s)=" << m_rv_soft_cut * sqrt(m_s)
                    << "\n";
    return 0;
  }
  Vec4D_Vector p(m_plab), pi(m_bornMomenta), pf(m_bornMomenta);
  double tot(0), sub(0);
  double norm = 2 * pow(2 * M_PI, 3);
  double flux(1);
  Vec4D kk = k;
  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);
  MapMomenta(p, k);
  double yfspole;
  p.push_back(k);
  CheckMasses(p, 1);
  Vec4D_Vector pp = p;
  pp.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, pp, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, pp, m_plab);
  p_nlodipoles->p_yfsFormFact->p_virt = p_realvirt->p_loop_me;
  double subloc = p_nlodipoles->CalculateRealVirtualSubEps(k);
  yfspole = p_nlodipoles->Get_E1();
  // Eikonal factor S(k) at the reduced kinematics.
  const double eikloc = p_nlodipoles->CalculateRealSub(k);
  const double aB = eikloc * m_oneloop;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k);
  else if (m_flux_mode == 2)
    flux = 0.5 *
           (p_nlodipoles->CalculateFlux(kk) + p_nlodipoles->CalculateFlux(k));
  else
    flux = p_dipoles->CalculateFlux(kk);
  // PRINT_VAR(yfspole);
  double subb;

  // Both arms of this were the same call; the fsrcount test did nothing.
  subb = p_dipoles->CalculateRealSubEEX(kk);
  if (p.size() != (m_flavs.size() + 1)) {
    msg_Error() << "Mismatch in " << METHOD << std::endl;
  }
  double r = p_realvirt->Calc(p, m_born) / norm;
  m_rv = r * flux;
  if (p_realvirt->FailCut())
    m_failcut = true;
  ;
  if (IsBad(r)) {
    m_zeroRV++;
    msg_Error() << "Real-Virtual is " << r << std::endl;
    return 0;
  }
  if (IsZero(r,1e-30)) {
    m_zeroRV++;
    msg_Error() << "Real-Virtual is " << r << std::endl;
    return 0;
  }
  // Tree-level real-emission ME at the same mapped kinematics; needed for
  // the YFS virtual subtraction of the real-virtual, beta_1^1 =
  // [RV - B_fin*R] - S(k)*[V - B_fin*Born].
  double rtree(0.);
  if (m_realtool)
    rtree = p_real->Calc_R(p) / norm;
  else
    msg_Error() << METHOD << ": no real-emission ME available, "
                << "RV YFS subtraction is incomplete.\n";
  m_rvsub = (subloc * rtree * flux + aB) / m_rescale_alpha;
  const double d1 = r * flux - m_rvsub;
  msg_Debugging() << METHOD << " r*flux=" << r * flux
                  << " rtree*flux=" << rtree * flux
                  << " B_fin=" << subloc << " S(k)=" << eikloc
                  << " oneloop=" << m_oneloop << " d1=" << d1 << "\n";
  if (m_submode == submode::local)
    tot = d1 / eikloc;
  else if (m_submode == submode::global)
    tot = d1 / subb;
  else if (m_submode == submode::off)
    tot = (r * flux) / subb;
  const double rvmax = ATOOLS::Max(fabs(m_rv), fabs(m_rvsub));
  const double C = (rvmax > 0.) ? fabs(d1) / rvmax : 0.;
  if (m_rv_cancel_hist) {
    if (C > 0.) {
      const double logC = ATOOLS::Max(-16., log10(C));
      m_histograms1d["RV_tot_by_logC_w"]->Insert(logC, tot);
      m_histograms1d["RV_tot_by_logC_n"]->Insert(logC, 1.);
    }
    const double efrac = k.E() / sqrt(m_s);
    m_histograms1d["RV_tot_by_Efrac_w"]->Insert(efrac, tot);
    m_histograms1d["RV_tot_by_Efrac_n"]->Insert(efrac, 1.);
    const double costh_beam = k.CosTheta();
    double maxcos_leg = -1.;
    for (size_t i = 0; i < m_flavs.size(); ++i) {
      if (m_flavs[i].Charge() == 0.) continue;
      maxcos_leg = ATOOLS::Max(maxcos_leg, k.CosTheta(m_plab[i]));
    }
    const bool hardwide = (efrac > 0.1) && (maxcos_leg < 0.9);
    const double subscale =
        ATOOLS::Max(ATOOLS::Max(fabs(m_rvsub), fabs(aB)), 1e-300);
    const double lr = log10(ATOOLS::Max(fabs(m_rv) / subscale, 1e-300));
    m_histograms1d["RV_MEstab_all"]->Insert(lr, 1.);
    m_histograms1d["RV_tot_by_MEstab_w"]->Insert(lr, tot);
    if (hardwide) m_histograms1d["RV_MEstab_hardwide"]->Insert(lr, 1.);
    if (C >= 1.) m_rvHiC++;
    const double bscale = ATOOLS::Max(fabs(m_born), 1e-30);
    if (fabs(tot) > 1.e3 * bscale) {
      m_rvBlowup++;
      if (rtree == 0.) m_rvBlowupRtree0++;
      if (C >= 1.) m_rvBlowupHiC++;
      if (efrac < 0.01) m_rvBlowupSoft++;
      if (hardwide) m_rvBlowupHardWide++;
      int rank = 0;
#ifdef USING__MPI
      if (mpi->Size() > 1) rank = mpi->Rank();
#endif
      std::ofstream bf(std::string(m_debugDIR_NLO) + "/RV_blowups_rank" +
                           std::to_string(rank) + ".txt",
                       std::ios_base::app);
      bf << std::setprecision(10)
         << "E=" << k.E() << " Efrac=" << efrac
         << " costh_beam=" << costh_beam << " maxcos_leg=" << maxcos_leg
         << " pT=" << k.PPerp() << " hardwide=" << (hardwide ? 1 : 0)
         << " tot=" << tot << " d1=" << d1 << " C=" << C
         << " | rv(r*flux)=" << (r * flux) << " rvsub=" << m_rvsub
         << " rtree*flux=" << (rtree * flux) << " subloc(Bfin)=" << subloc
         << " aB(eik*oneloop)=" << aB << " eikloc(S_k)=" << eikloc
         << " subb=" << subb << " flux=" << flux << " oneloop=" << m_oneloop
         << " rescale=" << m_rescale_alpha << " submode=" << (int)m_submode
         << "\n";
    }
  }
  if (m_rv_me_max_ratio > 0.) {
    const double subscale_g =
        ATOOLS::Max(ATOOLS::Max(fabs(m_rvsub), fabs(aB)), 1e-300);
    if (fabs(m_rv) > m_rv_me_max_ratio * subscale_g) {
      m_rvUnstable++;
      msg_Debugging() << METHOD << " skipping RV, unstable loop ME: |rv|="
                      << fabs(m_rv) << " > " << m_rv_me_max_ratio
                      << "*max(|rvsub|,|aB|)=" << (m_rv_me_max_ratio * subscale_g)
                      << " (E=" << k.E() << ", d1=" << d1 << ")\n";
      return 0;
    }
  }
  if (m_rv_cancel_eps > 0. && fabs(d1) < m_rv_cancel_eps * rvmax) {
    m_softRV++;
    msg_Debugging() << METHOD << " skipping RV, cancellation C=" << C
                    << " < RV_CANCEL_EPS=" << m_rv_cancel_eps
                    << " (d1=" << d1 << ", max(|rv|,|rvsub|)=" << rvmax << ")\n";
    return 0;
  }
  if (m_check_poles == 1 && r != 0) {
    double pr1 =
        p_realvirt->p_loop_me->ME_E1() * p_realvirt->m_factor * flux / norm;
    double pr2 = p_realvirt->p_loop_me->ME_E1() * p_realvirt->m_factor;
    // p_nlodipoles->CalculateRealSub(k)*
    const double correctdigit =
        ::countMatchingDigits(pr2, -p_nlodipoles->Get_E1());
    m_histograms1d["RVSinglePoleCD"]->Insert(correctdigit);
    m_histograms1d["RealLoopEpsLP"]->Insert(log10(fabs(pr1)));
    m_histograms1d["RealLoopEpsYFS"]->Insert(log10(fabs(yfspole)));
    if (!IsEqual(pr2, -yfspole, 1e-4)) {
      msg_Out() << "Poles do not cancel in YFS Real-Virtuals" << std::endl
                << "Process =  " << p_realvirt->p_loop_me->Name() << std::endl
                << "Correct Digits =  " << correctdigit << std::endl
                << "One-Loop Provider RV eps^{-1}  = " << pr2 << std::endl
                << "Sherpa RV eps^{-1} = " << -yfspole << std::endl
                << "Sherpa/One-Loop = " << -yfspole / pr2 << std::endl;
      return 0;
    } else {
      msg_Debugging() << std::setprecision(16)
                      << "Poles cancel in YFS Real-Virtuals" << std::endl
                      << "Process =  " << p_realvirt->p_loop_me->Name()
                      << std::endl
                      << "Correct Digits =  " << correctdigit << std::endl
                      << "One-Loop Provider RV eps^{-1}  = " << pr2 << std::endl
                      << "Sherpa RV eps^{-1} = " << p_nlodipoles->Get_E1()
                      << std::endl;
    }
  }
  if (IsZero(tot)) m_zeroRV++;
  else m_nonZeroRV++;
  return tot;
}

double NLO_Base::CalculateRealReal() {
  if (!m_rrtool)
    return 0;
  double rr(0);
  m_rr_hard2 = 0.;
  const YFS::Photon_Vector &photons(m_photons);
  if (photons.size() == 0)
    return 0;
  // The two hardest photons, so the pair they form can be captured inline
  // below without a second CalculateRealReal call. Was a hand-written two-max
  // scan over indices plus Min/Max juggling to recover the ordered pair;
  // pointer identity says the same thing without the arithmetic.
  const YFS::Photon *h1 = nullptr, *h2 = nullptr;
  for (const YFS::Photon &g : photons) {
    if (!h1 || g.E() > h1->E())      { h2 = h1; h1 = &g; }
    else if (!h2 || g.E() > h2->E()) { h2 = &g; }
  }
  m_photon_pairs.clear();
  m_photon_pairs.reserve(photons.size() * (photons.size() - 1) / 2);
  for (int i = 0; i < photons.size(); ++i) {
    for (int j = i + 1; j < photons.size(); ++j) {
      Vec4D k = photons[i].K();
      Vec4D kk = photons[j].K();
      // Origin no longer selects a recoil here: both photons come off the
      // beams regardless. YFS::Photon still carries it, for the callers that
      // do care.
      double contrib = CalculateRealReal(k, kk);
      rr += contrib;
      // beta20 belongs to the pair, not to either photon: splitting it between
      // them would be an accounting convention rather than physics.
      m_photon_pairs.emplace_back(photons[i], photons[j]);
      m_photon_pairs.back().m_beta20 = contrib;
      const YFS::Photon *pi = &photons[i], *pj = &photons[j];
      if ((pi == h1 && pj == h2) || (pi == h2 && pj == h1))
        m_rr_hard2 = contrib;
      if (m_check_rr_sub == 2) {
        // accumulating scatter: record each photon of the pair with the pair
        // residual, so energetic collinear photons in large residuals show up
        RecordSubScatter(k,  contrib, "rr", m_rr_eik);
        RecordSubScatter(kk, contrib, "rr", m_rr_eik);
      }
      if (m_check_rr_sub == 1) {
        // k*=2;
        // kk*=2;
        if (k.E() < 0.2 * sqrt(m_s))
          continue;
        if (kk.E() < 0.2 * sqrt(m_s))
          continue;
        if (!m_failcut)
          CheckRealRealSub(k, kk);
      }
    }
  }
  return rr;
}

double NLO_Base::CalculateRealReal(Vec4D k1, Vec4D k2) {
  if (m_rr_soft_cut > 0.) {
    const double emin = Min(k1.E(), k2.E());
    if (emin < m_rr_soft_cut * sqrt(m_s)) {
      m_softRR++;
      msg_Debugging() << METHOD << " skipping soft photon pair: min(E1,E2)="
                      << emin << " < RR_SOFT_CUT*sqrt(s)="
                      << m_rr_soft_cut * sqrt(m_s) << "\n";
      return 0;
    }
  }
  const double norm = 2. * pow(2 * M_PI, 6);// * (m_borngamma.size()==1?1:m_borngamma.size());
  m_rr_eik = 0.;  // reset; set once the crude eikonals are computed below
  Vec4D_Vector p(m_plab);
  Vec4D_Vector pp = p;
  Vec4D kk1 = k1, kk2 = k2;

  msg_Debugging() << METHOD
                  << " k1=" << k1 << " E1=" << k1.E() << " pt1=" << k1.PPerp()
                  << " k2=" << k2 << " E2=" << k2.E() << " pt2=" << k2.PPerp() << "\n";

  // Real-real recoils both photons off the beams, the same treatment
  // CalculateReal() gives a single real emission. There is no final-state
  // branch: a photon's ISR/FSR origin does not change how the pair is
  // recoiled here.
  //
  // Three FSR branches used to sit here, selected by per-photon fsr1/fsr2
  // flags, adding the photons to the FF dipoles and recoiling with
  // BoostNLO() plus a MapInitial() beam rebuild. They were introduced in
  // error. Kept below for posterity, and because one of them - the k2-only
  // case - was missing the Dip->ClearPhotons() its two siblings had, so it
  // recoiled against photons left over from an earlier call. That is worth
  // remembering if the NNLO instability is ever traced back through here.
  //
  // if (fsr1 && !fsr2) {
  // msg_Debugging() << METHOD << " FSR branch: k1 only\n";
  // if (!HasFSR())
  // msg_Error() << "Wrong dipole type in " << METHOD << "\n";
  // YFS::DipoleView Dip_view3(p_nlodipoles->GetDipoleFF());
  // for (auto Dip = Dip_view3.begin(); Dip != Dip_view3.end(); ++Dip) {
  // Dip->ClearPhotons();
  // double scalek = p_fsr->ScalePhoton(k1);
  // Dip->SetPhotonScale(scalek);
  // Dip->AddPhotonToDipole(k1);
  // if (!Dip->BoostNLO()) {
  // msg_Debugging() << METHOD << " BoostNLO failed (fsr1), returning 0\n";
  // msg_Error() << "NLO boost failed\n";
  // return 0;
  // }
  // int i(0);
  // for (auto f : Dip->m_flavs)
  // p[p_nlodipoles->m_flav_label[f]] = Dip->GetNewMomenta(i++);
  // }
  // }
  //
  // if (!fsr1 && fsr2) {
  // msg_Debugging() << METHOD << " FSR branch: k2 only\n";
  // if (!HasFSR())
  // msg_Error() << "Wrong dipole type in " << METHOD << "\n";
  // YFS::DipoleView Dip_view4(p_nlodipoles->GetDipoleFF());
  // for (auto Dip = Dip_view4.begin(); Dip != Dip_view4.end(); ++Dip) {
  // double scalek = p_fsr->ScalePhoton(k2);
  // Dip->SetPhotonScale(scalek);
  // Dip->AddPhotonToDipole(k2);
  // if (!Dip->BoostNLO()) {
  // msg_Debugging() << METHOD << " BoostNLO failed (fsr2), returning 0\n";
  // msg_Error() << "NLO boost failed\n";
  // return 0;
  // }
  // int i(0);
  // for (auto f : Dip->m_flavs)
  // p[p_nlodipoles->m_flav_label[f]] = Dip->GetNewMomenta(i++);
  // }
  // }
  //
  // if (fsr1 && fsr2) {
  // msg_Debugging() << METHOD << " FSR branch: k1+k2\n";
  // if (!HasFSR())
  // msg_Error() << "Wrong dipole type in " << METHOD << "\n";
  // YFS::DipoleView Dip_view5(p_nlodipoles->GetDipoleFF());
  // for (auto Dip = Dip_view5.begin(); Dip != Dip_view5.end(); ++Dip) {
  // Dip->ClearPhotons();
  // double scale = p_fsr->ScalePhoton(k1) + p_fsr->ScalePhoton(k2);
  // Dip->SetPhotonScale(scale);
  // Dip->AddPhotonToDipole(k1);
  // Dip->AddPhotonToDipole(k2);
  // if (!Dip->BoostNLO()) {
  // msg_Debugging() << METHOD << " BoostNLO failed (fsr1+fsr2), returning 0\n";
  // msg_Error() << "NLO boost failed\n";
  // return 0;
  // }
  // int i(0);
  // for (auto f : Dip->m_flavs)
  // p[p_nlodipoles->m_flav_label[f]] = Dip->GetNewMomenta(i++);
  // }
  // }
  //
  // if (!fsr1 && !fsr2) {
  // msg_Debugging() << METHOD << " ISR branch: mapping momenta\n";
  // MapMomenta(p, k1, k2);
  // }
  //
  // if (fsr1 || fsr2) MapInitial(p);
  //
  MapMomenta(p, k1, k2);

  p.push_back(k1);
  p.push_back(k2);

  Vec4D_Vector _p = p;
  _p.pop_back();
  _p.pop_back();
  p_nlodipoles->MakeDipolesII(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, _p, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, _p, m_plab);

  const double subloc1 = p_nlodipoles->CalculateRealSub(k1);
  const double subloc2 = p_nlodipoles->CalculateRealSub(k2);

  double flux;
  if (m_flux_mode == 1)
    flux = p_nlodipoles->CalculateFlux(k1 + k2);
  else
    flux = p_dipoles->CalculateFlux(k1) * p_dipoles->CalculateFlux(k2);

  msg_Debugging() << METHOD << " subloc1=" << subloc1 << " subloc2=" << subloc2
                  << " flux=" << flux << " (mode=" << m_flux_mode << ")\n";

  if (!CheckMomentumConservation(p)) {
    msg_Debugging() << METHOD << " momentum conservation failed, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  double r = p_realreal->Calc_R(p) / norm;
  if (p_realreal->FailCut()) {
    msg_Debugging() << METHOD << " FailCut triggered, returning 0\n";
    m_failcut = true;
    m_zeroRR++;
    return 0;
  }
  if (IsBad(r) || IsBad(flux)) {
    msg_Debugging() << METHOD << " bad point: r=" << r << " flux=" << flux << "\n";
    m_zeroRR++;
    return 0;
  }

  p_nlodipoles->MakeDipolesII(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipolesIF(m_flavs, m_plab, m_plab);
  p_nlodipoles->MakeDipoles(m_flavs, m_plab, m_plab);


  const double sub1  = p_dipoles->CalculateRealSubEEX(kk1);
  const double sub2  = p_dipoles->CalculateRealSubEEX(kk2);
  // crude eikonal product this pair is divided by (see tot below); exposed so
  // the accumulating scatter can separate the 1/S~ blow-up from the physical
  // contribution (residual*m_rr_eik = r*flux + fullsub).
  m_rr_eik = sub1 * sub2;
  // Vec4D_Vector r1(m_plab),r2(m_plab);
  // MapMomenta(r1, kk1);
  // MapMomenta(r2, kk2);
  // r1.push_back(kk1);
  // r2.push_back(kk2);
  // const double norm1 =  2. * pow(2 * M_PI, 3);
  // const double fl1 = p_nlodipoles->CalculateFlux(kk1);
  // const double fl2 = p_nlodipoles->CalculateFlux(kk2);
  // raw: real-real subtracts for itself below, so take the unsubtracted value.
  const double real1 = CalculateReal(kk1, /*raw*/true);
  const double real2 = CalculateReal(kk2, /*raw*/true);
  // const double real1 = fl1*p_real->Calc_R(r1)/norm1-m_born*p_nlodipoles->CalculateRealSub(k1);// CalculateReal(kk1, 3 + fsr1);
  // const double real2 = fl2*p_real->Calc_R(r2)/norm1-m_born*p_nlodipoles->CalculateRealSub(k2); //CalculateReal(kk2, 3 + fsr2);
  m_recola_evts += 1;

  msg_Debugging() << METHOD << " r=" << r
                  << " sub1=" << sub1 << " sub2=" << sub2
                  << " real1=" << real1 << " real2=" << real2
                  << " born=" << m_born << "\n";

  if (IsZero(real1) || IsZero(real2)) {
    msg_Debugging() << METHOD << " real1 or real2 is zero, returning 0\n";
    m_zeroRR++;
    return 0;
  }

  const double fullsub = -subloc2 * real1 - subloc1 * real2 - subloc1 * subloc2 * m_born;
  const double tot     = (r * flux + fullsub / m_rescale_alpha) / sub1 / sub2;

  msg_Debugging() << METHOD << " fullsub=" << fullsub
                  << " r*flux=" << r * flux
                  << " tot=" << tot << "\n";

  if (IsBad(tot))
    msg_Error() << METHOD << " NNLO RR is NaN: r=" << r << " flux=" << flux
                << " fullsub=" << fullsub << " sub1=" << sub1 << " sub2=" << sub2 << "\n";

  if (!IsZero(tot)) m_nonZeroRR++;
  return tot;
}

double NLO_Base::CalculateVV() {
  if (!m_vvtool)
    return 0;
  if (m_eex_virt) {
    return p_dipoles->CalculateEEXVirtual() * m_born - m_born;
  }
  double virt;
  double sub;
  // CheckMassReg();
  if (!HasISR())
    virt = p_vv->Calc(m_bornMomenta, m_born);
  else
    virt = p_vv->Calc(m_plab, m_born);
  if (m_check_virt_born) {
    // the provider's Born is pointlike, m_born is dressed with the pion form
    // factor, so compare against the dressed provider Born
    if (!IsEqual(m_born, p_virt->p_loop_me->ME_Born()
                         * ExternalFormFactor(m_plab, m_flavs), 1e-6)) {
      msg_Error() << METHOD
                  << "\n Warning! Loop provider's born is different! YFS "
                     "Subtraction likely fails\n"
                  << "Loop Provider " << ":  " << p_virt->p_loop_me->ME_Born()
                  << "\nSherpa" << ":  " << m_born << std::endl
                  << "PhaseSpace Point = ";
      for (auto _p : m_plab)
        msg_Error() << _p << std::endl;
    }
  }
  if (p_vv->FailCut())
    return 0;
  if (m_virt_sub && p_virt->p_loop_me->Mode() != 1)
    sub = p_dipoles->CalculateVirtualSub();
  else
    sub = 0;
  double sub2 = p_dipoles->CalculateVVSubEps();
  // m_oneloop = (virt- sub * m_born/m_rescale_alpha );
  m_oneloop = (virt - sub * CalculateVirtual() / m_rescale_alpha -
               0.5 * sub * sub * m_born / m_rescale_alpha);
  if (p_virt->p_loop_me->Mode() == 1) {
    m_oneloop /= m_rescale_alpha;
    // PRINT_VAR(m_rescale_alpha);
  }
  if (IsBad(m_oneloop) || IsBad(sub)) {
    msg_Error() << "YFS Virtual is NaN" << std::endl
                << "Virtual:  " << virt << std::endl
                << "Subtraction: " << sub * m_born << std::endl
                << "PhaseSpace Point: " << std::endl
                << m_plab << std::endl;
  }
  double loope1 =
      p_vv->p_loop_me->ME_E1() *
      p_vv->m_factor; //*p_vv->m_factor;//+p_virt->p_loop_me->ME_E1()*p_virt->m_factor;;
  double loope2 =
      2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
  double yfse1 = p_dipoles->Get_E1();
  double yfse2 = p_dipoles->GetVV_E2();
  // PRINT_VAR(::countMatchingDigits(loope1, yfse1));
  // PRINT_VAR(::countMatchingDigits(loope2, -yfse2));
  // PRINT_VAR(loope1);
  // PRINT_VAR(loope1);
  // PRINT_VAR(loope1/yfse1);
  // PRINT_VAR(loope2);
  // PRINT_VAR(yfse2);
  // PRINT_VAR(loope2/yfse2);
  // PRINT_VAR(m_born);
  // PRINT_VAR(1./m_born);
  // PRINT_VAR(p_dipoles->Get_E1()+0.5*pow(p_dipoles->Get_E1(),2));
  if (m_check_poles == 1) {
    if (m_virt_sub == 0)
      sub = p_dipoles->CalculateVirtualSub();
    const double p1 = p_vv->p_loop_me->ME_E1() * p_vv->m_factor;
    const double p2 =
        2. * p_vv->p_loop_me->ME_E2() * p_vv->m_factor * p_vv->m_factor;
    const double yfspole1 = (p_dipoles->Get_E1());
    const double yfspole2 = p_dipoles->GetVV_E2();
    PRINT_VAR(p1 / yfspole1);
    int ncorrect1 = ::countMatchingDigits(p1, yfspole1, 32);
    int ncorrect2 = ::countMatchingDigits(p2, -yfspole2, 32);
    if (!IsEqual(p2, -yfspole2, 1e-6) || ncorrect1 < 10) {
      msg_Error() << "Poles do not cancel in YFS Double Virtuals" << std::endl
                  << "Correct digits \epsion^{-1} =  " << ncorrect1 << std::endl
                  << "Correct digits \epsion^{-2} =  " << ncorrect2
                  << std::endl;
      return 0;
    } else {
      int i = 0;
      msg_Debugging() << std::setprecision(32);
      msg_Out() << "Poles cancel in YFS double Virtuals to " << ncorrect2
                << " digits" << std::endl;
      m_histograms1d["SinglePoleVV"]->Insert(ncorrect1);
      m_histograms1d["DoublePoleVV"]->Insert(ncorrect2);
    }
  }
  return 0;
}

void NLO_Base::RandomRotate(Vec4D &p) {
  Vec4D t1 = p;
  // rotate around x
  p[2] = cos(m_ranTheta) * t1[2] - sin(m_ranTheta) * t1[3];
  p[3] = sin(m_ranTheta) * t1[2] + cos(m_ranTheta) * t1[3];
  t1 = p;
  // rotate around z
  p[1] = cos(m_ranPhi) * t1[1] - sin(m_ranPhi) * t1[2];
  p[2] = sin(m_ranPhi) * t1[1] + cos(m_ranPhi) * t1[2];
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k) {
  Vec4D Q;
  Vec4D QQ, PP;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  for (int i = 2; i < p.size(); ++i) {
    Q += p[i];
  }
  Q += k;
  double sq = Q.Abs2();
  Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  for (int i = 0; i < p.size(); ++i) {
    pRot.RotateBack(p[i]);
    boostQ.Boost(p[i]);
  }
  pRot.RotateBack(k);
  boostQ.Boost(k);
  double qx(0), qy(0), qz(0);
  for (int i = 2; i < p.size(); ++i) {
    qx += p[i][1];
    qy += p[i][2];
    qz += p[i][3];
  }
  if (!IsEqual(k[1], -qx, 1e-5) || !IsEqual(k[2], -qy, 1e-5) ||
      !IsEqual(k[3], -qz, 1e-5)) {
    if (k[1] > 1e-6 && k[2] > 1e-6 && k[3] > 1e-6) {
      msg_Error() << "YFS Mapping has failed for ISR\n";
      msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx << std::endl;
      msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy << std::endl;
      msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz << std::endl;
    }
  }
  for (int i = 2; i < p.size(); ++i) {
    QQ += p[i];
  }
  QQ += k;
  double sqq = QQ.Abs2();
  if (!IsEqual(sqq, sq, 1e-6)) {
    msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD
                << std::endl;
  }
  // if(m_is_isr) QQ = p[0]+p[1];
  // double zz = sqrt(sqq) / 2.;
  // double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) *
  // (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
  // Anchor the beam-0 z-orientation to the FIXED Born beam, not to p[0] after
  // boostQ. boostQ is built from Q = p[2..]+k, and p[2..] are still the Born
  // back-to-back pair here, so Q carries the photon's full transverse momentum:
  // the boost is fully 3D and p[0][3] becomes a smooth function of the photon
  // direction, flipping sign once the recoil against beam 0 is hard enough
  // (p_z' = gamma*(p_z - beta*E)). That silently swapped beams 0/1 for those
  // events, mirroring the real ME's forward-backward asymmetry and showing up
  // as an excess AFB in cos(theta) of the outgoing leptons. Matches the
  // convention Dipole::BoostNLO() uses for the subtraction terms
  // (Dipole.C:245), which anchors to m_bornmomenta[0][3] with no boost
  // involved -- so numerator and subtraction now agree on beam labelling.
  double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  // p[0] = {zz, 0, 0, z};
  // p[1] = {zz, 0, 0, -z};
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamRaw = sqq*sqq + sqr(m1*m1) + sqr(m2*m2)
                  - 2.*sqq*m1*m1 - 2.*sqq*m2*m2 - 2.*m1*m1*m2*m2;
  if (lamRaw < 0.)
    msg_Error()<<METHOD<<"(): below-threshold Kaellen argument = "<<lamRaw
               <<" (sqq = "<<sqq<<", threshold = "<<sqr(m1+m2)<<")"<<std::endl;
  double lamCM = 0.5 * sqrt(Lambda(sqq, m1 * m1, m2 * m2) / sqq);
  double E1 = lamCM * sqrt(1 + m1 * m1 / sqr(lamCM));
  double E2 = lamCM * sqrt(1 + m2 * m2 / sqr(lamCM));
  p[0] = {E1, 0, 0, sign_z * lamCM};
  p[1] = {E2, 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (int i = 0; i < p.size(); ++i) {
    pRot2.Rotate(p[i]);
    boostLab.BoostBack(p[i]);
  }
  pRot2.Rotate(k);
  boostLab.BoostBack(k);
}

void NLO_Base::MapMomenta(Vec4D_Vector &p, Vec4D &k1, Vec4D &k2) {
  Vec4D Q;
  Vec4D QQ, PP;
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  for (int i = 2; i < p.size(); ++i) {
    Q += p[i];
  }
  Q += k1 + k2;
  double sq = Q.Abs2();
  Poincare boostQ(Q);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  for (int i = 0; i < p.size(); ++i) {  
    pRot.RotateBack(p[i]);
    boostQ.Boost(p[i]);
  }
  pRot.RotateBack(k1);
  boostQ.Boost(k1);
  pRot.RotateBack(k2);
  boostQ.Boost(k2);
  double qx(0), qy(0), qz(0);
  for (int i = 2; i < p.size(); ++i) {
    qx += p[i][1];
    qy += p[i][2];
    qz += p[i][3];
  }
  // if (!IsEqual(k[1], -qx, 1e-5) || !IsEqual(k[2], -qy, 1e-5) ||
  // !IsEqual(k[3], -qz, 1e-5) ) { 	if( k[1]> 1e-6 && k[2]> 1e-6 && k[3]> 1e-6 ){
  // 		msg_Error() << "YFS Mapping has failed for ISR\n";
  // 		msg_Error() << " Photons px = " << k[1] << "\n Qx = " << -qx <<
  // std::endl; 		msg_Error() << " Photons py = " << k[2] << "\n Qy = " << -qy <<
  // std::endl; 		msg_Error() << " Photons pz = " << k[3] << "\n Qz = " << -qz <<
  // std::endl;
  // 	}
  // 	}
  for (int i = 2; i < p.size(); ++i) {
    QQ += p[i];
  }
  QQ += k1 + k2;
  double sqq = QQ.Abs2();
  if (!IsEqual(sqq, sq, 1e-6)) {
    msg_Error() << "YFS Real mapping not conserving momentum in " << METHOD
                << std::endl;
  }
  // if(m_is_isr) QQ = p[0]+p[1];
  // double zz = sqrt(sqq) / 2.;
  // double z = zz * sqrt((sqq - sqr(m_flavs[0].Mass() - m_flavs[1].Mass())) *
  // (sqq - sqr(m_flavs[0].Mass() + m_flavs[1].Mass()))) / sqq;
  // Anchor the beam-0 z-orientation to the FIXED Born beam, not to p[0] after
  // boostQ. boostQ is built from Q = p[2..]+k, and p[2..] are still the Born
  // back-to-back pair here, so Q carries the photon's full transverse momentum:
  // the boost is fully 3D and p[0][3] becomes a smooth function of the photon
  // direction, flipping sign once the recoil against beam 0 is hard enough
  // (p_z' = gamma*(p_z - beta*E)). That silently swapped beams 0/1 for those
  // events, mirroring the real ME's forward-backward asymmetry and showing up
  // as an excess AFB in cos(theta) of the outgoing leptons. Matches the
  // convention Dipole::BoostNLO() uses for the subtraction terms
  // (Dipole.C:245), which anchors to m_bornmomenta[0][3] with no boost
  // involved -- so numerator and subtraction now agree on beam labelling.
  double sign_z = (m_bornMomenta[0][3] < 0 ? -1 : 1);
  // p[0] = {zz, 0, 0, z};
  // p[1] = {zz, 0, 0, -z};
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamRaw = sqq*sqq + sqr(m1*m1) + sqr(m2*m2)
                  - 2.*sqq*m1*m1 - 2.*sqq*m2*m2 - 2.*m1*m1*m2*m2;
  if (lamRaw < 0.)
    msg_Error()<<METHOD<<"(): below-threshold Kaellen argument = "<<lamRaw
               <<" (sqq = "<<sqq<<", threshold = "<<sqr(m1+m2)<<")"<<std::endl;
  double lamCM = 0.5 * sqrt(Lambda(sqq, m1 * m1, m2 * m2) / sqq);
  double E1 = lamCM * sqrt(1 + m1 * m1 / sqr(lamCM));
  double E2 = lamCM * sqrt(1 + m2 * m2 / sqr(lamCM));
  p[0] = {E1, 0, 0, sign_z * lamCM};
  p[1] = {E2, 0, 0, -sign_z * lamCM};
  Poincare pRot2(m_bornMomenta[0], Vec4D(0., 0., 0, 1.));
  for (int i = 0; i < p.size(); ++i) {
    pRot2.Rotate(p[i]);
    boostLab.BoostBack(p[i]);
  }
  pRot2.Rotate(k1);
  pRot2.Rotate(k2);
  boostLab.BoostBack(k1);
  boostLab.BoostBack(k2);
}

void NLO_Base::MapInitial(Vec4D_Vector &p) {
  Vec4D QQ;
  Vec4D_Vector born = p;
  for (int i = 2; i < p.size(); ++i) {
    QQ += p[i];
  }
  double sqq = QQ.Abs2();
  double sign_z = (p[0][3] > 0 ? -1 : 1);
  double m1 = m_flavs[0].Mass();
  double m2 = m_flavs[1].Mass();
  double lamCM = 0.5 * sqrt(Lambda(sqq, m1 * m1, m2 * m2) / sqq);
  double E1 = lamCM * sqrt(1 + m1 * m1 / sqr(lamCM));
  double E2 = lamCM * sqrt(1 + m2 * m2 / sqr(lamCM));
  p[0] = {E1, 0, 0, sign_z * lamCM};
  p[1] = {E2, 0, 0, -sign_z * lamCM};
  Poincare boostLab(QQ);
  Poincare pRot = Poincare(p[0], Vec4D(0., 0., 0., 1.));
  for (int i = 0; i < 2; ++i) {
    pRot.Rotate(p[i]);
    boostLab.BoostBack(p[i]);
    // pRot2.Rotate(p[i]);
  }
}

void NLO_Base::CheckMasses(Vec4D_Vector &p, int realmode) {
  bool allonshell = true;
  std::vector<double> masses;
  Flavour_Vector flavs = m_flavs;
  if (realmode >= 1)
    flavs.push_back(Flavour(kf_photon));
  if (realmode >= 2)
    flavs.push_back(Flavour(kf_photon));
  if (p.size() != flavs.size())
    msg_Error() << "Mismatch between mass and flavour vectors in " << METHOD
                << std::endl;
  for (int i = 0; i < p.size(); ++i) {
    masses.push_back(flavs[i].Mass());
    if (!IsEqual(p[i].Mass(), flavs[i].Mass()) && flavs[i].Mass() != 0) {
      // msg_Debugging() << "Wrong particle masses in YFS Mapping" << std::endl
      //                 << "Flavour = " << flavs[i]
      //                 << ", with mass = " << flavs[i].Mass() << std::endl
      //                 << "Four momentum = " << p[i]
      //                 << ", with mass = " << p[i].Mass() << std::endl;
      allonshell = false;
    }
  }
  if (!allonshell) {
    m_stretcher.StretchMomenta(p, masses);
    // for (int i = 0; i < p.size(); ++i) {
      // msg_Debugging() << "Flavour = " << flavs[i]
      //                 << ", with mass = " << flavs[i].Mass() << std::endl
      //                 << "Four momentum = " << p[i]
      //                 << ", with new mass = " << p[i].Mass() << std::endl;
    // }
  }
}

void NLO_Base::RescaleMasses(Vec4D_Vector &p, std::vector<double> masses) {
  bool allonshell = true;
  if (p.size() != masses.size())
    msg_Error() << "Mismatch between mass and vectors in " << METHOD
                << std::endl;
  m_stretcher.StretchMomenta(p, masses);
  // return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k) {
  for (int i = 0; i < m_plab.size(); ++i) {
    if (m_flavs[i].IsChargedLepton()) {
      double sik = (k + m_plab[i]).Abs2();
      if (sik  < m_hardmin*m_plab[i].Abs2()) {
        msg_Out() << "Rejecting photon k = " << k << std::endl
                  << "sik = " << sik << std::endl;
        return false;
      }
    }
  }
  return true;
}

bool NLO_Base::CheckPhotonForReal(const Vec4D &k, const Vec4D_Vector &p) {
  for (int i = 0; i < p.size(); ++i) {
    if (m_flavs[i].IsChargedLepton()) {
      double sik = (k + p[i]).Abs2();
      if (sik < m_hardmin * p[i].Abs2()) {
        msg_Out() << "Rejecting photon k = " << k << std::endl
                  << "sik = " << sik << std::endl;
        return false;
      }
      // if(p[i].PPerp() < m_hardmin) return false;
    }
  }
  // if(k.PPerp() < m_hardmin) return false;
  return true;
}

bool NLO_Base::CheckMomentumConservation(Vec4D_Vector p) {
  Vec4D incoming = p[0] + p[1];
  Vec4D outgoing;
  for (int i = 2; i < p.size(); ++i) {
    if (p[i].E() < 0 || IsBad(p[i].E())) {
      msg_Error() << "Energy less than zero!: " << p[i] << std::endl;
      return false;
    }
    outgoing += p[i];
  }
  Vec4D diff = incoming - outgoing;
  if (!IsEqual(incoming, outgoing, 1e-8)) {
    msg_Error() << METHOD << std::endl
                << "Momentum not conserverd in YFS NLO" << std::endl
                << "Incoming momentum = " << incoming << std::endl
                << "Outgoing momentum = " << outgoing << std::endl
                << "Difference = " << diff << std::endl
                << "Vetoing Event " << std::endl;
    return false;
  }
  return true;
}



namespace {
// Accumulates a soft-photon subtraction scan (energy vs |residual|/Born) and
// prints a convergence summary - the residual must vanish as the photon(s)
// soften; this is the pass/fail signal the raw per-point dump doesn't give
// you without plotting it first.
}









Vec4D NLO_Base::MostEnergeticPhoton() const {
  Vec4D hardest;
  for (const auto &k : m_ISRPhotons)
    if (k.E() > hardest.E()) hardest = k;
  for (const auto &k : m_FSRPhotons)
    if (k.E() > hardest.E()) hardest = k;
  return hardest;
}

// Deterministic stand-in for MostEnergeticPhoton(), for CHECK_RV validation:
// hard photons are rare in the YFS spectrum, so selecting on
// MostEnergeticPhoton() means waiting many events for one hard enough to
// exercise CalculateRealVirtual()/CheckRealVirtualSub(). This instead builds
// a photon of fixed energy fraction (RV_TEST_PHOTON_X, of sqrt(s)/2) and
// direction (RV_TEST_PHOTON_THETA/PHI) in the Born CMS, then rotates/boosts
// it into the same frame as m_bornMomenta, matching the tail of MapMomenta.
// One-shot dump of beta_0 and beta_1 at a single, fully specified phase-space
// point, for the number-for-number comparison against KKMC's CEEX
// (Test/SherpaCompare/kkmc_ceex_crosscheck.cxx in the KKMC repo). Enabled with
// YFS: CEEX_Compare: 1.
//
// Why a hook inside a running Sherpa rather than a standalone driver: unlike
// YFS_Form_Factor and Dipole (pure functions of their arguments, hence the
// FSR/IFI harnesses), beta_1 needs p_real and p_virt wired to real ME
// providers, which means the whole process/model/generator stack. Reusing the
// live wiring is far cheaper and cannot drift from what production does.
//
// The Born configuration is whatever the phase-space generator produced for
// this event, so the momenta are PRINTED - feed them to the KKMC driver so both
// sides evaluate at exactly the same point. Only the photon is deterministic
// (FixedTestPhoton, from RV_TEST_PHOTON_X/THETA/PHI).


Vec4D NLO_Base::FixedTestPhoton() const {
  double E = m_rv_test_x * sqrt(m_s) / 2.;
  double st = sin(m_rv_test_theta), ct = cos(m_rv_test_theta);
  Vec4D k(E, E * st * cos(m_rv_test_phi), E * st * sin(m_rv_test_phi),
         E * ct);
  Poincare pRot(m_bornMomenta[0], Vec4D(0., 0., 0., 1.));
  Poincare boostLab(m_bornMomenta[0] + m_bornMomenta[1]);
  pRot.Rotate(k);
  boostLab.BoostBack(k);
  return k;
}

