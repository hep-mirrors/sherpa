/*!
  \file NLO_LoopHelicity.C

  YFS: ACRAIC_Virtual: helicity (ACRAIC in the paper). The one-loop virtual of
  the ACRAIC column resolved by helicity, from the loop provider's
  per-helicity amplitudes (Recola 2: get_amplitude_r1_rcl after the
  compute_process_rcl call YFS.NLO already makes, so no extra loop call).

  This file builds the IR-finite complex factor delta^h per helicity
  configuration of the provider, at the Born point of the virtual, and runs
  the per-event identity check. The label map to ACRAIC's index and the
  squared weight sum_h |A_1^h + delta^h A_0^h|^2 live in
  YFS/ACRAIC/Acraic_LoopHelicity.C. Design and standalone checks:
  NOTES-acraic-loop-helicity-2026-09-26.md.
*/

#include "YFS/NLO/NLO_Base.H"
#include "YFS/Main/Define_Dipoles.H"
#include "YFS/NLO/RealVirtual.H"
#include "PHASIC++/Process/Virtual_ME2_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/My_MPI.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Math/MathTools.H"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <map>

using namespace YFS;
using namespace ATOOLS;

namespace {

  /*!
    Acraic_Base::CnuA (KKMC KKbvir::CnuA), exact in the masses, s- and
    t-channel. Duplicated because NLO_Base has no ACRAIC object to call; the
    two must stay identical.
  */
  Complex LoopCnuA(double svar, double m1, double m2)
  {
    const double mm(m1*m2);
    const Complex Mas12(mm, 0.), Nu((-svar + m1*m1 + m2*m2)/2., 0.);
    const Complex xlam(sqrt((Nu - Mas12)*(Nu + Mas12)));
    const Complex z(svar > 0. ? Mas12/(Nu - xlam) : (Nu + xlam)/Mas12);
    const Complex lz((z.imag() == 0. && z.real() <= 0.)
                     ? log(-z) + Complex(0., M_PI) : log(z));
    return Nu/xlam * lz;
  }

  Settings &Main() { return Settings::GetMainSettings(); }

}


bool NLO_Base::LoopHelicityIR(const Vec4D_Vector &p, double alpha, double mu2,
                              Complex &K, double &imsh) const
{
  return LoopHelicityIR(m_flavs, p, alpha, mu2, K, imsh);
}

bool NLO_Base::LoopHelicityIR(const Flavour_Vector &flavs, const Vec4D_Vector &p,
                              double alpha, double mu2,
                              Complex &K, double &imsh) const
{
  static const double c(Main()["ACRAIC"]["VIRT_HEL_IM_C"].SetDefault(0.5).Get<double>());
  K = Complex(0., 0.);
  imsh = 0.;
  // a neutral leg (the emitted photon of the (n+1)-point) has no pair term
  const size_t n(Min(p.size(), flavs.size()));
  for (size_t i(0); i < n; ++i) {
    const double Qi(flavs[i].Charge());
    if (Qi == 0.) continue;
    const double mi(flavs[i].Mass());
    if (!(mi > 0.)) return false;
    const double Zi((i < 2 ? 1. : -1.)*Qi);
    for (size_t j(i+1); j < n; ++j) {
      const double Qj(flavs[j].Charge());
      if (Qj == 0.) continue;
      const double mj(flavs[j].Mass());
      if (!(mj > 0.)) return false;
      const double Zj((j < 2 ? 1. : -1.)*Qj);
      const bool same((i < 2) == (j < 2));
      const double sij(same ? (p[i]+p[j]).Abs2() : (p[i]-p[j]).Abs2());
      const Complex Kij(alpha/(2.*M_PI)*(-Zi*Zj)*(LoopCnuA(sij, mi, mj) - 1.));
      K += Kij;
      if (Kij.imag() != 0. && sij != 0.)
        imsh += Kij.imag()*log(mu2*exp(c)/std::abs(sij));
    }
  }
  return true;
}


double NLO_Base::VirtualSubPhotonMass()
{
  // The photon-mass branch of Define_Dipoles::CalculateVirtualSub, whatever
  // Dim_Reg says (ISR+FSR mode).
  double sub(0.);
  YFS_Form_Factor *ff(p_dipoles->p_yfsFormFact.get());
  const double kmax(sqrt(m_s)/2.);
  for (auto &D : p_dipoles->m_set.ByType(dipoletype::initial))
    sub += D.ChargeNorm()*ff->BVirtGeneral(D, kmax);
  for (auto &D : p_dipoles->m_set.FF())
    sub += D.ChargeNorm()*ff->BVirtGeneral(D, kmax);
  for (auto &D : p_dipoles->m_set.IF())
    sub += D.ChargeNorm()*ff->BVirtGeneral(D, kmax);
  return sub;
}


void NLO_Base::BuildLoopHelicityFactors(double virt, double sub)
{
  m_lhelok = false;
  static const int probe(Main()["ACRAIC"]["VIRT_HEL_PROBE"].SetDefault(0).Get<int>());
  static const double minborn(Main()["ACRAIC"]["VIRT_HEL_MINBORN"].SetDefault(1e-3).Get<double>());
  static long nprobe(0);
  PHASIC::Virtual_ME2_Base *lme(p_virt->p_loop_me.get());
  std::vector<Complex> a0, a1;
  std::vector<std::vector<int> > hel;
  if (!lme->HelicityAmplitudes(a0, a1, hel) || a0.empty() || a0.size() != a1.size()) {
    if (m_lhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): the loop provider returned no helicity"
                 <<" amplitudes (only Recola 2 can); ACRAIC_Virtual: helicity"
                 <<" falls back to delta = v/2, i.e. to external.\n";
    return;
  }
  const double alpha(p_virt->m_factor*2.*M_PI);
  const double mu(lme->IRscale()), mu2(mu*mu);
  Complex K;
  double imsh(0.);
  if (!LoopHelicityIR(m_plab, alpha, mu2, K, imsh)) {
    if (m_lhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): a charged leg is massless; the IR coefficient"
                 <<" needs Massive: true. Falling back to delta = v/2.\n";
    return;
  }
  // V/B of the provider exactly as Virtual::Calc_V formed it, the uniform
  // running-alpha correction it subtracted, and the YFS subtraction in the
  // loop's coupling.
  const double vraw(p_virt->m_factor*lme->ME_Finite());
  const double run(vraw - virt/m_born);
  const double subs(sub/m_rescale_alpha);
  const double v(m_oneloop/m_born);
  double s0(0.), s1(0.), amax(0.);
  for (size_t c(0); c < a0.size(); ++c) {
    s0 += std::norm(a0[c]);
    s1 += 2.*std::real(std::conj(a0[c])*a1[c]);
    amax = Max(amax, std::abs(a0[c]));
  }
  if (!(s0 > 0.)) { ++m_lhel_fail; return; }
  const double clo(s1/s0 - vraw);
  const Complex common(-0.5*run - 0.5*subs, -imsh);
  std::vector<Complex> delta(a0.size()), ratio(a0.size(), Complex(0., 0.));
  std::vector<char> ratiook(a0.size(), 0);
  long nfall(0);
  double sd(0.);
  for (size_t c(0); c < a0.size(); ++c) {
    if (std::abs(a0[c]) > minborn*amax) {
      ratio[c] = a1[c]/a0[c]; ratiook[c] = 1;
      delta[c] = ratio[c] + common;
    }
    else { delta[c] = Complex(0.5*v, 0.); ++nfall; }
    if (IsBad(delta[c].real()) || IsBad(delta[c].imag())) {
      if (m_lhel_fail++ == 0)
        msg_Error()<<METHOD<<"(): NaN in the helicity amplitudes (a non-zero"
                   <<" Higgs width with external photons does this in Recola"
                   <<" 2.3.0). Falling back to delta = v/2.\n";
      return;
    }
    sd += std::norm(a0[c])*2.*delta[c].real();
  }
  const double res(sd/s0 - v);
  m_la0.swap(a0);
  m_lhel.swap(hel);
  m_ldelta.swap(delta);
  m_lratio.swap(ratio);
  m_lratio_ok.swap(ratiook);
  m_lrun = run; m_limsh = imsh;
  m_lhelmom = m_plab;
  m_lhelv = v;
  m_lhelok = true;
  ++m_lhel_n;
  m_lhel_nfall += nfall;
  m_lhel_ncfg += m_la0.size();
  m_lhel_res_sum += res; m_lhel_res_sq += res*res;
  m_lhel_res_max = Max(m_lhel_res_max, std::abs(res));
  m_lhel_clo_max = Max(m_lhel_clo_max, std::abs(clo));
  // the residual's decade, 1e-17 .. 1 in bins 1 .. 18 (0: exactly zero)
  const int decade(res != 0. ? Max(1, Min(19, (int)std::floor(log10(std::abs(res))) + 18)) : 0);
  ++m_lhel_reshist[decade];
  /*
    IR diagnostic, Dim_Reg only: the eps-form subtraction YFS.NLO used against
    the photon-mass form moved from lambda to mu_IR with the analytic K. Zero
    when BVV_full_eps with epsloop = 4 pi is the photon-mass B at lambda = mu_IR
    and the dipoles cover every charged pair. Informative only: delta itself
    uses the subtraction YFS.NLO used, so the identity above is exact.
  */
  double ir(0.);
  if (m_dim_reg && m_photonMass > 0.) {
    const double spm(VirtualSubPhotonMass()/m_rescale_alpha);
    ir = subs - spm - 2.*K.real()*log(mu2/sqr(m_photonMass));
    m_lhel_ir_sum += ir; m_lhel_ir_sq += ir*ir;
    m_lhel_ir_max = Max(m_lhel_ir_max, std::abs(ir));
  }
  if (probe > 0 && nprobe < probe) {
    ++nprobe;
    std::ostringstream o;
    o<<std::setprecision(10)
     <<"@@@ VHEL v="<<v<<" vraw="<<vraw<<" run="<<run<<" sub="<<subs
     <<" K="<<K.real()<<","<<K.imag()<<" imsh="<<imsh<<" muIR="<<mu
     <<" closure="<<clo<<" identity="<<res<<" irdiag="<<ir
     <<" nfall="<<nfall<<" ncfg="<<m_la0.size()<<"\n";
    for (size_t c(0); c < m_la0.size(); ++c) {
      if (std::abs(m_la0[c]) <= minborn*amax) continue;
      o<<"@@@ VHEL   h=(";
      for (size_t l(0); l < m_lhel[c].size(); ++l)
        o<<(l?",":"")<<(m_lhel[c][l] > 0 ? "+" : (m_lhel[c][l] < 0 ? "-" : "0"));
      o<<") w="<<std::norm(m_la0[c])/s0<<" delta="<<m_ldelta[c].real()
       <<","<<m_ldelta[c].imag()<<"\n";
    }
    std::cerr<<o.str();
  }
}


void NLO_Base::ReportLoopHelicity()
{
  // same decision on every rank, so the reductions below stay matched
  if (!m_useacraic || m_acraicvirtsrc != acraicvirt::helicity) return;
  double buf[10] = {(double)m_lhel_n, (double)m_lhel_fail, (double)m_lhel_nfall,
                    (double)m_lhel_ncfg, m_lhel_res_sum, m_lhel_res_sq,
                    m_lhel_ir_sum, m_lhel_ir_sq, 0., 0.};
  double mx[3] = {m_lhel_res_max, m_lhel_clo_max, m_lhel_ir_max};
  double hist[20];
  for (int b(0); b < 20; ++b) hist[b] = (double)m_lhel_reshist[b];
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 10, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 3, MPI_DOUBLE, MPI_MAX);
    mpi->Allreduce(hist, 20, MPI_DOUBLE, MPI_SUM);
  }
#endif
  if (buf[0] + buf[1] <= 0.) return;
  const double n(buf[0] > 0. ? buf[0] : 1.);
  const double rmean(buf[4]/n), rrms(sqrt(Max(0., buf[5]/n - rmean*rmean)));
  const double imean(buf[6]/n), irms(sqrt(Max(0., buf[7]/n - imean*imean)));
  msg_Out()<<std::setprecision(6)
           <<"ACRAIC_Virtual: helicity - "<<(long)buf[0]<<" events with helicity"
           <<" factors, "<<(long)buf[1]<<" fell back to delta = v/2 entirely.\n"
           <<"  helicity configurations below VIRT_HEL_MINBORN (delta = v/2): "
           <<(long)buf[2]<<" of "<<(long)buf[3]<<"\n"
           <<"  identity  <2 Re delta>_Born - v : mean "<<rmean<<"  rms "<<rrms
           <<"  max |.| "<<mx[0]<<"\n"
           <<"  closure   sum 2Re(A0*A1)/sum|A0|^2 - V/B : max |.| "<<mx[1]<<"\n"
           <<"  IR diag   sub_eps - sub_lambda - 2ReK ln(mu^2/lambda^2) : mean "
           <<imean<<"  rms "<<irms<<"  max |.| "<<mx[2]<<"\n"
           <<"  identity |residual| histogram (0 exactly, then decades from 1e-17):\n   ";
  for (int b(0); b < 20; ++b) msg_Out()<<" "<<(long)hist[b];
  msg_Out()<<std::endl;
}


/*
  ACRAIC: REAL_VIRTUAL helicity is requested when ACRAIC runs and the
  setting says so; the remainder needs the Born point's per-helicity ratios
  (BuildLoopHelicityFactors) whatever ACRAIC_Virtual places the Born
  virtual with.
*/
bool NLO_Base::RealVirtualHelicityRequested() const
{
  static const bool on(Main()["ACRAIC"]["REAL_VIRTUAL"].SetDefault(acraicrv::automatic)
                       .Get<acraicrv::code>() == acraicrv::helicity);
  return on && m_useacraic;
}


bool NLO_Base::BuildRealVirtualHelicityFactors(const Vec4D &k, double dv,
                                               const Vec4D_Vector &plab,
                                               RVHelicityFactors &out)
{
  out = RVHelicityFactors();
  out.k = k; out.dv = dv; out.plab = plab;
  static const int probe(Main()["ACRAIC"]["RV_HEL_PROBE"].SetDefault(0).Get<int>());
  static const double minborn(Main()["ACRAIC"]["VIRT_HEL_MINBORN"].SetDefault(1e-3).Get<double>());
  static long nprobe(0);
  if (p_realvirt == NULL || !m_lhelok || !m_vborn_own_ok) {
    if (m_rvhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): no per-helicity Born factors (the Born loop"
                 <<" provider must give helicity amplitudes, i.e. Recola);"
                 <<" ACRAIC: REAL_VIRTUAL helicity falls back to the averaged"
                 <<" remainder.\n";
    return false;
  }
  PHASIC::Virtual_ME2_Base *lme(p_realvirt->p_loop_me.get());
  const Flavour_Vector &flavs(p_realvirt->m_flavs);
  const size_t nl(flavs.size());
  if (nl != m_flavs.size() + 1 || !flavs.back().IsPhoton()
      || m_rvLastLoop.size() != nl) {
    if (m_rvhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): the real-virtual process "<<flavs
                 <<" is not the Born plus one photon last; no helicity"
                 <<" real-virtual.\n";
    return false;
  }
  std::vector<Complex> a0, a1;
  std::vector<std::vector<int> > hel;
  if (!lme->HelicityAmplitudes(a0, a1, hel) || a0.empty() || a0.size() != a1.size()) {
    if (m_rvhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): the real-virtual provider returned no helicity"
                 <<" amplitudes (only Recola 2 can); ACRAIC: REAL_VIRTUAL"
                 <<" helicity falls back to the averaged remainder.\n";
    return false;
  }
  const double alpha(p_realvirt->m_factor*2.*M_PI);
  const double mu(lme->IRscale()), mu2(mu*mu);
  Complex K;
  double imsh(0.);
  if (!LoopHelicityIR(flavs, m_rvLastLoop, alpha, mu2, K, imsh)) {
    if (m_rvhel_fail++ == 0)
      msg_Error()<<METHOD<<"(): a charged leg of the (n+1)-point is massless;"
                 <<" the IR coefficient needs Massive: true.\n";
    return false;
  }
  // the provider's loop over tree as RealVirtual::LoopOverTree formed it,
  // and what it added beyond the amplitudes' own finite part: the
  // massification constant Z (YFS: VIRTUAL_MASSIFICATION), helicity
  // independent, so 1/2 of it on every amplitude
  const double ltfin(p_realvirt->m_factor*lme->ME_Finite());
  const double mshift(m_rvLastLt - ltfin);
  // the Born configuration of each (n+1) configuration: the same labels on
  // the first n legs (the photon is last)
  std::map<std::vector<int>, size_t> bornidx;
  for (size_t c(0); c < m_lhel.size(); ++c) bornidx[m_lhel[c]] = c;
  const Complex commonB(-0.5*m_lrun - 0.5*m_vborn_own_sub, -m_limsh);
  const Complex commonR(-0.5*m_rvLastSub + 0.5*m_rvct + 0.5*mshift,
                        -imsh + 0.5*m_rvct_im);
  const double v1(m_rvLastLt - m_rvLastSub + m_rvct);   // v_{n+1} + c
  double s0(0.), amax(0.), sd(0.), sr(0.);
  for (size_t c(0); c < a0.size(); ++c) {
    s0 += std::norm(a0[c]);
    amax = Max(amax, std::abs(a0[c]));
  }
  if (!(s0 > 0.)) { ++m_rvhel_fail; return false; }
  std::vector<Complex> rem(a0.size());
  long nfall(0), nbfall(0);
  for (size_t c(0); c < a0.size(); ++c) {
    Complex dnp1;
    if (std::abs(a0[c]) > minborn*amax) dnp1 = a1[c]/a0[c] + commonR;
    else { dnp1 = Complex(0.5*v1, 0.); ++nfall; }
    Complex dB(0.5*m_vborn_own, 0.);
    const std::vector<int> hb(hel[c].begin(), hel[c].end() - 1);
    const std::map<std::vector<int>, size_t>::const_iterator it(bornidx.find(hb));
    if (it != bornidx.end() && m_lratio_ok[it->second]) dB = m_lratio[it->second] + commonB;
    else if (std::abs(a0[c]) > minborn*amax) ++nbfall;
    rem[c] = dnp1 - dB;
    if (IsBad(rem[c].real()) || IsBad(rem[c].imag())) {
      if (m_rvhel_fail++ == 0)
        msg_Error()<<METHOD<<"(): NaN in the (n+1)-point helicity amplitudes"
                   <<" (a non-zero Higgs width with external photons does this"
                   <<" in Recola 2.3.0). Falling back to the averaged remainder.\n";
      return false;
    }
    sd += std::norm(a0[c])*2.*dnp1.real();
    sr += std::norm(a0[c])*2.*rem[c].real();
  }
  const double res(sd/s0 - v1);        // the (n+1)-point identity
  const double dres(sr/s0 - dv);       // helicity resolution of the Born's v
  out.hel.swap(hel); out.a0.swap(a0); out.rem.swap(rem);
  out.ok = true;
  ++m_rvhel_n;
  m_rvhel_nfall += nfall; m_rvhel_nbornfall += nbfall;
  m_rvhel_ncfg += out.a0.size();
  m_rvhel_res_sum += res; m_rvhel_res_sq += res*res;
  m_rvhel_res_max = Max(m_rvhel_res_max, std::abs(res));
  m_rvhel_dv_sum += dres; m_rvhel_dv_sq += dres*dres;
  m_rvhel_dv_max = Max(m_rvhel_dv_max, std::abs(dres));
  if (probe > 0 && nprobe < probe) {
    ++nprobe;
    double rmax(0.), rsum(0.), wsum(0.);
    for (size_t c(0); c < out.a0.size(); ++c) {
      if (std::abs(out.a0[c]) <= minborn*amax) continue;
      const double w(std::norm(out.a0[c])/s0);
      rmax = Max(rmax, std::abs(out.rem[c]));
      rsum += w*std::norm(out.rem[c]); wsum += w;
    }
    // the photon's angle to the nearest beam, in the loop frame
    double thmin(M_PI);
    for (size_t i(0); i < 2; ++i) {
      const Vec3D a(m_rvLastLoop[i]), b(m_rvLastLoop.back());
      const double cth(a*b/(a.Abs()*b.Abs()));
      thmin = Min(thmin, acos(Max(-1., Min(1., cth))));
    }
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ RVHEL x="<<2.*k.E()/sqrt(m_s)
     <<" fsr="<<(PhotonIsFSR(k)?1:0)<<" nph="<<m_photons.size()
     <<" thbeam="<<thmin<<" dv="<<dv<<" v1="<<v1<<" vB="<<m_vborn_own
     <<" mshift="<<mshift<<" identity="<<res<<" dvdiff="<<dres
     <<" remmax="<<rmax<<" remrms="<<(wsum > 0. ? sqrt(rsum/wsum) : 0.)
     <<" imsh="<<imsh<<" imshB="<<m_limsh<<" nfall="<<nfall<<" nbfall="<<nbfall
     <<" ncfg="<<out.a0.size()<<"\n";
    for (size_t c(0); c < out.a0.size(); ++c) {
      if (std::abs(out.a0[c]) <= minborn*amax) continue;
      o<<"@@@ RVHEL   h=(";
      for (size_t l(0); l < out.hel[c].size(); ++l)
        o<<(l?",":"")<<(out.hel[c][l] > 0 ? "+" : (out.hel[c][l] < 0 ? "-" : "0"));
      o<<") w="<<std::norm(out.a0[c])/s0<<" rem="<<out.rem[c].real()
       <<","<<out.rem[c].imag()<<"\n";
    }
    std::cerr<<o.str();
  }
  return true;
}


void NLO_Base::ReportRealVirtualHelicity()
{
  // same decision on every rank, so the reductions below stay matched
  if (!RealVirtualHelicityRequested()) return;
  double buf[10] = {(double)m_rvhel_n, (double)m_rvhel_fail, (double)m_rvhel_nfall,
                    (double)m_rvhel_ncfg, (double)m_rvhel_nbornfall,
                    m_rvhel_res_sum, m_rvhel_res_sq, m_rvhel_dv_sum, m_rvhel_dv_sq, 0.};
  double mx[2] = {m_rvhel_res_max, m_rvhel_dv_max};
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 10, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 2, MPI_DOUBLE, MPI_MAX);
  }
#endif
  if (buf[0] + buf[1] <= 0.) return;
  const double n(buf[0] > 0. ? buf[0] : 1.);
  const double rmean(buf[5]/n), rrms(sqrt(Max(0., buf[6]/n - rmean*rmean)));
  const double dmean(buf[7]/n), drms(sqrt(Max(0., buf[8]/n - dmean*dmean)));
  msg_Out()<<std::setprecision(6)
           <<"ACRAIC: REAL_VIRTUAL helicity - "<<(long)buf[0]<<" photons with"
           <<" helicity remainders, "<<(long)buf[1]<<" without (averaged dv).\n"
           <<"  (n+1) configurations below VIRT_HEL_MINBORN (rem = dv/2): "
           <<(long)buf[2]<<" of "<<(long)buf[3]<<"; live ones whose Born"
           <<" configuration is below it (delta_B = v_B/2): "<<(long)buf[4]<<"\n"
           <<"  identity  <2 Re delta_{n+1}>_tree - (v_{n+1} + c) : mean "<<rmean
           <<"  rms "<<rrms<<"  max |.| "<<mx[0]<<"\n"
           <<"  <2 Re rem>_tree - dv (Born delta resolved on the real's"
           <<" helicity mix) : mean "<<dmean<<"  rms "<<drms<<"  max |.| "<<mx[1]
           <<std::endl;
}
