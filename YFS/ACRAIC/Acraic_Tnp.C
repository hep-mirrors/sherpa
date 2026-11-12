/*!
  \file Acraic_Tnp.C

  ACRAIC: TNP_STREAMS - theory-nuisance templates (arXiv:2411.18606) in the
  YFS beta-hat expansion, B_c -> B_c (1 + theta_c N_c). Each term added to
  A_1 or the beta_2 sum carries its class (TnpAdd), so sum_c B_c^F = F and
  sum_c T_c = 2|F|^2/4; ReportTnp checks this per event. The external virtual
  enters as real v/2 on every helicity (an O(alpha) relative effect on T_c).

  Classes: a photon's whole beta_1 is "isr" when its production-stage
  eikonal is the larger, else "fsr"; beta_2 (and beta_m under TREE_BETA_MAX)
  split into all-initial, all-final and mixed.

  ACRAIC: TNP_IFI adds the initial-final interference template, the part
  odd in the final-state charge beyond beta_0 (TnpIfiTemplate):
      T_IFI = Odd[w_k] - Odd[w_0] ,   Odd[w] = (w - w^flip)/2 ,
  w^flip with s_F -> -s_F per partition, M_1^F -> -M_1^F, v_IF -> -v_IF and
  Y_IF -> -Y_IF (KKMC's KeyINT on minus off at leading order). beta_2 and
  ACRAIC's own boxes are taken as charge even.
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#ifdef USING__MPI
#include "ATOOLS/Org/My_MPI.H"
#endif

using namespace YFS;
using namespace ATOOLS;

const char *Acraic_Base::TnpClassName(int c)
{
  static const char *names[tnp_nclass] =
    {"B0", "ISR", "FSR", "MULTI", "B2I", "B2F", "B2M"};
  return c >= 0 && c < tnp_nclass ? names[c] : "?";
}

int Acraic_Base::TnpPhotonClass(size_t iphot, bool) const
{
  /*
    By the photon, not by its stage in the current partition: a per-partition
    tag would split M_1 from the eikonal subtraction, not initial from final
    radiation. A photon goes to the class of its larger stage eikonal, the same
    in every partition.
  */
  double si(0.), sf(0.);
  for (int g(0); g < m_nstages; ++g) {
    const double n(std::norm(m_Sfac[g][iphot]));
    if (g == m_initstage) si += n;
    else sf += n;
  }
  return si >= sf ? tnp_isr : tnp_fsr;
}

int Acraic_Base::TnpSetClass(const std::vector<int> &photons) const
{
  size_t ni(0);
  for (int j : photons) {
    const bool fixed(m_fixedstage.size() == m_stage.size() && m_fixedstage[j]);
    if (TnpPhotonClass(j, fixed) == tnp_isr) ++ni;
  }
  if (ni == photons.size()) return tnp_b2i;
  if (ni == 0) return tnp_b2f;
  return tnp_b2m;
}

bool Acraic_Base::TnpTemplates(int order, double vhalf, std::vector<double> &t) const
{
  t.assign(tnp_nclass, 0.);
  if (!m_tnp) return false;
  if (!m_ifi_coherent) {
    // the incoherent sum squares each partition; the class amplitudes are summed
    THROW(fatal_error, "ACRAIC: TNP_STREAMS needs the coherent partition sum "
          "(ACRAIC: IFI and YFS: IFI true).");
  }
  if (order == 2 && m_order != 2) return false;
  const int nh(Amplitude::NHel());
  // closure of the tagging, on every event
  double d1(0.), d2(0.), n1(0.), nterm(0.);
  for (int f(0); f < nh; ++f) {
    Complex s1(0., 0.), s2(0., 0.);
    for (int c(tnp_b0); c <= tnp_multi; ++c) s1 += m_tnpamp[c].m_A[f];
    for (int c(tnp_b2i); c <= tnp_b2m; ++c) s2 += m_tnpamp[c].m_A[f];
    nterm += m_tnpabs[f]*m_tnpabs[f];
    d1 += std::norm(s1 - m_AmpExpo1.m_A[f]);
    if (m_order == 2) d2 += std::norm(s2 - m_AmpBeta2.m_A[f]);
    n1 += std::norm(m_AmpExpo1.m_A[f]);
  }
  if (!(n1 > 0.)) return false;
  const double err(sqrt((d1 + d2)/n1));
  /*
    The same residual against sum |term| over every term added (TnpAdd): a
    tagging bug shows at O(1), rounding at ~1e-15.
  */
  const double errterm(nterm > 0. ? sqrt((d1 + d2)/nterm) : 0.);
  ++m_tnpn;
  if (err > m_tnpworst) m_tnpworst = err;
  if (err > 1e-8) ++m_tnpbad;
  if (errterm > m_tnpworstterm) m_tnpworstterm = errterm;
  if (errterm > 1e-12) ++m_tnpbadterm;
  /*
    The composition F the nominal squares, and how much of it each class is.
    order 1: ACRAIC_O1V1, the paper's O(alpha) column, F = (1 + vhalf) A_1;
    at ORDER 1 with REAL_VIRTUAL off the nominal is A_1 + vhalf A_0 instead.
    order 2: F = A_2 + vhalf A_v, A_v = A_1 (VirtualOnBeta1) or A_0.
  */
  const bool vona1(order == 2 ? m_vonbeta1 : (m_order == 2 || m_vonbeta1));
  const double sa1(1. + (vona1 ? vhalf : 0.));
  const double sa0(vona1 ? 0. : vhalf);
  const int ncl(order == 2 ? tnp_nclass : tnp_multi + 1);
  for (int f(0); f < nh; ++f) {
    const Complex ax(order == 2 ? m_AmpExpo2.m_A[f] : m_AmpExpo1.m_A[f]);
    const Complex F(ax + (vona1 ? vhalf*m_AmpExpo1.m_A[f] : vhalf*m_AmpExpo0.m_A[f]));
    for (int c(0); c < ncl; ++c) {
      Complex bf(c <= tnp_multi ? sa1*m_tnpamp[c].m_A[f] : m_tnpamp[c].m_A[f]);
      if (c == tnp_b0) bf += sa0*m_AmpExpo0.m_A[f];
      t[c] += 2.*std::real(std::conj(bf)*F);
    }
  }
  for (double &x : t) x /= 4.;
  return true;
}

void Acraic_Base::TnpIfiPartitionFlip()
{
  static bool checked(false);
  if (!checked) {
    checked = true;
    std::string why;
    if (m_nstages != 2)
      why = "it needs one initial and one final stage (e+e- -> f fbar), this process has "
            + ToString(m_nstages) + " stages";
    else if (m_bhabha || BornHasExchangeLine())
      why = "the Born has a t-channel exchange line, where the split by emitting stage "
            "is not gauge invariant";
    else if (WStagesActive()) why = "W stages are active";
    if (!why.empty())
      THROW(fatal_error, "ACRAIC: TNP_IFI is not available: " + why
            + ". Set ACRAIC: {TNP_IFI: false}.");
  }
  const size_t n(m_allphotons.size());
  m_tnpflipph.assign(n, Complex(1., 0.));
  Complex prod(1., 0.);
  for (size_t j(0); j < n; ++j) {
    Complex fj(1., 0.);
    if (m_fixedstage.size() == n && m_fixedstage[j]) {
      // a collapsed photon carries the total eikonal sum_g s_g
      Complex s(0., 0.), sb(0., 0.);
      for (int g(0); g < m_nstages; ++g) {
        s  += m_Sfac[g][j];
        sb += (g == m_initstage ? 1. : -1.)*m_Sfac[g][j];
      }
      if (std::abs(s) > 0.) fj = sb/s;
    }
    else if (j < m_stage.size() && m_stage[j] != m_initstage) fj = Complex(-1., 0.);
    m_tnpflipph[j] = fj;
    prod *= fj;
  }
  m_tnppartflip = prod;
}

Acraic_Base::PropShifts Acraic_Base::TnpIfiFlippedResonance(const PropShifts &sh) const
{
  PropShifts out(sh);
  for (auto &e : out)
    if (e.first & PHASIC::Process_Base::s_propshiftresonance) e.second[0] = -e.second[0];
  return out;
}

bool Acraic_Base::TnpIfiInitialM1(int iphot, const Vec4D_Vector &pp, int hel,
                                  const PropShifts &shifts, const Amplitude &M1,
                                  Amplitude &M1I)
{
  if (!m_ifi || SingleSided()) return false;
  // with the photon on the production stage, the vetoes remove its final
  // placement: the diagrams in which it comes off the final stage
  PropShifts shI(shifts);
  AddPlacementVetoes(std::vector<int>(1, iphot), shI);
  if (!ComixRealShifted(pp, hel, M1I, shI)) return false;
  if (m_tnpifichk > 0) {
    // M_1^F: the initial placement vetoed instead
    const int keep(m_stage[iphot]);
    m_stage[iphot] = 1 - m_initstage;
    PropShifts shF(shifts);
    AddPlacementVetoes(std::vector<int>(1, iphot), shF);
    m_stage[iphot] = keep;
    Amplitude M1F;
    if (ComixRealShifted(pp, hel, M1F, shF)) {
      double d(0.), nm(0.);
      for (int f(0); f < Amplitude::NHel(); ++f) {
        d  += std::norm(M1I.m_A[f] + M1F.m_A[f] - M1.m_A[f]);
        nm += std::norm(M1.m_A[f]);
      }
      if (nm > 0.) m_tnpifichkworst = Max(m_tnpifichkworst, sqrt(d/nm));
      ++m_tnpifichkn;
    }
    --m_tnpifichk;
  }
  return true;
}

bool Acraic_Base::TnpIfiTemplate(int order, double vhalf, double vif, double eflip,
                                 double &t, double *odd) const
{
  t = 0.;
  if (odd) *odd = 0.;
  if (!TnpIfiOn() || !m_ifi_coherent) return false;
  if (order == 2 && m_order != 2) return false;
  const int nh(Amplitude::NHel());
  // the composition of TnpTemplates: F = A_k + vhalf A_v, A_v = A_1 or A_0
  const bool vona1(order == 2 ? m_vonbeta1 : (m_order == 2 || m_vonbeta1));
  const double vhb(vhalf - vif);
  double w(0.), wb(0.), w0(0.), w0b(0.);
  for (int f(0); f < nh; ++f) {
    // beta_2 taken as charge even: A_2^flip = A_1^flip + beta_2
    const Complex ak(order == 2 ? m_AmpExpo2.m_A[f] : m_AmpExpo1.m_A[f]);
    const Complex akb(order == 2 ? m_tnpflip1.m_A[f] + m_AmpBeta2.m_A[f]
                                 : m_tnpflip1.m_A[f]);
    const Complex F(ak + vhalf*(vona1 ? m_AmpExpo1.m_A[f] : m_AmpExpo0.m_A[f]));
    const Complex Fb(akb + vhb*(vona1 ? m_tnpflip1.m_A[f] : m_tnpflip0.m_A[f]));
    w   += std::norm(F);
    wb  += std::norm(Fb);
    w0  += std::norm(m_AmpExpo0.m_A[f]);
    w0b += std::norm(m_tnpflip0.m_A[f]);
  }
  t = 0.5*((w - eflip*wb) - (w0 - eflip*w0b))/4.;
  if (odd) *odd = 0.5*(w - eflip*wb)/4.;
  if (order == 1) ++m_tnpifin;
  return !IsBad(t);
}

void Acraic_Base::ReportTnp() const
{
  if (!m_tnp) return;
  if (m_tnpifi) {
    double b[5] = {(double)m_tnpifin, (double)m_tnpifiunsplit, (double)m_tnpifichkn,
                   (double)m_tnpifievenbox, (double)m_tnpifievenres};
    double mx(m_tnpifichkworst);
#ifdef USING__MPI
    if (mpi->Size() > 1) {
      mpi->Allreduce(b, 5, MPI_DOUBLE, MPI_SUM);
      mpi->Allreduce(&mx, 1, MPI_DOUBLE, MPI_MAX);
    }
#endif
    msg_Out()<<"ACRAIC: TNP_IFI: template on "<<(long)b[0]<<" events; M_1 taken as "
             <<"charge even (not split) for "<<(long)b[1]<<" photon(s); own boxes taken as "
             <<"even on "<<(long)b[3]<<" partition(s)"<<std::endl;
    if (b[4] > 0.)
      msg_Out()<<"ACRAIC: TNP_IFI: greco exponent taken as charge even on "<<(long)b[4]
               <<" partition Born(s) (hand-coded Born)"<<std::endl;
    if (b[2] > 0.)
      msg_Out()<<"ACRAIC: TNP_IFI_CHECK: |M_1^I + M_1^F - M_1|/|M_1| worst "<<mx
               <<" on "<<(long)b[2]<<" photon(s)"<<std::endl;
  }
  double buf[3] = {(double)m_tnpn, (double)m_tnpbad, (double)m_tnpbadterm};
  double mx[2] = {m_tnpworst, m_tnpworstterm};
#ifdef USING__MPI
  if (mpi->Size() > 1) {
    mpi->Allreduce(buf, 3, MPI_DOUBLE, MPI_SUM);
    mpi->Allreduce(mx, 2, MPI_DOUBLE, MPI_MAX);
  }
#endif
  if (buf[0] <= 0.) return;
  msg_Out()<<"ACRAIC: TNP_STREAMS: templates on "<<(long)buf[0]<<" events; "
           <<"class tagging closure |sum_c B_c - A|/|A_1| worst "<<mx[0]
           <<", above 1e-8 on "<<(long)buf[1]<<" events"<<std::endl;
  // a separate line, so the gate's parse of the line above is unchanged
  msg_Out()<<"ACRAIC: TNP_STREAMS: closure against sum|term| (rounding scale) worst "
           <<mx[1]<<", above 1e-12 on "<<(long)buf[2]<<" events"<<std::endl;
}
