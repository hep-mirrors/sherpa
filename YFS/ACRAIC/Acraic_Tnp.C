/*!
  \file Acraic_Tnp.C

  ACRAIC: TNP_STREAMS - templates for theory nuisance parameters
  (Tackmann, arXiv:2411.18606) in the YFS beta-hat expansion.

  The YFS exponent is exact to all orders; what a truncated ACRAIC lacks is
  the next order of the IR-finite residues beta-hat. Each residue class c
  gets a nuisance parameter theta_c that rescales it at the next order,

      B_c -> B_c (1 + theta_c N_c) ,

  with N_c the expansion parameter of that class (chosen offline, e.g.
  (alpha/pi) L_e for the production stage). The weight is linear in
  theta_c to the order that matters,

      w(theta) = w + sum_c theta_c N_c T_c ,   T_c = 2 Re<B_c^F, F>/4 ,

  where F is the amplitude the nominal squares and B_c^F the part of F that
  B_c contributes. Every term the partition loop adds to A_1 or to the
  beta_2 sum carries its class (TnpAdd at the write sites), so sum_c B_c^F = F
  and sum_c T_c = 2|F|^2/4 exactly; ReportTnp checks the first per event.

  Approximations: the templates use the external virtual as real vhalf = v/2
  on every helicity (the helicity-resolved delta^h, beta_0^(2), Dyson VP and
  the REAL_VIRTUAL remainder are left out). They change T_c at relative
  O(alpha), which is beyond the accuracy a nuisance template needs.

  Classes: the emission stages are charge neutral, so each stage's beta_1 is
  gauge invariant on its own, but a photon's per-partition terms are not
  separately small. Each photon (all its partitions together, i.e. its whole
  beta_1) is "isr" when its production-stage eikonal is the larger, else
  "fsr". beta_2 sets (and beta_m, m > 2, under TREE_BETA_MAX) split into
  all-initial, all-final and mixed.
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
    By the photon, not by its stage in the current partition: a photon on the
    production stage carries M_1 - s_I B there and only -s_F B on a final
    stage, so a per-partition tag splits "M_1" from "eikonal subtraction" (two
    O(1) templates that cancel, measured +0.66/-0.63 of sigma at the Z), not
    initial from final radiation. Summed over the partitions, a photon's terms
    are its whole IR-finite beta_1; it goes to the class of its larger stage
    eikonal (beam- or final-collinear), the same in every partition.
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
    The same residual against sum |term| over every term added (TnpAdd).
    Summing one set of terms by class rather than all at once can differ by
    rounding times the cancellation sum|term| / |A_1|, large on rare events
    (gamma gamma multi-photon: |A_0|^2/|A_1|^2 up to ~1e7); a tagging bug
    would show here at O(1), rounding at ~1e-15.
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

void Acraic_Base::ReportTnp() const
{
  if (!m_tnp) return;
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
