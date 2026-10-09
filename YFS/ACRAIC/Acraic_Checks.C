/*!
  \file Acraic_Checks.C

  Internal ACRAIC diagnostics: soft factor and U/V matrices against
  hep-ph/0006359, Born shape, closures and partition probes.
*/

#include "YFS/ACRAIC/Acraic_Base.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "METOOLS/Main/Spin_Structure.H"
#include "PHASIC++/Process/Process_Base.H"

#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "MODEL/Main/Running_AlphaQED.H"
#include "EXTAMP/External_ME_Interface.H"
#include "PHASIC++/Process/External_ME_Args.H"

using namespace YFS;


// |Sfactor|^2 against the closed form -e^2 Q^2 (pa/k.pa - pb/k.pb)^2 / 2.
void Acraic_Base::CheckSoftFactor(const Vec4D &k)
{
  const Vec4D &pa(m_momenta[0]), &pb(m_momenta[1]);
  const double kpa(k*pa), kpb(k*pb);
  if (IsZero(kpa) || IsZero(kpb)) return;
  const Vec4D d(pa/kpa - pb/kpb);
  const double rhs(-0.5 * sqr(m_e * m_qe) * d.Abs2());
  if (rhs <= 0.) return;
  for (int h(-1); h <= 1; h += 2) {
    const Complex sf(Sfactor(pa, pb, k, h));
    const double lhs(std::norm(sf));
    const double rel(std::abs(lhs - rhs) / (std::abs(lhs) + std::abs(rhs)));
    if (rel > m_sfacworst) m_sfacworst = rel;
    ++m_sfacn;
  }
}



// U/V at equal legs against diag(b_sigma), and the helicity flip of eq.(222).
void Acraic_Base::CheckUVDiagonality(const Vec4D &k)
{
  const Vec4D &p(m_momenta[0]);
  for (int h(-1); h <= 1; h += 2) {
    Amplitude AU, AV;
    UGamma(p, p, k, h, AU);
    VGamma(p, p, k, h, AV);
    // b_sigma(k,p) = sqrt(2) * Xi(p,k) * s_sigma(k,phat), as Sfactor inlines it
    const Complex b(sqrt(2.) * Xi(p, k) *
                    (h < 0 ? Sminus(k, p) : Splus(k, p)));
    const double scale(std::abs(b) > 0. ? std::abs(b) : 1.);
    for (int a(0); a <= 1; ++a)
      for (int c(0); c <= 1; ++c) {
        const Complex u(AU.m_U[a][c]), v(AV.m_V[a][c]);
        if (a == c) {
          m_uvdiag = Max(m_uvdiag, std::abs(u - b)/scale);
          m_uvdiag = Max(m_uvdiag, std::abs(v - b)/scale);
        } else {
          m_uvoff = Max(m_uvoff, (std::abs(u) + std::abs(v))/scale);
        }
      }
    if (m_uvn == 0)
      std::cerr<<"@@@ UVDIAG sigma="<<h
               <<" U00="<<AU.m_U[0][0]<<" U11="<<AU.m_U[1][1]
               <<" U01="<<AU.m_U[0][1]<<" U10="<<AU.m_U[1][0]
               <<" V00="<<AV.m_V[0][0]<<" V11="<<AV.m_V[1][1]
               <<" b="<<b<<std::endl;
    ++m_uvn;
  }
  // b_{-sigma} = -(b_sigma)*
  const Complex bp(sqrt(2.) * Xi(p, k) * Splus(k, p));
  const Complex bm(sqrt(2.) * Xi(p, k) * Sminus(k, p));
  const double sc(std::abs(bp) > 0. ? std::abs(bp) : 1.);
  m_bsig = Max(m_bsig, std::abs(bm + conj(bp))/sc);

  const Vec4D &q(m_momenta.size() > 2 ? m_momenta[2] : m_momenta[1]);
  const double m1(m_mass_I), m2(m_mass_F);
  if (!IsZero(m_zeta*p) && !IsZero(m_zeta*q)) {
    const Complex flip(sqrt(2.) * (m2*Xi(p, q) - m1*Xi(q, p)));
    const double fs(std::abs(flip) > 0. ? std::abs(flip) : 1.);
    for (int h(-1); h <= 1; h += 2) {
      Amplitude AU, AV;
      UGamma(p, q, k, h, AU, m1, m2);
      VGamma(p, q, k, h, AV, m1, m2);
      // (-+) is m_U[1][0] for sigma=+1, (+-) is m_U[0][1] for sigma=-1
      const int a(h > 0 ? 1 : 0), c(h > 0 ? 0 : 1);
      m_uvflip = Max(m_uvflip, std::abs(AU.m_U[a][c] - flip)/fs);
      m_uvflip = Max(m_uvflip, std::abs(AV.m_V[c][a] + flip)/fs);
      // the OTHER off-diagonal slot must vanish in both
      m_uvoff  = Max(m_uvoff, (std::abs(AU.m_U[c][a])
                               + std::abs(AV.m_V[a][c]))/fs);
      ++m_uvflipn;
    }
  }
}


// Born shape fit a + b cos + c cos^2 for e+e- -> f fbar: c/a must be 1.
void Acraic_Base::AccumulateBornShape(double cth, double val)
{
  const double x[3] = {1., cth, cth*cth};
  for (int i(0); i < 3; ++i) {
    m_fitb[i] += x[i]*val;
    for (int j(0); j < 3; ++j) m_fit[i][j] += x[i]*x[j];
  }
  m_fity += val; m_fity2 += val*val; ++m_fitn;
}


void Acraic_Base::ReportBornShape()
{
  // No shape fit for Bhabha (t-channel pole); the other reports still run.
  double a(0.), b(0.), c(0.);
  bool doshape(!m_bhabha && m_fitn >= 10);
  if (doshape) {
    double M[3][4];
    for (int i(0); i < 3; ++i) {
      for (int j(0); j < 3; ++j) M[i][j] = m_fit[i][j];
      M[i][3] = m_fitb[i];
    }
    for (int i(0); i < 3 && doshape; ++i) {
      int pv(i);
      for (int r(i); r < 3; ++r) if (dabs(M[r][i]) > dabs(M[pv][i])) pv = r;
      for (int cc(0); cc < 4; ++cc) std::swap(M[i][cc], M[pv][cc]);
      if (M[i][i] == 0.) { doshape = false; break; }
      for (int r(0); r < 3; ++r) if (r != i) {
        const double f(M[r][i]/M[i][i]);
        for (int cc(0); cc < 4; ++cc) M[r][cc] -= f*M[i][cc];
      }
    }
    if (doshape) {
      a = M[0][3]/M[0][0]; b = M[1][3]/M[1][1]; c = M[2][3]/M[2][2];
      if (a == 0.) doshape = false;
    }
  }
  if (m_uvn > 0)
    msg_Out()<<"ACRAIC U/V vs eq.(diagonality): worst off-diagonal "<<m_uvoff
             <<", worst |diag - b_sigma| "<<m_uvdiag
             <<", worst |b_-s + conj(b_s)| "<<m_bsig
             <<"  over "<<m_uvn<<" evaluations"<<std::endl;
  bool any(false);
  for (size_t n(0); n < 16; ++n) if (m_pzn[n]) any = true;
  if (any) {
    msg_Out()<<"ACRAIC partition sum, propagator spread max|propZ|/min|propZ| "
             <<"across the partitions of one event:\n";
    for (size_t n(0); n < 16; ++n)
      if (m_pzn[n])
        msg_Out()<<"    n="<<n<<"  events "<<m_pzn[n]
                 <<"  mean "<<m_pzsum[n]/m_pzn[n]
                 <<"  worst "<<m_pzspread[n]<<"\n";
    msg_Out()<<"    (spread 1 => B is partition-independent and the sum "
             <<"factorises to prod_j (Sini+Sfin), O(n) not O(2^n))\n";
    msg_Out()<<"  if photons below x_cut were factorised instead of "
             <<"partitioned, the propagator they are evaluated at moves by:\n";
    for (int c(0); c < 4; ++c)
      if (m_pzsoftn[c])
        msg_Out()<<"    x_cut="<<m_pzxcut[c]<<"  events "<<m_pzsoftn[c]
                 <<"  mean "<<m_pzsoftsum[c]/m_pzsoftn[c]
                 <<"  worst "<<m_pzsoft[c]
                 <<"  partitions saved/event "
                 <<m_pzsoftsav[c]/m_pzsoftn[c]<<"\n";
  }
  if (m_maxnphot > 0)
    msg_Out()<<"ACRAIC partition sum: highest photon multiplicity seen "
             <<m_maxnphot<<" over "<<m_nstages<<" stages (soft photons "
             <<"collapsed, not enumerated); partition count wrong on "
             <<m_partbad<<" of "<<m_partn<<" events"<<std::endl;
  if (m_uvflipn > 0)
    msg_Out()<<"ACRAIC U/V helicity-flip vs eq.(222) at DISTINCT legs: worst "
             <<"relative deviation "<<m_uvflip<<" over "<<m_uvflipn
             <<" evaluations"<<std::endl;
  if (m_sfacn > 0)
    msg_Out()<<"ACRAIC soft factor vs hep-ph/0006359 eq.(soft-fac-isr): worst "
             <<"relative deviation "<<m_sfacworst<<" over "<<m_sfacn
             <<" evaluations"<<std::endl;
  if (m_rho0sum != 0.) {
    msg_Out()<<"ACRAIC O(alpha) correction <rho1/rho0 - 1> = "
             <<(m_rho1sum/m_rho0sum - 1.)
             <<"  (vertex + boxes; expect a few percent, not 0)"<<std::endl;
    // virtual and real collinear logs must come with opposite signs
    msg_Out()<<"  split:  virtual-only <rhoBV/rho0 - 1> = "
             <<(m_rhobvsum/m_rho0sum - 1.)
             <<" ,  real-only <rhoBR/rho0 - 1> = "
             <<(m_rhobrsum/m_rho0sum - 1.)<<std::endl;
    // the weight divides by RhoCrud (incoherent), not rho0 (coherent)
    if (m_rhocrudsum > 0.)
      msg_Out()<<"  ACRAIC weight denominator: <rho0>/<RhoCrud> = "
               <<(m_rho0sum/m_rhocrudsum)
               <<"   (1 => the two denominators agree)"<<std::endl;
  }
  if (!doshape) return;
  msg_Out()<<"ACRAIC Born angular shape over "<<m_fitn<<" points: "
           <<"b/a = "<<b/a<<" (forward-backward), c/a = "<<c/a
           <<" (must be 1 for the (1+cos^2) structure)"<<std::endl;
  if (!IsEqual(c/a, 1., 1e-3))
    msg_Error()<<METHOD<<"(): ACRAIC Born cos^2 coefficient is "<<c/a
               <<", not 1 - the Born spin amplitude is wrong."<<std::endl;
}

// ACRAIC: SOFT_NORM_CHECK, on the softest photon
void Acraic_Base::RunSoftNormCheck()
{
  static const bool sn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["SOFT_NORM_CHECK"].Get<int>() != 0);
  if (sn && !m_allphotons.empty()) {
    size_t js(0);
    for (size_t i(1); i < m_allphotons.size(); ++i)
      if (m_allphotons[i][0] < m_allphotons[js][0]) js = i;
    SoftNormCheck(m_allphotons[js]); }
}

void Acraic_Base::RunSoftLimitTest()
{
  static const bool sl(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["SOFT_LIMIT_TEST"].Get<int>() != 0);
  if (sl) SoftLimitTest();
}

void Acraic_Base::RunBeta1Closure()
{
  static const bool bc(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_CLOSURE"].Get<int>() != 0);
  if (bc) Beta1Closure();
}

// Partition count must be m_nstages^(photons not on a fixed stage).
void Acraic_Base::CheckPartitionCount(const long nparts)
{
  const size_t n(m_allphotons.size());
  size_t nfree(0);
  if (m_fixedstage.size() == n)
    for (size_t j(0); j < n; ++j) if (!m_fixedstage[j]) ++nfree;
  long want(1);
  if (HasFSR() && m_nstages > 1) {
    for (size_t j(0); j < nfree && want <= (1L << 40); ++j) want *= m_nstages;
    if (want > (1L << 40)) want = nparts;       // overflow guard only
  }
  if (nparts != want) {
    ++m_partbad;
    static bool warned(false);
    if (!warned) {
      warned = true;
      msg_Error()<<METHOD<<"(): partition counter enumerated "<<nparts
                 <<" partitions for "<<n<<" photons ("<<nfree
                 <<" enumerated over "<<m_nstages<<" stages), expected "
                 <<want<<". Reported once."<<std::endl;
    }
  }
  ++m_partn;
}

// ACRAIC: VIRT_PARTITION_CHECK: per-event spread of V(h) over partitions.
void Acraic_Base::CheckVirtualPartitions(const int nparts)
{
  static const bool vp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["VIRT_PARTITION_CHECK"].Get<int>() != 0);
  static int nv(0);
  if (vp && nv < 8 && !m_vpmin.empty()) {
    ++nv;
    double worst(0.); int nl(0);
    for (size_t f(0); f < m_vpmin.size(); ++f) {
      if (m_vpmax[f] < m_vpmin[f]) continue;
      ++nl;
      const double mid(0.5*(m_vpmax[f] + m_vpmin[f]));
      if (mid > 0.) worst = Max(worst, (m_vpmax[f] - m_vpmin[f])/mid);
    }
    // also the range: a spread of 2 means the minimum reached zero
    double gmin(1e30), gmax(-1e30);
    for (size_t f(0); f < m_vpmin.size(); ++f) {
      if (m_vpmax[f] < m_vpmin[f]) continue;
      gmin = Min(gmin, m_vpmin[f]); gmax = Max(gmax, m_vpmax[f]);
    }
    std::cerr<<"@@@ VIRTPART nparts="<<nparts<<" nphot="<<m_allphotons.size()
             <<" nlive="<<nl<<" worst_spread="<<worst
             <<" |V| in ["<<gmin<<","<<gmax<<"]"
             <<" sQ="<<m_svarQ<<std::endl;
  }
  m_vpmin.clear(); m_vpmax.clear();
}

// ACRAIC: BETA1_CHECK: Comix's beta_1 against the partition-summed m_realphot[j].
void Acraic_Base::CheckBeta1Ratio()
{
  static const bool b1(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_CHECK"].Get<int>() != 0);
  static int nd(0);
  if (b1 && nd < 14)
    for (size_t j(0); j < m_allphotons.size() && nd < 14; ++j) {
      if (j >= m_cxratiook.size() || !m_cxratiook[j]) continue;
      ++nd;
      double nr(0.), ncs(0.);
      const int nh(Amplitude::NHel());
      for (int f = 0; f < nh; ++f) {
        nr  += std::norm(m_realphot[j].m_A[f]);
        ncs += std::norm(m_cxcsub[j].m_A[f]);
      }
      std::cerr<<"@@@ BETA1REF Ek="<<m_allphotons[j][0]
               <<" |realphot|="<<sqrt(nr)
               <<" |Csub|="<<sqrt(ncs)
               <<" ratio="<<(nr>0.? sqrt(ncs/nr) : -1.)<<std::endl;
    }
}

// ACRAIC: BETA1_PARTITION_CHECK: spread of |beta_1| over partitions
void Acraic_Base::ReportBeta1Partitions()
{
  static const bool b1chk(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["BETA1_PARTITION_CHECK"].Get<int>() != 0);
  if (b1chk && m_b1n > 1 && m_b1min > 0.)
    { const double mean(m_b1sum/m_b1n);
      const double var(m_b1sq/m_b1n - mean*mean);
      std::cerr<<"@@@ B1PART nphot="<<m_allphotons.size()
               <<" nparts="<<m_b1n
               <<" spread="<<(m_b1max/m_b1min)
               <<" maxovmean="<<(mean>0.?m_b1max/mean:-1.)
               <<" cv="<<(mean>0.?sqrt(var>0.?var:0.)/mean:-1.)
               <<std::endl; }
}

// ACRAIC: CLOSURE_1PHOT: m_AmpExpo1 against Comix's exact one-photon amplitude.
void Acraic_Base::CheckClosureOnePhoton()
{
  static const bool cl(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["CLOSURE_1PHOT"].Get<int>() != 0);
  if (cl && m_allphotons.size() == 1 && m_cxbalignok) {
    Amplitude Mex;
    Vec4D_Vector pp(m_pacraic);
    pp.push_back(m_allphotons[0]);
    const double rn(RealNorm());
    if (rn > 0. && ComixRealAt(pp, m_PhoHel[0], Mex, -1.)) {
      double nd(0.), ne(0.);
      for (int f(0); f < Amplitude::NHel(); ++f) {
        const Complex ex(m_cxbalign.m_A[f]*Mex.m_A[f]/rn);
        nd += std::norm(m_AmpExpo1.m_A[f] - ex);
        ne += std::norm(ex);
      }
      if (ne > 0.) {
        double nb(0.), nr(0.), n1(0.), nbr(0.);
        for (int f(0); f < Amplitude::NHel(); ++f) {
          nb += std::norm(m_snapBorn.m_A[f]); nr += std::norm(m_snapReal.m_A[f]);
          n1 += std::norm(m_AmpExpo1.m_A[f]);
          nbr += std::norm(m_snapBorn.m_A[f] + m_snapReal.m_A[f]); }
        std::cerr<<"@@@ CLOS1 xg="
                 <<(m_s>0.?2.*m_allphotons[0][0]/sqrt(m_s):-1.)
                 <<" rel="<<sqrt(nd/ne)
                 <<" |Born|2="<<nb<<" |Real|2="<<nr<<" |Born+Real|2="<<nbr
                 <<" |A1|2="<<n1<<" |Mex|2="<<ne<<" nparts="<<m_nparts<<std::endl; }
    }
  }
}

// ACRAIC: KKMC_FLUX_EMULATION: the FF-pair flux terms 2 k_j.k_l/Q^2 times B
void Acraic_Base::EmulateKKMCFlux()
{
  static const int kkflux(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["KKMC_FLUX_EMULATION"].Get<int>());
  if (kkflux && m_order == 2 && m_flavs.size() == 4) {
    const Vec4D Q(m_pacraic[m_if1] + m_pacraic[m_if2]);
    const double Q2(Q.Abs2());
    double rall(0.), rsoft(0.);
    for (size_t j(0); j < m_allphotons.size(); ++j)
      for (size_t l(j + 1); l < m_allphotons.size(); ++l)
        if (!m_stagereduces[m_stage[j]] && !m_stagereduces[m_stage[l]]) {
          const double t(2.*(m_allphotons[j]*m_allphotons[l])/Q2);
          rall += t;
          if (!Beta2PairEligible((int)j, (int)l)) rsoft += t;
        }
    for (int f(0); f < Amplitude::NHel(); ++f) {
      m_AmpFluxAll.m_A[f]  += m_Sprod*rall*m_partborn0.m_A[f];
      m_AmpFluxSoft.m_A[f] += m_Sprod*rsoft*m_partborn0.m_A[f];
    }
  }
}

// ACRAIC: BETA1_PARTITION_CHECK: accumulate |beta_1| of photon 0 per partition
void Acraic_Base::AccumulateBeta1Partition()
{
  static const bool b1chk(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                          ["BETA1_PARTITION_CHECK"].Get<int>() != 0);
  if (b1chk && !m_allphotons.empty()) {
    Amplitude B1;
    if (ComixBeta1At(m_allphotons[0], m_PhoHel[0], m_sp, B1)) {
      double n2(0.);
      for (int f(0); f < Amplitude::NHel(); ++f) n2 += std::norm(B1.m_A[f]);
      const double nb(sqrt(n2));
      if (nb > 0.) {
        if (m_b1n == 0) { m_b1min = m_b1max = nb; }
        else { m_b1min = Min(m_b1min, nb); m_b1max = Max(m_b1max, nb); }
        m_b1sum += nb; m_b1sq += nb*nb;
        ++m_b1n;
      }
    }
  }
}

void Acraic_Base::RunBeta1Scan()
{
  static const bool bs(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_SCAN"].Get<int>() != 0);
  if (bs) Beta1Scan();
}

// ACRAIC: BETA1_TRACE: trace beta_1 on the event the CHECK_XS dump writes.
void Acraic_Base::SetBeta1Trace()
{
  static const bool tr(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_TRACE"].Get<int>() != 0);
  static const double dx(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["DUMP_XMIN"].Get<double>());
  static const size_t dn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                         ["DUMP_NPHOT"].Get<int>());
  const double xg(!m_allphotons.empty() && m_momenta.size() >= 2 ?
                  2.*m_allphotons[0][0]/(m_momenta[0]+m_momenta[1]).Mass()
                  : 0.);
  static const double heavy(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                            ["TRACE_FACTOR_ABOVE"].Get<double>());
  m_b1trace = (tr && m_checkxs && !m_acraicdumped
               && m_allphotons.size() == dn && xg > dx)
              || (heavy > 0. && !m_allphotons.empty());
}

// ACRAIC: FRAME_PROBE: momentum balance of ACRAIC's input, per photon family.
void Acraic_Base::ProbeFrames()
{
  static const int fp(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                      ["FRAME_PROBE"].SetDefault(0).Get<int>());
  static int nfp(0);
  // FRAME_PROBE: 1 prints the first 40 events, N > 1 the first N
  if (fp && nfp < (fp > 1 ? fp : 40) && m_pacraic.size() >= 4) {
    ++nfp;
    Vec4D P(m_pacraic[0] + m_pacraic[1]), KI, KF;
    for (size_t i(2); i < m_pacraic.size(); ++i) P -= m_pacraic[i];
    for (const Vec4D &k : m_isrphotons) KI += k;
    for (const Vec4D &k : m_fsrphotons) KF += k;
    const Vec4D r(P - KI - KF);
    std::cerr<<std::setprecision(5)<<"@@@ FRAME nI="<<m_isrphotons.size()
             <<" nF="<<m_fsrphotons.size()
             <<" |KI|="<<KI[0]<<" |KF|="<<KF[0]
             <<" resid=("<<r[0]<<","<<r[1]<<","<<r[2]<<","<<r[3]<<")"
             <<" resid_noF3="<<(P-KI)[3]<<" resid_noI3="<<(P-KF)[3]
             <<" Pfin3="<<(m_pacraic[0]+m_pacraic[1]-P)[3]
             <<" KI3="<<KI[3]<<" KF3="<<KF[3]<<std::endl;
    // off-shellness per leg: Comix puts leaves on shell without re-balancing
    std::cerr<<std::setprecision(3)<<"@@@ FRAME shell";
    for (size_t i(0); i < m_pacraic.size(); ++i)
      std::cerr<<" ["<<i<<"] "<<(m_pacraic[i].Abs2() - sqr(m_flavs[i].Mass()));
    for (const Vec4D &k : m_isrphotons) std::cerr<<" I "<<k.Abs2();
    for (const Vec4D &k : m_fsrphotons) std::cerr<<" F "<<k.Abs2();
    std::cerr<<" |resid|="<<std::setprecision(3)
             <<sqrt(sqr(r[0])+sqr(r[1])+sqr(r[2])+sqr(r[3]))
             <<" before_repair="<<m_repairresid<<std::endl;
    // one photon: generator Born point against ACRAIC's BornLegsAt(P - k)
    if (m_isrphotons.size() + m_fsrphotons.size() == 1 && m_bornmomenta.size() >= 4) {
      Vec4D_Vector pb;
      const Vec4D k(m_isrphotons.empty() ? m_fsrphotons[0] : m_isrphotons[0]);
      if (BornLegsAt(m_pacraic[0] + m_pacraic[1] - k, pb)) {
        std::cerr<<std::setprecision(6)<<"@@@ FRAME born1";
        for (size_t i(0); i < 4 && i < pb.size(); ++i)
          std::cerr<<" gen["<<i<<"]="<<m_bornmomenta[i]<<" red["<<i<<"]="<<pb[i];
        std::cerr<<std::endl;
        // Comix Born squared at both points, against m_born
        Amplitude Ag, Ar;
        Vec4D_Vector pg(m_bornmomenta.begin(), m_bornmomenta.begin() + Min(m_bornmomenta.size(), m_pacraic.size()));
        double ng(0.), nr(0.);
        if (pg.size() == m_pacraic.size() && ComixBornAmplitude(pg, Ag, NULL, -1., -1.))
          for (int f(0); f < Amplitude::NHel(); ++f) ng += std::norm(Ag.m_A[f]);
        if (ComixBornAmplitude(pb, Ar, NULL, -1., -1.))
          for (int f(0); f < Amplitude::NHel(); ++f) nr += std::norm(Ar.m_A[f]);
        const double x(2.*k[0]/sqrt(m_s));
        // and at the generator's s' point m_plabmom (lab frame, invariant)
        double np(-1.);
        if (m_plabmom.size() == m_pacraic.size()) {
          Amplitude Ap;
          if (ComixBornAmplitude(m_plabmom, Ap, NULL, -1., -1.)) {
            np = 0.; for (int f(0); f < Amplitude::NHel(); ++f) np += std::norm(Ap.m_A[f]); }
          std::cerr<<"@@@ FRAME plab";
          for (size_t i(0); i < m_plabmom.size(); ++i) std::cerr<<" ["<<i<<"]="<<m_plabmom[i];
          std::cerr<<std::endl;
        }
        std::cerr<<"@@@ FRAME bornsq x="<<x<<" m_born="<<m_born
                 <<" B2plab/m_born="<<(m_born!=0.?np/m_born:-1.)
                 <<" B2gen="<<ng<<" B2red="<<nr
                 <<" B2gen/m_born="<<(m_born!=0.?ng/m_born:-1.)
                 <<" B2red/m_born="<<(m_born!=0.?nr/m_born:-1.)
                 <<" B2red/m_born*(1-x)="<<(m_born!=0.?nr/m_born*(1.-x):-1.)<<std::endl;
      }
    }
  }
}

// ACRAIC: BORNNORM_CHECK: ACRAIC's Born at the physical pair's s' against m_born.
void Acraic_Base::CheckBornNorm()
{
  static const bool bn(ATOOLS::Settings::GetMainSettings()["ACRAIC"]["BORNNORM_CHECK"].Get<int>() != 0);
  if (bn && m_bornmomenta.size() >= 4 && m_pacraic.size() >= 4) {
    const Vec4D Qb(m_bornmomenta[m_if1] + m_bornmomenta[m_if2]);
    const Vec4D Qp(m_pacraic[m_if1] + m_pacraic[m_if2]);
    Vec4D q2(m_pacraic[m_if1]), b0(m_pacraic[0]);
    Poincare cm(Qp);
    cm.Boost(q2); cm.Boost(b0);
    const double n1(Vec3D(q2).Abs()), n2(Vec3D(b0).Abs());
    const double cth(n1 > 0. && n2 > 0. ? (Vec3D(q2)*Vec3D(b0))/(n1*n2) : 0.);
    const double sp_save(m_sp);
    m_sp = Qp.Abs2();
    MakeProp();
    Amplitude AB;
    BornAmplitude(m_pacraic, AB, -1., -1., -1);
    double sum(0.);
    const int nh(Amplitude::NHel());
    for (int f = 0; f < nh; ++f) {
      const Complex a(m_e * m_e * AB.m_A[f]);
      sum += std::real(a * conj(a));
    }
    m_sp = sp_save;
    MakeProp();
    const double avg(sum / 4.);
    std::cerr << "@@@ BORNANG m_born=" << m_born
              << " sphys=" << Qp.Abs2() << " sborn=" << Qb.Abs2()
              << " cth=" << cth
              << " acraicborn=" << avg
              << " ratio=" << (m_born != 0. ? avg/m_born : 0.) << std::endl;
  }
}

// ACRAIC: BETA1_CHECK: Comix beta_1 against the hand-coded one and soft limits.
void Acraic_Base::Beta1Check(const size_t j, const int nh, const Amplitude &C,
                             const Amplitude &Csub, const Amplitude &H,
                             const Amplitude &Bc, const Complex &Sk)
{
  static const bool d1(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BETA1_CHECK"].Get<int>() != 0);
  static int nd(0);
  if (d1 && nd < 14) { ++nd;
    double nc(0.), nh2(0.), ncs(0.);
    for (int f = 0; f < nh; ++f) {
      nc  += std::norm(C.m_A[f]);
      ncs += std::norm(Csub.m_A[f]);
      nh2 += std::norm(H.m_A[f]);
    }
    // soft theorem, phase-free: |M_1|^2 / (S~ |B|^2) -> 1
    double snum(0.), sden(0.);
    for (int f = 0; f < nh; ++f) {
      snum += std::norm(C.m_A[f]);
      sden += std::norm(Sk*Bc.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
    }
    // |M_1|^2 over the closed-form initial eikonal: must be flat when soft
    const Vec4D q1(m_pacraic[0]), q2(m_pacraic[1]), kk(m_allphotons[j]);
    const double mel(m_flavs[0].Mass());
    const double pk1(q1*kk), pk2(q2*kk), p12(q1*q2);
    const double Scl(pk1 != 0. && pk2 != 0.
                     ? 2.*p12/(pk1*pk2) - mel*mel/(pk1*pk1)
                                        - mel*mel/(pk2*pk2) : 0.);
    double nb(0.);
    for (int f = 0; f < nh; ++f)
      nb += std::norm(Bc.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
    msg_Error()<<"SOFTRAT Ek="<<m_allphotons[j][0]
             <<" x="<<(2.*m_allphotons[j][0]/sqrt(m_s))
             <<" |M1|^2/(S~|B|^2)="<<(sden>0.? snum/sden : -1.)
             <<" eik_flat="<<((Scl != 0. && nb > 0.)
                              ? snum/(Scl*nb) : -1.)<<std::endl;
    // direct subtraction M_1 - S B with Comix's Born
    double nsb(0.), ndir(0.);
    for (int f = 0; f < nh; ++f) {
      // Bc is raw Comix, C already carries m_comixnorm
      const Complex SB(RealNorm()*Sk*Bc.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
      nsb  += std::norm(SB);
      ndir += std::norm(C.m_A[f] - SB);
    }
    msg_Error()<<"BETA1DIR Ek="<<m_allphotons[j][0]
             <<" |M1|="<<sqrt(nc)<<" |S*B|="<<sqrt(nsb)
             <<" |M1|/|S*B|="<<(nsb>0.? sqrt(nc/nsb) : -1.)
             <<" |M1-S*B|/|H|="<<(nh2>0.? sqrt(ndir/nh2) : -1.)<<std::endl;
    msg_Error()<<"BETA1 Ek="<<m_allphotons[j][0]
             <<" |M1_comix|="<<sqrt(nc)
             <<" |Csub|="<<sqrt(ncs)
             <<" |H_hand|="<<sqrt(nh2)
             <<" Csub/H="<<(nh2>0.? sqrt(ncs/nh2) : -1.)<<std::endl; }
}

// YFS: ME_PROBE: the momenta ACRAIC hands Comix for its one-photon M_1.
void Acraic_Base::ProbeComixMomenta(const Vec4D_Vector &p,
                                    const std::vector<METOOLS::Spin_Amplitudes> *amps,
                                    const size_t ng, const size_t nl,
                                    const double worst)
{
  static const bool mp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
  if (mp && ng == 1 && amps && !amps->empty()) {
    const Vec4D &g(p.back());
    double best(1e99), el(0.);
    for (size_t i(2); i < nl && i < m_flavs.size(); ++i) {
      if (!m_flavs[i].IsChargedLepton()) continue;
      const double ct(Vec3D(p[i])*Vec3D(g)/(Vec3D(p[i]).Abs()*Vec3D(g).Abs()));
      const double th(acos(Max(-1., Min(1., ct)))/(m_flavs[i].Mass()/p[i][0]));
      if (th < best) { best = th; el = p[i][0]; }
    }
    double all(0.);
    for (size_t j(0); j < (*amps)[0].size(); ++j) all += std::norm((*amps)[0][j]);
    std::ostringstream o;
    o<<std::setprecision(6)<<"@@@ ACRAICPP x="<<2.*g[0]/sqrt(m_s)
     <<std::setprecision(10)<<" th_pp="<<best<<" El_pp="<<el
     <<" bal="<<worst<<" sumA2="<<all<<" g="<<g<<"\n";
    std::cerr<<o.str();
  }
}

// ACRAIC: PARTITION_BORN_CHECK: partition Born against the pair rebuilt at m_sp.
void Acraic_Base::CheckPartitionBorn(const Amplitude &AmpBorn)
{
  static const bool bx(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["PARTITION_BORN_CHECK"].Get<int>() != 0);
  static int nb(0);
  if (bx && nb < 400 && m_pacraic.size() >= 4 && m_sp > 0.) {
    const double m3(m_flavs[m_if1].Mass()), m4(m_flavs[m_if2].Mass());
    const double rs(sqrt(m_sp));
    if (rs > m3 + m4) {
      // the pair rebuilt at this partition's invariant, same direction
      Vec4D d3(m_pacraic[m_if1]);
      Poincare c0(m_pacraic[m_if1] + m_pacraic[m_if2]);
      c0.Boost(d3);
      const double n3(d3.PSpat());
      if (n3 > 0.) {
        ++nb;
        const double E3((m_sp + m3*m3 - m4*m4)/(2.*rs));
        const double pm(sqrt(Max(0., E3*E3 - m3*m3)));
        Vec4D q3(E3, pm*d3[1]/n3, pm*d3[2]/n3, pm*d3[3]/n3);
        Vec4D q4(rs - E3, -q3[1], -q3[2], -q3[3]);
        const double lcm(0.5*rs);
        const double sgn(m_bornmomenta[0][3] < 0 ? -1. : 1.);
        Vec4D_Vector rp{Vec4D(lcm,0.,0., sgn*lcm), Vec4D(lcm,0.,0.,-sgn*lcm),
                        q3, q4};
        Amplitude Br;
        BornAmplitude(rp, Br, -1., -1., -1);
        double na(0.), nb2(0.), nd(0.);
        const int nh(Amplitude::NHel());
        for (int f = 0; f < nh; ++f) {
          na  += std::norm(AmpBorn.m_A[f]);
          nb2 += std::norm(Br.m_A[f]);
          nd  += std::norm(AmpBorn.m_A[f] - Br.m_A[f]);
        }
        // Comix on ACRAIC's non-conserving legs, at this partition's scale
        Amplitude Bc;
        const bool okc(ComixBornAmplitude(m_pacraic, Bc, NULL, m_sp));
        double nc(0.); bool fin(true);
        for (int f = 0; f < Amplitude::NHel(); ++f) {
          nc += std::norm(Bc.m_A[f]);
          if (!std::isfinite(Bc.m_A[f].real())
              || !std::isfinite(Bc.m_A[f].imag())) fin = false;
        }
        const Vec4D bal(m_pacraic[0]+m_pacraic[1]-m_pacraic[m_if1]-m_pacraic[m_if2]);
        // |A_comix| = 2 e^2 |A_hand|, so 1 iff Comix used this partition's s
        double na2(0.);
        for (int f = 0; f < Amplitude::NHel(); ++f)
          na2 += std::norm(AmpBorn.m_A[f]);
        const double pred(2.*m_e*m_e*sqrt(na2));
        std::cerr<<std::setprecision(14)
                 <<"@@@ NONCONS ok="<<okc<<" finite="<<fin
                 <<" sp="<<m_sp<<" sQ="<<m_svarQ
                 <<" |B_comix|="<<sqrt(nc)
                 <<" 2e2|B_acraic|="<<pred
                 <<" ratio="<<(pred>0.? sqrt(nc)/pred : -1.)
                 <<" E_imbal="<<bal[0]
                 <<std::setprecision(6)<<std::endl;
        std::cerr<<"@@@ PARTBORN sbeam="<<(m_pacraic[0]+m_pacraic[1]).Abs2()
                 <<" sp="<<m_sp<<" sQ="<<m_svarQ
                 <<" |B_acraic|="<<sqrt(na)<<" |B_reduced|="<<sqrt(nb2)
                 <<" |diff|="<<sqrt(nd)
                 <<" rel="<<(na>0.? sqrt(nd/na) : -1.)<<std::endl;
      }
    }
  }
}

// ACRAIC: BORN_ALIGN_CHECK: is the hand/Comix alignment the same at s as at sQ?
void Acraic_Base::CheckBornAlignment(const Amplitude &H, const Amplitude &C,
                                     const int nh, const double sp_save)
{
  static const bool ck(ATOOLS::Settings::GetMainSettings()["ACRAIC"]
                       ["BORN_ALIGN_CHECK"].Get<int>() != 0);
  static int nd(0);
  if (ck && nd < 60) { ++nd;
    const double s2(m_s);           // the full beam invariant
    const double sv(sp_save);
    m_sp = s2; MakeProp(false);
    Amplitude H2, C2;
    BornAmplitude(m_pacraic, H2, -1., -1., -1);
    const bool ok2(ComixBornAmplitude(m_pacraic, C2, NULL, s2));
    m_sp = sv; MakeProp();
    if (ok2) {
      // cut each amplitude against its own norm, not the one at sQ
      double nc2n(0.), nh2n(0.);
      for (int f = 0; f < nh; ++f) {
        nc2n += std::norm(C2.m_A[f]);
        nh2n += std::norm(m_e*m_e*H2.m_A[f]);
      }
      nc2n = sqrt(nc2n); nh2n = sqrt(nh2n);
      double amin(-1.), amax(-1.), pmin(0.), pmax(0.);
      int nl(0);
      for (int f = 0; f < nh; ++f) {
        const Complex C2f(C2.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))]);
        // 1e-3 excludes the mass-suppressed slots (ratio is a mass artefact)
        if (std::abs(C2f) <= 1e-3*nc2n || std::abs(m_cxbalign.m_A[f]) == 0.
            || std::abs(m_e*m_e*H2.m_A[f]) <= 1e-3*nh2n)
          continue;
        const Complex a2((m_e*m_e*H2.m_A[f])/C2f);
        const Complex r(a2/m_cxbalign.m_A[f]);
        const double aa(std::abs(r)), pp(std::arg(r));
        if (nl++ == 0) { amin = amax = aa; pmin = pmax = pp; }
        else { amin = Min(amin,aa); amax = Max(amax,aa);
               pmin = Min(pmin,pp); pmax = Max(pmax,pp); }
      }
      // per helicity, at the physical scale, with and without the mask
      static int nr(0);
      if (nr < 2) { ++nr;
        for (int f = 0; f < nh; ++f)
          std::cerr<<"@@@ BORNHEL f="<<f
                   <<" |e2H|="<<std::abs(m_e*m_e*H.m_A[f])
                   <<" |C_masked|="
                   <<std::abs(C.m_A[FlipLegs(f, ComixLegFlips(Amplitude::s_nlegs))])
                   <<" |C_raw|="<<std::abs(C.m_A[f])<<std::endl;
        // which mask lines the moduli up here?
        int bm(-1); double bmet(1e30);
        double sh(0.), sc(0.);
        for (int f = 0; f < nh; ++f) { sh += std::norm(m_e*m_e*H.m_A[f]);
                                       sc += std::norm(C.m_A[f]); }
        sh = sqrt(sh); sc = sqrt(sc);
        for (int m = 0; m < nh; ++m) {
          double met(0.);
          for (int f = 0; f < nh; ++f)
            met += sqr(std::abs(m_e*m_e*H.m_A[f])/sh
                       - std::abs(C.m_A[FlipLegs(f, m)])/sc);
          if (bm < 0 || met < bmet) { bmet = met; bm = m; }
        }
        std::cerr<<"@@@ BORNMASK inuse="<<ComixLegFlips(Amplitude::s_nlegs)
                 <<" best="<<bm<<" metric="<<sqrt(bmet)<<std::endl; }
      std::cerr<<"@@@ BORNALIGN sQ="<<m_svarQ<<" s="<<s2<<" nlive="<<nl
               <<" |a(s)/a(sQ)|=["<<amin<<","<<amax<<"]"
               <<" arg=["<<pmin<<","<<pmax<<"]"<<std::endl;
    } }
}

// YFS: ME_PROBE: one-photon rho_1, rho_0, sum_f |M_1|^2 and the crude.
void Acraic_Base::ProbeMakeRho(const int nh)
{
  static const bool mp(ATOOLS::Settings::GetMainSettings()["YFS"]
                       ["ME_PROBE"].SetDefault(0).Get<int>() != 0);
  if (mp && NPhot() == 1) {
    double m1(0.);
    for (int f = 0; f < nh; ++f) m1 += std::norm(m_comixM1.m_A[f]);
    std::ostringstream o;
    o<<std::setprecision(10)<<"@@@ ACRAICMP x="<<std::setprecision(6)
     <<2.*m_allphotons[0][0]/sqrt(m_s)<<std::setprecision(10)
     <<" rho1="<<m_result<<" rho0="<<m_result0<<" M1sq="<<m1/4.
     <<" rhocr="<<m_rhocrud;
    // final-stage |s_F|^2 (physical and pre-emission legs) and pseudo-flux
    for (int g(0); g < m_nstages && g < (int)m_Sfac.size(); ++g) {
      if (g == m_initstage || m_Sfac[g].empty()) continue;
      Vec4D q; Complex spre(0., 0.);
      bool ok(m_prefsr.size() == m_flavs.size());
      for (size_t l(0); l < m_stagelegs[g].size(); ++l) {
        const int lg(m_stagelegs[g][l].leg);
        if (lg < (int)m_pacraic.size()) q += m_pacraic[lg];
        if (ok && lg < (int)m_prefsr.size())
          spre += m_stagelegs[g][l].w * SfactorLeg(m_prefsr[lg], m_allphotons[0], m_PhoHel[0]);
        else ok = false;
      }
      const double q2(q.Abs2()), pf(q2 > 0. ? (q + m_allphotons[0]).Abs2()/q2 : 0.);
      o<<" g="<<g<<" sFpost2="<<std::norm(m_Sfac[g][0])<<" pflux="<<pf
       <<" sFpre2="<<(ok ? std::norm(spre) : -1.);
    }
    if (m_initstage < (int)m_Sfac.size() && !m_Sfac[m_initstage].empty())
      o<<" sI2="<<std::norm(m_Sfac[m_initstage][0]);
    o<<"\n";
    std::cerr<<o.str();
  }
}
