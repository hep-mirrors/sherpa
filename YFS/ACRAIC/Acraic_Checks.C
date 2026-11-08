/*!
  \file Acraic_Checks.C

  Internal ACRAIC diagnostics: soft factor and U/V matrices against
  hep-ph/0006359, the partition count and the end-of-run summaries.
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


void Acraic_Base::ReportBornShape()
{
  if (m_uvn > 0)
    msg_Debugging()<<"ACRAIC U/V vs eq.(diagonality): worst off-diagonal "<<m_uvoff
             <<", worst |diag - b_sigma| "<<m_uvdiag
             <<", worst |b_-s + conj(b_s)| "<<m_bsig
             <<"  over "<<m_uvn<<" evaluations"<<std::endl;
  if (m_maxnphot > 0)
    msg_Debugging()<<"ACRAIC partition sum: highest photon multiplicity seen "
             <<m_maxnphot<<" over "<<m_nstages<<" stages (soft photons "
             <<"collapsed, not enumerated); partition count wrong on "
             <<m_partbad<<" of "<<m_partn<<" events"<<std::endl;
  if (m_uvflipn > 0)
    msg_Debugging()<<"ACRAIC U/V helicity-flip vs eq.(222) at DISTINCT legs: worst "
             <<"relative deviation "<<m_uvflip<<" over "<<m_uvflipn
             <<" evaluations"<<std::endl;
  if (m_sfacn > 0)
    msg_Debugging()<<"ACRAIC soft factor vs hep-ph/0006359 eq.(soft-fac-isr): worst "
             <<"relative deviation "<<m_sfacworst<<" over "<<m_sfacn
             <<" evaluations"<<std::endl;
  if (m_rho0sum != 0.) {
    msg_Debugging()<<"ACRAIC O(alpha) correction <rho1/rho0 - 1> = "
             <<(m_rho1sum/m_rho0sum - 1.)
             <<"  (vertex + boxes; expect a few percent, not 0)"<<std::endl;
    // virtual and real collinear logs must come with opposite signs
    msg_Debugging()<<"  split:  virtual-only <rhoBV/rho0 - 1> = "
             <<(m_rhobvsum/m_rho0sum - 1.)
             <<" ,  real-only <rhoBR/rho0 - 1> = "
             <<(m_rhobrsum/m_rho0sum - 1.)<<std::endl;
    // the weight divides by RhoCrud (incoherent), not rho0 (coherent)
    if (m_rhocrudsum > 0.)
      msg_Debugging()<<"  ACRAIC weight denominator: <rho0>/<RhoCrud> = "
               <<(m_rho0sum/m_rhocrudsum)
               <<"   (1 => the two denominators agree)"<<std::endl;
  }
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
