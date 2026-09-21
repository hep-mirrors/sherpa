/*!
  \file Ceex_Checks.C

  Self-checks that need no second generator: the soft factor and the U/V
  matrices against the closed forms of hep-ph/0006359, and the Born angular
  shape. Each reports once at end of run.

  Split out of Ceex_Base.C, which had reached 1867 lines. The pieces are
  grouped by what they compute, not by call order.
*/

#include "YFS/CEEX/Ceex_Base.H"
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


void Ceex_Base::CheckSoftFactor(const Vec4D &k)
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



void Ceex_Base::CheckUVDiagonality(const Vec4D &k)
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


void Ceex_Base::AccumulateBornShape(double cth, double val)
{
  const double x[3] = {1., cth, cth*cth};
  for (int i(0); i < 3; ++i) {
    m_fitb[i] += x[i]*val;
    for (int j(0); j < 3; ++j) m_fit[i][j] += x[i]*x[j];
  }
  m_fity += val; m_fity2 += val*val; ++m_fitn;
}


void Ceex_Base::ReportBornShape()
{
  // The (1+cos^2) fit is an s-channel statement. Bhabha's t-channel pole makes
  // it ill-conditioned (it comes back NaN), so skip it there - but do not let
  // that skip the checks below, which apply to both.
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
    msg_Out()<<"CEEX U/V vs eq.(diagonality): worst off-diagonal "<<m_uvoff
             <<", worst |diag - b_sigma| "<<m_uvdiag
             <<", worst |b_-s + conj(b_s)| "<<m_bsig
             <<"  over "<<m_uvn<<" evaluations"<<std::endl;
  {
    bool any(false);
    for (size_t n(0); n < 16; ++n) if (m_pzn[n]) any = true;
    if (any) {
      msg_Out()<<"CEEX partition sum, propagator spread max|propZ|/min|propZ| "
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
  }
  if (m_maxnphot > 0)
    msg_Out()<<"CEEX partition sum: highest photon multiplicity seen "
             <<m_maxnphot<<" ("<<(1u<<Min(m_maxnphot,size_t(30)))
             <<" partitions); 2^n enumeration wrong on "<<m_partbad
             <<" of "<<m_partn<<" events"<<std::endl;
  if (m_uvflipn > 0)
    msg_Out()<<"CEEX U/V helicity-flip vs eq.(222) at DISTINCT legs: worst "
             <<"relative deviation "<<m_uvflip<<" over "<<m_uvflipn
             <<" evaluations"<<std::endl;
  if (m_sfacn > 0)
    msg_Out()<<"CEEX soft factor vs hep-ph/0006359 eq.(soft-fac-isr): worst "
             <<"relative deviation "<<m_sfacworst<<" over "<<m_sfacn
             <<" evaluations"<<std::endl;
  if (m_rho0sum != 0.) {
    msg_Out()<<"CEEX O(alpha) correction <rho1/rho0 - 1> = "
             <<(m_rho1sum/m_rho0sum - 1.)
             <<"  (vertex + boxes; expect a few percent, not 0)"<<std::endl;
    // Split by beta. The virtual carries the collinear logs ln(s/m^2); the real
    // must carry them back with the opposite sign. Same-sign entries mean the
    // cancellation is not happening and the total is not trustworthy.
    msg_Out()<<"  split:  virtual-only <rhoBV/rho0 - 1> = "
             <<(m_rhobvsum/m_rho0sum - 1.)
             <<" ,  real-only <rhoBR/rho0 - 1> = "
             <<(m_rhobrsum/m_rho0sum - 1.)<<std::endl;
    // The CEEX event weight divides by RhoCrud (incoherent), not rho0
    // (coherent). They are equal for one partition only. This ratio is how much
    // that choice is worth once FSR opens the 2^n sum.
    if (m_rhocrudsum > 0.)
      msg_Out()<<"  CEEX weight denominator: <rho0>/<RhoCrud> = "
               <<(m_rho0sum/m_rhocrudsum)
               <<"   (1 => the two denominators agree)"<<std::endl;
  }
  if (!doshape) return;
  msg_Out()<<"CEEX Born angular shape over "<<m_fitn<<" points: "
           <<"b/a = "<<b/a<<" (forward-backward), c/a = "<<c/a
           <<" (must be 1 for the (1+cos^2) structure)"<<std::endl;
  if (!IsEqual(c/a, 1., 1e-3))
    msg_Error()<<METHOD<<"(): CEEX Born cos^2 coefficient is "<<c/a
               <<", not 1 - the Born spin amplitude is wrong."<<std::endl;
}
