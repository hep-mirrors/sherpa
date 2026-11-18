#include "AHADIC++/Tools/Z_Pair_Selector.H"
#include "ATOOLS/Math/Random.H"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace AHADIC;
using namespace ATOOLS;

namespace {
  // envelope: panels per parton to start with, at most maxpanels after
  // halving those where l = log f varies by more than dlog
  const int    s_npanels   = 8;
  const size_t s_maxpanels = 48;
  const double s_dlog      = 4.;
  // pair draw: no pair after this many trials, the cluster decays instead
  const int    s_maxtrials = 100000;
}

double Z_Pair_Selector::Exponent(const double gamma, const double kt02,
                                 const double scale) {
  return dabs(gamma) > 5.e-3 ? gamma * scale / kt02 : 0.;
}

double Z_Pair_Selector::Peak(const double alpha, const double beta,
                             const double c,
                             const double zmin, const double zmax) {
  // stationary point of alpha ln z + beta ln(1-z) - c/z:
  // (alpha+beta) z^2 - (alpha-c) z - c = 0
  const double A = alpha + beta;
  if (std::abs(A) > 1.e-10) {
    const double disc = sqr(alpha-c) + 4.*A*c;
    if (disc >= 0.) {
      const double root = std::sqrt(disc);
      for (const double sign : {1., -1.}) {
        const double z = ((alpha-c) + sign*root)/(2.*A);
        if (z > zmin && z < zmax) return z;
      }
    }
  }
  else if (std::abs(alpha-c) > 1.e-10) {
    const double z = c/(c-alpha);
    if (z > zmin && z < zmax) return z;
  }
  // none inside: the larger end
  const double flo = alpha*std::log(zmin) + beta*std::log1p(-zmin) - c/zmin;
  const double fhi = alpha*std::log(zmax) + beta*std::log1p(-zmax) - c/zmax;
  return fhi >= flo ? zmax : zmin;
}

Z_Pair_Selector::Z_Pair_Selector() : m_giveup(false) {}

void Z_Pair_Selector::SetRegion(const double Q2, const double kt2,
                                const double * mdec2,
                                const double * zmin, const double * zmax) {
  m_Q2  = Q2;
  m_kt2 = kt2;
  for (int k=0;k<2;++k) {
    m_mdec2[k] = mdec2[k];
    m_zmin[k]  = zmin[k];
    m_zmax[k]  = zmax[k];
    m_a[k]     = (mdec2[k]+2.*kt2)/Q2;
  }
  // live interval of z_i: the z_i that leave an allowed z_j in
  // [zmin_j,zmax_j], i.e. z_i^2 + 2 ph z_i + a_i <= 0
  for (int i=0;i<2;++i) {
    const int j = 1-i;
    const double ph = 0.5*(m_a[j]-m_a[i]-1.), disc = ph*ph-m_a[i];
    m_live[i] = false;
    if (disc<0.) continue;
    m_lo[i] = Max(Max(m_zmin[i],-ph-std::sqrt(disc)),m_a[i]/(1.-m_zmin[j]));
    m_hi[i] = Min(Min(m_zmax[i],-ph+std::sqrt(disc)),1.-m_a[j]/m_zmax[j]);
    m_live[i] = m_lo[i]<m_hi[i];
  }
  m_giveup = false;
}

void Z_Pair_Selector::Envelope(const int k, const double lo, const double hi,
                               const double alpha, const double beta,
                               const double c) {
  // panels of z_k on [lo,hi] with the largest l = log f_k on each; f_k is
  // log-concave or monotonic there, so the largest value of a panel is at an
  // end or at the peak
  auto l = [&](const double z) {
    return alpha*std::log(z)+beta*std::log1p(-z)-c/z;
  };
  const double zpeak = Peak(alpha,beta,c,lo,hi);
  const double lpeak = (zpeak>lo && zpeak<hi) ? l(zpeak) : 0.;
  std::vector<double> & edge = m_edge[k], & lmax = m_lmax[k];
  edge.clear();
  m_lz.clear();
  for (int n=0;n<=s_npanels;++n) {
    edge.push_back(n<s_npanels ? lo+(hi-lo)*n/s_npanels : hi);
    m_lz.push_back(l(edge.back()));
  }
  double top;
  for (bool split=true;;) {
    lmax.resize(edge.size()-1);
    top = -std::numeric_limits<double>::max();
    for (size_t n=0;n+1<edge.size();++n) {
      lmax[n] = Max(m_lz[n],m_lz[n+1]);
      if (zpeak>edge[n] && zpeak<edge[n+1]) lmax[n] = Max(lmax[n],lpeak);
      top = Max(top,lmax[n]);
    }
    if (!split || edge.size()>s_maxpanels) break;
    // halve the panels where l varies by more than dlog, so a narrow peak
    // costs panels, not trials; those far below the top stay
    split = false;
    m_newedge.assign(1,edge[0]);
    m_newlz.assign(1,m_lz[0]);
    for (size_t n=0;n+1<edge.size();++n) {
      if (lmax[n]-Min(m_lz[n],m_lz[n+1])>s_dlog && lmax[n]>top-30. &&
          edge.size()+m_newedge.size()<2*s_maxpanels) {
        m_newedge.push_back(0.5*(edge[n]+edge[n+1]));
        m_newlz.push_back(l(m_newedge.back()));
        split = true;
      }
      m_newedge.push_back(edge[n+1]);
      m_newlz.push_back(m_lz[n+1]);
    }
    edge = m_newedge;
    m_lz = m_newlz;
  }
  // envelope height of each panel, relative to the top
  m_fmax[k].resize(lmax.size());
  for (size_t n=0;n<lmax.size();++n) m_fmax[k][n] = std::exp(lmax[n]-top);
}

bool Z_Pair_Selector::operator()(const double * alpha, const double * beta,
                                 const double * c, double * z) {
  // z0 runs over its live interval, z1 over the range the region reaches
  // above it; a cell, i.e. a pair of panels the region crosses, is picked by
  // its envelope mass, a point in it uniformly, and accepted with f0 f1 over
  // the envelope
  if (!m_live[0]) return false;
  const double lo[2] = {m_lo[0],Lower(1,m_lo[0])};
  const double hi[2] = {m_hi[0],Upper(1,m_hi[0])};
  if (!(lo[1]<hi[1])) return false;
  for (int k=0;k<2;++k) Envelope(k,lo[k],hi[k],alpha[k],beta[k],c[k]);
  // cells: panel q of z1 for panel p of z0, cut to the window of z1 over p
  const std::vector<double> & e0 = m_edge[0], & e1 = m_edge[1];
  m_cells.clear();
  m_cellsum.clear();
  double sum = 0.;
  for (size_t p=0;p+1<e0.size();++p) {
    const double z1lo = Lower(1,e0[p]), z1hi = Upper(1,e0[p+1]);
    // first panel of z1 that ends above z1lo
    size_t q = std::upper_bound(e1.begin()+1,e1.end(),z1lo)-e1.begin()-1;
    for (;q+1<e1.size() && e1[q]<z1hi;++q) {
      const double clo = Max(e1[q],z1lo), chi = Min(e1[q+1],z1hi);
      if (!(clo<chi)) continue;
      const double w = m_fmax[0][p]*m_fmax[1][q]*(e0[p+1]-e0[p])*(chi-clo);
      if (!(w>0.)) continue;
      m_cells.push_back({p,q});
      m_cellsum.push_back(sum += w);
    }
  }
  if (m_cells.empty()) return false;
  for (int trials=0;trials<s_maxtrials;++trials) {
    const double r = ran->Get()*sum;
    const size_t k =
      Min(size_t(std::upper_bound(m_cellsum.begin(),m_cellsum.end(),r)-
                 m_cellsum.begin()),m_cells.size()-1);
    const size_t p = m_cells[k].first, q = m_cells[k].second;
    const double z0  = e0[p]+ran->Get()*(e0[p+1]-e0[p]);
    const double clo = Max(e1[q],Lower(1,e0[p]));
    const double chi = Min(e1[q+1],Upper(1,e0[p+1]));
    const double z1  = clo+ran->Get()*(chi-clo);
    if (!Allowed(z0,z1)) continue;
    const double l =
      alpha[0]*std::log(z0)+beta[0]*std::log1p(-z0)-c[0]/z0-m_lmax[0][p]+
      alpha[1]*std::log(z1)+beta[1]*std::log1p(-z1)-c[1]/z1-m_lmax[1][q];
    if (std::log(ran->Get())<l) { z[0] = z0; z[1] = z1; return true; }
  }
  m_giveup = true;
  return false;
}
