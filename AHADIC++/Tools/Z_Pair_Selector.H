#ifndef AHADIC_Tools_Z_Pair_Selector_H
#define AHADIC_Tools_Z_Pair_Selector_H

#include "ATOOLS/Math/MathTools.H"
#include <utility>
#include <vector>

namespace AHADIC {
  // Light-cone fractions (z0,z1) of a cluster splitting.  The new clusters
  // have m_k^2 = z_k (1-z_{1-k}) Q^2 - kt^2 and must lie above their lightest
  // decay, m_k^2 >= m_dec,k^2 + kt^2.  Together with the kinematic limits
  // zmin_k <= z_k <= zmax_k this is the allowed region
  //   z0 (1-z1) >= a0,  z1 (1-z0) >= a1,  a_k = (m_dec,k^2 + 2 kt^2)/Q^2.
  //
  // Z_Pair_Selector holds this region for one splitting and provides
  //  - the live interval of z_i, i.e. the z_i that leave an allowed z_j, and
  //    the window of z_j that a given z_i leaves, for the sequential draw;
  //  - a draw of (z0,z1) ~ f0(z0) f1(z1) on the region, by accept/reject
  //    against a piecewise constant envelope.
  // f_k = z^alpha (1-z)^beta exp(-c/z) is the fragmentation function of
  // parton k.
  class Z_Pair_Selector {
  private:
    double m_Q2, m_kt2, m_mdec2[2], m_zmin[2], m_zmax[2], m_a[2];
    double m_lo[2], m_hi[2];
    bool   m_live[2], m_giveup;

    // envelope of the pair draw: panels of z0 and z1 with the largest
    // l = log f on each, and the cells, the pairs of panels the region
    // crosses, with their summed envelope mass
    std::vector<double> m_edge[2], m_lmax[2], m_fmax[2];
    std::vector<double> m_lz, m_newedge, m_newlz, m_cellsum;
    std::vector<std::pair<size_t,size_t> > m_cells;

    void Envelope(const int k, const double lo, const double hi,
                  const double alpha, const double beta, const double c);
  public:
    Z_Pair_Selector();

    // c = gamma (kt^2 + M^2)/kT_0^2 with scale = kt^2 + M^2, switched off
    // for |gamma| <= 5e-3
    static double Exponent(const double gamma, const double kt02,
                           const double scale);
    // maximum of log f on [zmin,zmax]: the stationary point inside, else
    // the larger end
    static double Peak(const double alpha, const double beta, const double c,
                       const double zmin, const double zmax);

    // region of a cluster with mass^2 Q2 splitting with kt2 into clusters
    // above mdec2[k] + kt2, with z_k in [zmin[k],zmax[k]]
    void SetRegion(const double Q2, const double kt2, const double * mdec2,
                   const double * zmin, const double * zmax);

    // live interval [Lower(i),Upper(i)] of z_i, empty if !Live(i)
    bool   Live(const int i)  const { return m_live[i]; }
    double Lower(const int i) const { return m_lo[i]; }
    double Upper(const int i) const { return m_hi[i]; }
    // window of z_j left by z_i, i = 1-j
    double Lower(const int j, const double zi) const {
      return ATOOLS::Max(m_zmin[j],m_a[j]/(1.-zi));
    }
    double Upper(const int j, const double zi) const {
      return ATOOLS::Min(m_zmax[j],1.-m_a[1-j]/zi);
    }
    bool Allowed(const double z0, const double z1) const {
      return !(z0*(1.-z1)*m_Q2-m_kt2<m_mdec2[0]+m_kt2 ||
               z1*(1.-z0)*m_Q2-m_kt2<m_mdec2[1]+m_kt2);
    }

    // (z0,z1) ~ f0 f1 on the region, with f_k given by alpha[k], beta[k],
    // c[k]; false if there is no pair, GaveUp() if no trial was accepted
    bool operator()(const double * alpha, const double * beta,
                    const double * c, double * z);
    bool GaveUp() const { return m_giveup; }
  };
}

#endif
