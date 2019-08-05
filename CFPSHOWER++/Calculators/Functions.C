#include "ATOOLS/Math/Vector.H"
#include <vector>
#include <algorithm>


namespace CFPSHOWER {
  class Functions {
  protected:
    double m_s, m_sijk;
    std::vecor<double>              m_zi;
    std::vecor<std::vecor<double> > m_sij, m_tij;
    
    void InitInvariants(std::vector<ATOOLS::Vec4D> & moms);
    
    inline const double & Pqq(const double & z) const { return (1.+sqr(z))/(1.-z); }
    inline const double & P_q2qpqpq(const size_t & i, const size_t & j) const {
      if (i==j || i<0 || i>2 || j<0 || j>2) return 0.;
      return 1./2. * m_sijk/m_sij[i][j] *
	( sqr(m_tij[i][j])/(m_sij[i][j]*sijk) +
	  (4.*m_zi[2]+sqr(m_zi[i]-m_zi[j]))/(m_zi[i]+m_zi[j]) +
	  (m_zi[i]+m_zi[j]-m_sij[i][j]/m_sijk)  );
    }
    inline const double & Soft_q2qpqpq(const size_t & i, const size_t & j) const {
      if (i==j || i<0 || i>2 || j<0 || j>2) return 0.;
      return 1./2. *
	( -sqr(m_tij[i][j]/m_sij[i][j] - (m_zi[i]-m_zi[j])/(m_zi[i]+m_zi[j]) ) -
	  4.*m_z[2]/(m_zi[i]+m_zi[j]) * (1.-sijk/m_sij[i][j]) );
    }
    inline const double & Id_q2qqq(const size_t & i, const size_t & j) const {
      if (i==j || i<0 || i>2 || j<0 || j>2) return 0.;
      size_t k = 3-i-j;
      return 1./m_Nc *
	( 2.m_sij[j][k]/m_sij[i][j] -
	  sqr(m_sijk)/(m_sij[i][j]*m_sij[i][k]) * m_zi[i]/2. *
	  ( 1.+sqr(m_zi[i])/((1.-m_zi[j])*(1.-m_zi[k])) -
	    (1.+2.*(1.-m_zi[j])/(1.-m_zi[k]) ) ) +
	  m_sijk/m_sij[i][j] *
	  ( Pqq(m_zi[j]) - 2.*m_zi[j]/(1.-m_zi[k]) ) );
    }
    
  public:
    inline const double &
    P_q2qp(std::vector<ATOOLS::Vec4D> & moms) const {
      if (!InitInvariants(moms)) return 0.;
      return P_q2qpqpq(0,1);
    }
    
    inline const double &
    Soft_q2qpqpq(std::vector<ATOOLS::Vec4D> & moms) const {
      if (!InitInvariants(moms)) return 0.;
      return Soft_q2qpqpq(0,1);
    }
    
    inline const double &    
    Id_q2qqq(std::vector<ATOOLS::Vec4D> & moms) const {
      if (!InitInvariants(moms)) return 0.;
      return Soft_q2qpqpq(0,1);
    }
  };

  bool Functions::InitInvariants(std::vector<ATOOLS::Vec4D> & moms) {
    if (moms.size()<5) return false;
    m_zi.resize(moms.size()-3);
    m_sij.resize(moms.size()-3);
    double m_s = (moms[0]+moms[1]).Abs2();
    for (size_t i=2;i<moms.size()-1;i++) {
      m_zi[i-2] = moms[i]*moms[1]/m_s;
      m_sij[i].resize(moms.size()-3);
      for (size_t j=i;j<moms.size()-1;j++) {
	if (i==j) continue;
	m_sij[i][j] = m_s[j][i] = (moms[i]+moms[j]).Abs2();
      }
    }
    if (moms.size()<6) return true;
    m_sijk = (moms[2]+moms[3]+moms[4]).Abs2();
    m_tij.resize(moms.size()-3);
    for (size_t i=2;i<moms.size()-1;i++) {
      m_tij[i].resize(moms.size()-3);
      for (size_t j=i;j<moms.size()-1;j++) {
	if (i==j) continue;
	size_t k     = 3-i-j;
	m_tij[i][j]  = m_tij[j][i] =
	  2.*(m_z[i]*m_sij[j][k]-m_z[j]*m_sij[i][k])/(m_zi[i]+m_zi[j]);
	double diff  = (m_zi[i]-m_zi[j])/(m_zi[i]+m_zi[j])*m_sij[i][j];
	m_tij[i][j] += diff;
	m_tij[j][i] -= diff;
      }
    }
  };
}
