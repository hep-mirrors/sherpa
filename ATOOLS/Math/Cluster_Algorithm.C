#include "ATOOLS/Math/Cluster_Algorithm.H"

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/My_Limits.H"

using namespace ATOOLS;

template <class PointType,class MeasureType,class RecombinationType>
Cluster_Algorithm<PointType,MeasureType,RecombinationType>::
Cluster_Algorithm(): m_recalc(false) {}

template <class PointType,class MeasureType,class RecombinationType>
bool Cluster_Algorithm<PointType,MeasureType,RecombinationType>::
ArrangePoints() 
{
  for (size_t i(0);i<m_n;++i) {
    m_p[i]=m_p[m_i[i]];
    for (size_t j(0);j<m_j.size();++j) if (m_j[j]==(int)m_i[i]) m_j[j]=i;
  }
  m_p.resize(m_n);
  return true;
}

template <class PointType,class MeasureType,class RecombinationType>
bool Cluster_Algorithm<PointType,MeasureType,RecombinationType>::
Cluster(const double &crit,const cs::code &code)
{
  m_n=m_p.size();
  if (m_n>m_i.size()) {
    m_i.resize(m_n);
    m_d.resize(m_n);
    for (size_t i(0);i<m_n;++i) m_d[i]=Double_Vector(m_n);
  }
  m_j.resize(m_n);
  m_dmin=std::numeric_limits<double>::max();
  for (size_t i(0);i<m_n;++i) m_j[i]=m_i[i]=i;
  for (size_t i(0);i<m_n;++i) {
    SetDMin(i,i,m_d[i][i]=m_measure(m_p[i]));
    for (size_t j(i+1);j<m_n;++j) 
      SetDMin(i,j,m_d[i][j]=m_measure(m_p[i],m_p[j]));
  }
  m_r.clear();
  m_lp.clear();
  if (m_n==0) return false;
  if (code==cs::num && m_n<=crit) return ArrangePoints();
  while (true) {
    m_r.push_back(m_dmin);
    size_t iimin(m_i[m_imin]), ijmin(m_i[m_jmin]);
    if (iimin!=ijmin) {
      for (size_t i(0);i<m_j.size();++i) if (m_j[i]==(int)ijmin) m_j[i]=iimin;
      m_p[iimin]=m_recom(m_p[iimin],m_p[ijmin]);
    }
    else {
      m_lp.push_back(m_p[iimin]);
      m_recom(m_p[iimin]);
      int pos(-m_lp.size());
      for (size_t i(0);i<m_j.size();++i) if (m_j[i]==(int)ijmin) m_j[i]=pos;
    }
    for (size_t i(m_jmin+1);i<m_n;++i) m_i[i-1]=m_i[i];
    --m_n;
    if ((code==cs::num && m_n<=crit) || m_n==0) return ArrangePoints();
    if (iimin!=ijmin) {
      if (!m_recalc) {
	for (size_t i(0);i<m_imin;++i) 
	  m_d[m_i[i]][iimin]=m_measure(m_p[m_i[i]],m_p[iimin]);
	m_d[iimin][iimin]=m_measure(m_p[iimin]);
	for (size_t j(m_imin+1);j<m_n;++j) 
	  m_d[iimin][m_i[j]]=m_measure(m_p[iimin],m_p[m_i[j]]);
      } 
      else {
	for (size_t i(0);i<m_n;++i) {
	  for (size_t j(i+1);j<m_n;++j) 
	    m_d[m_i[i]][m_i[j]]=m_measure(m_p[m_i[i]],m_p[m_i[j]]);
	}
      }
    }
    m_dmin=std::numeric_limits<double>::max();
    for (size_t i(0);i<m_n;++i) {
      size_t ii(m_i[i]);
      SetDMin(i,i,m_d[ii][ii]);
      for (size_t j(i+1);j<m_n;++j) SetDMin(i,j,m_d[ii][m_i[j]]);
    }
    if (m_dmin==std::numeric_limits<double>::max()) break;
    if (code==cs::dist && m_dmin>=crit) return ArrangePoints();
  }
  ArrangePoints();
  return false;
}

