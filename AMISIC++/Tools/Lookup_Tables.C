#include "AMISIC++/Tools/Lookup_Tables.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"

using namespace AMISIC;
using namespace std;

axis::axis(const size_t & nbins,const double & xmin,const double & xmax,
           const axis_mode& mode)
    : m_nbins(nbins), m_xmin(xmin), m_xmax(xmax), m_mode(mode)
{
  if (m_nbins == 1) m_xstep = 1.;
  else if (m_mode==axis_mode::linear) {
    m_xstep = (m_xmax-m_xmin)/double(m_nbins-1);
  }
  else if (m_mode==axis_mode::log) {
    m_xstep = log(m_xmax/m_xmin)/double(m_nbins-1);
  }
}

double axis::x(const size_t & bin) const {
  if (m_nbins==1) return m_mode==axis_mode::linear ? (m_xmin+m_xmax)/2. : sqrt(m_xmin*m_xmax);
  if (bin>=m_nbins) {
    THROW(normal_exit,"Wrong bin called");
  }
  if (m_mode==axis_mode::linear)   return m_xmin + (double)bin*m_xstep;
  else if (m_mode==axis_mode::log) return m_xmin * exp(m_xstep*(double)bin);
  else return 0.;
}

size_t axis::bin(const double & x) const {
  if (x>=m_xmin && x<=m_xmax) {
    if (m_mode == axis_mode::linear)
      return static_cast<size_t>((x - m_xmin) / m_xstep);
    else if (m_mode == axis_mode::log)
      return static_cast<size_t>(log(x / m_xmin) / m_xstep);
  }
  if (x<m_xmin) return 0;
  return m_nbins-1;
}

//////////////////////////////////////////////////////////////////////////////
// One-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
OneDim_Table::OneDim_Table(const axis & xbins) :
  m_x(xbins)
{
  m_values.resize(m_x.m_nbins, 0.);
}

void OneDim_Table::Fill(const size_t & xbin,const double & value) {
  if (xbin<m_x.m_nbins) m_values[xbin] = value;
}

double OneDim_Table::operator()(const double & x) const {
  if (m_x.m_nbins==1)                 return m_values[0];
  if (x>=m_x.m_xmax || x<=m_x.m_xmin) return 0.;
  size_t bin = m_x.bin(x);
  if (bin>=m_x.m_nbins) return 0;
  double x1 = m_x.x(bin), x2 = m_x.x(bin+1);
  return ( m_values[bin]*(x2-x) + m_values[bin+1]*(x-x1) ) / (x2-x1);
}


//////////////////////////////////////////////////////////////////////////////
// Two-dimensional look-up table: will assume y axis has more than one bin.
//////////////////////////////////////////////////////////////////////////////
TwoDim_Table::TwoDim_Table(const axis & xbins,const axis & ybins) :
  m_x(xbins), m_y(ybins)
{
  if (m_y.m_nbins<=1) THROW(fatal_error,"Only one or less bins in y direction.")
  m_values.resize(m_x.m_nbins);
  for (auto& val : m_values) val.resize(m_y.m_nbins, 0.);
}

void TwoDim_Table::Fill(const size_t & xbin,const size_t & ybin,const double & value) {
  if (xbin<m_x.m_nbins && ybin<m_y.m_nbins)
    m_values[xbin][ybin] = value;
}

double TwoDim_Table::operator()(const double & x,const double & y) const {
  if (m_x.m_nbins==1) {
    if (y<m_y.m_xmin || y>=m_y.m_xmax) return 0.;
    if (m_y.m_nbins==1) return m_values[0][0];
    size_t ybin = m_y.bin(y);
    if (ybin>=m_y.m_nbins) return 0.;
    double y1 = m_y.x(ybin), y2 = m_y.x(ybin+1);
    return ( m_values[0][ybin]   * (y2-y) +
	     m_values[0][ybin+1] * (y-y1) ) / (y2-y1);
  }
  if (x<m_x.m_xmin || x>=m_x.m_xmax || y<m_y.m_xmin || y>=m_y.m_xmax) return 0.;
  size_t xbin = m_x.bin(x),             ybin = m_y.bin(y);
  if (xbin>=m_x.m_nbins || ybin>=m_y.m_nbins) return 0.;
  double x1   = m_x.x(xbin),            x2   = m_x.x(xbin+1);
  double y1   = m_y.x(ybin),            y2   = m_y.x(ybin+1);
  return ( ((y2-y) * (m_values[xbin][ybin]     * (x2-x) +
		      m_values[xbin+1][ybin]   * (x-x1)) +
	    (y-y1) * (m_values[xbin][ybin+1]   * (x2-x) +
		      m_values[xbin+1][ybin+1] * (x-x1)) ) /
	   ((x2-x1)*(y2-y1)) );
}



//////////////////////////////////////////////////////////////////////////////
// Three-dimensional look-up table: will assume y axis has more than one bin.
//////////////////////////////////////////////////////////////////////////////
ThreeDim_Table::ThreeDim_Table(const axis & xbins,const axis & ybins,const axis & zbins) :
  m_x(xbins), m_y(ybins), m_z(zbins)
{
  if (m_z.m_nbins<=1) THROW(fatal_error,"Only one or less bins in z direction.")
  m_values.resize(m_x.m_nbins);
  for (auto& valy : m_values) {
    valy.resize(m_y.m_nbins);
    for (auto& valz : valy) valz.resize(m_z.m_nbins,0.);
  }
}

void ThreeDim_Table::Fill(const size_t & xbin,const size_t & ybin,const size_t & zbin,
			  const double & value) {
  if (xbin<m_x.m_nbins && ybin<m_y.m_nbins && zbin<m_z.m_nbins)
    m_values[xbin][ybin][zbin] = value;
}

double ThreeDim_Table::operator()(const double & x,const double & y,const double & z) const {
  if (m_x.m_nbins==1 && m_y.m_nbins==1) {
    if (m_z.m_nbins==1) return m_values[0][0][0];
    if (z<m_z.m_xmin || z>=m_z.m_xmax) return 0.;
    size_t zbin = m_z.bin(y);
    if (zbin>=m_z.m_nbins) return 0.;
    double z1 = m_z.x(zbin), z2 = m_z.x(zbin+1);
    return ( m_values[0][0][zbin]   * (z2-z) +
	     m_values[0][0][zbin+1] * (z-z1) ) / (z2-z1);
  }
  if (m_x.m_nbins==1 && m_y.m_nbins>1) {
    if (y<m_y.m_xmin || y>=m_y.m_xmax ||
	z<m_z.m_xmin || z>=m_z.m_xmax) return 0.;
    size_t ybin = m_y.bin(y), zbin = m_z.bin(z);
    if (ybin>=m_y.m_nbins || zbin>=m_z.m_nbins) return 0.;
    double y1 = m_y.x(ybin), y2 = m_y.x(ybin+1);
    double z1 = m_z.x(zbin), z2 = m_z.x(zbin+1);
    return ( ((z2-z) * (m_values[0][ybin][zbin]     * (y2-y) +
			m_values[0][ybin+1][zbin]   * (y-y1)) +
	      (z-z1) * (m_values[0][ybin][zbin+1]   * (y2-y) +
			m_values[0][ybin+1][zbin+1] * (y-y1)) ) /
	     ((y2-y1)*(z2-z1)) );
  }
  if (m_x.m_nbins>1 && m_y.m_nbins==1) {
    if (x<m_x.m_xmin || x>=m_x.m_xmax ||
	z<m_z.m_xmin || z>=m_z.m_xmax) return 0.;
    size_t xbin = m_x.bin(x), zbin = m_z.bin(z);
    if (xbin>=m_x.m_nbins || zbin>=m_z.m_nbins) return 0.;
    double x1 = m_x.x(xbin), x2 = m_x.x(xbin+1);
    double z1 = m_z.x(zbin), z2 = m_z.x(zbin+1);
    return ( ((z2-z) * (m_values[xbin][0][zbin]     * (x2-x) +
			m_values[xbin+1][0][zbin]   * (x-x1)) +
	      (z-z1) * (m_values[xbin][0][zbin+1]   * (x2-x) +
			m_values[xbin+1][0][zbin+1] * (x-x1)) ) /
	     ((x2-x1)*(z2-z1)) );
  }  
  if (x<m_x.m_xmin || x>=m_x.m_xmax ||
      y<m_y.m_xmin || y>=m_y.m_xmax ||
      z<m_z.m_xmin || z>=m_z.m_xmax) return 0.;
  size_t xbin = m_x.bin(x), ybin = m_y.bin(y), zbin = m_z.bin(z);
  if (xbin>=m_x.m_nbins || ybin>=m_y.m_nbins || zbin>=m_z.m_nbins) return 0.;
  double x1 = m_x.x(xbin), x2   = m_x.x(xbin+1);
  double y1 = m_y.x(ybin), y2   = m_y.x(ybin+1);
  double z1 = m_z.x(zbin), z2   = m_z.x(zbin+1);
  return ( ((z2-z) * ( (y2-y) * ( m_values[xbin][ybin][zbin]       * (x2-x) +
				  m_values[xbin+1][ybin][zbin]     * (x-x1) ) +
		       (y-y1) * ( m_values[xbin][ybin+1][zbin]     * (x2-x) +
				  m_values[xbin+1][ybin+1][zbin]   * (x-x1) ) ) +
	    (z-z1) * ( (y2-y) * ( m_values[xbin][ybin][zbin+1]     * (x2-x) +
				  m_values[xbin+1][ybin][zbin+1]   * (x-x1) ) +
		       (y-y1) * ( m_values[xbin][ybin+1][zbin+1]   * (x2-x) +
				  m_values[xbin+1][ybin+1][zbin+1] * (x-x1) ) ) )/
	   ((x2-x1)*(y2-y1)*(z2-z1)) );
}



