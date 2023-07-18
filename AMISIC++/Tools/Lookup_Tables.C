#include "AMISIC++/Tools/Lookup_Tables.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"

using namespace AMISIC;
using namespace std;

axis::axis(const size_t & nbins,const double & xmin,const double & xmax,
	   const axis_mode::code & mode) :
  m_nbins(nbins), m_xmin(xmin), m_xmax(xmax), m_mode(mode)
{
  if (m_mode==axis_mode::linear) {
    m_xstep = (m_xmax-m_xmin)/double(m_nbins-1);
  }
  else if (m_mode==axis_mode::log) {
    m_xstep = log(m_xmax/m_xmin)/double(m_nbins);
  }
}

const double axis::x(const size_t & bin) const {
  if (m_nbins==1) return m_mode==axis_mode::linear ? (m_xmin+m_xmax)/2. : sqrt(m_xmin*m_xmax);
  if (bin<m_nbins) {
    if (m_mode==axis_mode::linear)   return m_xmin + bin*m_xstep;
    else if (m_mode==axis_mode::log) return m_xmin * exp(m_xstep*bin);
  }
  return 0.;
}

const size_t axis::bin(const double & x) const {
  if (x>=m_xmin && x<=m_xmax) {
    if (m_mode==axis_mode::linear)   return size_t((x-m_xmin)/m_xstep);
    else if (m_mode==axis_mode::log) return size_t(log(x/m_xmin)/m_xstep);
  }
  if (x<=m_xmin) return 0;
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
  if (xbin>=0 && xbin<m_x.m_nbins) m_values[xbin] = value;
}

const double OneDim_Table::operator()(const double & x) const {
  if (m_x.m_nbins==1)                 return m_values[0];
  if (x>=m_x.m_xmax || x<=m_x.m_xmin) return 0.;
  size_t bin = m_x.bin(x);
  double x1  = m_x.x(bin),    x2 = m_x.x(bin+1);
  double y1  = m_values[bin], y2 = m_values[bin+1];
  return ( y1*(x2-x) + y2*(x-x1) ) / (x2-x1);
}


//////////////////////////////////////////////////////////////////////////////
// Two-dimensional look-up table
//////////////////////////////////////////////////////////////////////////////
TwoDim_Table::TwoDim_Table(const axis & xbins,const axis & ybins) :
  m_x(xbins), m_y(ybins)
{
  m_values.resize(m_x.m_nbins);
  for (size_t i=0;i<m_values.size();i++) m_values[i].resize(m_y.m_nbins, 0.);
}

void TwoDim_Table::Fill(const size_t & xbin,const size_t & ybin,const double & value) {
  if (xbin>=0 && xbin<m_x.m_nbins && ybin>=0 && xbin<m_y.m_nbins)
    m_values[xbin][ybin] = value;
}

const double TwoDim_Table::operator()(const double & x,const double & y) const {
  if (m_x.m_nbins==1) {
    if (y<m_y.m_xmin || y>=m_y.m_xmax) return 0.;
    if (m_y.m_nbins==1) return m_values[0][0];
    size_t ybin = m_y.bin(y);
    double y1   = m_y.x(ybin),       y2 = m_y.x(ybin+1);
    double z1   = m_values[0][ybin], z2 = m_values[0][ybin+1];
    return ( z1*(y2-y) + z2*(y-y1) ) / (y2-y1);
  }
  if (m_y.m_nbins==1) {
    if (x<m_x.m_xmin || x>=m_x.m_xmax) return 0.;
    size_t xbin = m_x.bin(x);
    double x1   = m_x.x(xbin),       x2 = m_x.x(xbin+1);
    double z1   = m_values[xbin][0], z2 = m_values[xbin+1][0];
    return ( z1*(x2-x) + z2*(x-x1) ) / (x2-x1);
  }
  if (x<m_x.m_xmin || x>=m_x.m_xmax || y<m_y.m_xmin || y>=m_y.m_xmax) return 0.;
  size_t xbin = m_x.bin(x),             ybin = m_y.bin(y);
  if (xbin+1>=m_x.m_nbins || ybin+1>=m_y.m_nbins)                     return 0.;
  double x1   = m_x.x(xbin),            x2   = m_x.x(xbin+1);
  double y1   = m_y.x(ybin),            y2   = m_y.x(ybin+1);
  double z11  = m_values[xbin][ybin],   z21  = m_values[xbin+1][ybin];
  double z12  = m_values[xbin][ybin+1], z22  = m_values[xbin+1][ybin+1];
  return ( ((z11*(x2-x)+z21*(x-x1)) * (y2-y) + (z12*(x2-x)+z22*(x-x1)) * (y-y1)) /
	   ((x2-x1)*(y2-y1)) );
}



