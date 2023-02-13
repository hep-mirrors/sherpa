#ifndef AMISIC_Tools_Lookup_Tables_H
#define AMISIC_Tools_Lookup_Tables_H

#include <cstddef>
#include <vector>

namespace AMISIC {
  struct axis_mode {
    enum code {
      unknown = 0,
      linear  = 1,
      log     = 2
    };
  };
  
  class axis {
  public:
    size_t m_nbins;
    double m_xmin, m_xmax, m_xstep;
    axis_mode::code m_mode;
  public:
    axis(const size_t & nbins,const double & xmin,const double & xmax,
	 const axis_mode::code & mode=axis_mode::linear);
    axis(const axis & old) :
      axis(old.m_nbins, old.m_xmin, old.m_xmax, old.m_mode) {}
    const double x(const size_t & bin) const;
    const size_t bin(const double & x) const;
  };
  
  class OneDim_Table {
  private:
    axis m_x;
    std::vector<double> m_values;
  public:
    OneDim_Table(const axis & xbins);
    void Fill(const size_t & xbin,const double & value);
    const double operator()(const double & x) const;
    inline const double Value(const size_t & xbin) const { return m_values[xbin]; }
    inline const double Size()                     const { return m_values.size(); }
  };

  class TwoDim_Table {
  private:
    axis m_x, m_y;
    std::vector<std::vector<double> > m_values;
  public:
    TwoDim_Table(const axis & xbins,const axis & ybins);
    void Fill(const size_t & xbin,const size_t & ybin,const double & value);
    const double operator()(const double & x,const double & y) const;
    const double Value(const size_t & xbin,const size_t & ybin) const {
      return m_values[xbin][ybin];
    }
  };
}
#endif
