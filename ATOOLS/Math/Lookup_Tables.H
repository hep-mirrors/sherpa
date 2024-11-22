#ifndef ATOOLS_Math_Lookup_Tables_H
#define ATOOLS_Math_Lookup_Tables_H

#include <cstddef>
#include <fstream>
#include <vector>

namespace ATOOLS {
  struct axis_mode {
    enum code { linear = 1, log = 2 };
  };

  class axis {
  public:
    size_t          m_nbins;
    double          m_xmin, m_xmax, m_xstep;
    axis_mode::code m_mode;
    axis(const size_t& nbins, const double& xmin, const double& xmax,
         const axis_mode::code& mode = axis_mode::linear);
    axis(const axis& old)
        : axis(old.m_nbins, old.m_xmin, old.m_xmax, old.m_mode)
    {}
    ~axis() = default;
    double x(const size_t& bin) const;
    size_t bin(const double& x) const;
  };

  class OneDim_Table {
  private:
    axis                m_x;
    std::vector<double> m_values;

  public:
    explicit OneDim_Table(const axis& xbins);
    void          Fill(const size_t& xbin, const double& value);
    void          Rescale(const double& factor);
    double        operator()(const double& x) const;
    const double  Inverse(const double& value) const;
    const double  Integral() const;
    OneDim_Table* Cumulative(const double& expo, double& integral) const;
    void          OutputToCSV(std::ofstream& outfile);

    inline double Value(const size_t& xbin) const { return m_values[xbin]; }
    inline size_t Size() const { return m_values.size(); }
    inline const axis& GetAxis() const { return m_x; }
  };

  class TwoDim_Table {
  private:
    axis                              m_x, m_y;
    std::vector<std::vector<double> > m_values;

  public:
    TwoDim_Table(const axis& xbins, const axis& ybins);
    void   Fill(const size_t& xbin, const size_t& ybin, const double& value);
    double operator()(const double& x, const double& y) const;
    const double  Integral(const size_t& naxis, const size_t& bin) const;
    OneDim_Table* Cumulative(const double& expo, const size_t& naxis,
                             const size_t& bin, double& integral) const;
    void          OutputToCSV(std::ofstream& outfile);

    double Value(const size_t& xbin, const size_t& ybin) const
    {
      return m_values[xbin][ybin];
    }
    inline const axis& GetAxis(const size_t& dim) const
    {
      return dim == 0 ? m_x : m_y;
    }
  };

  class ThreeDim_Table {
  private:
    axis                                            m_x, m_y, m_z;
    std::vector<std::vector<std::vector<double> > > m_values;

  public:
    ThreeDim_Table(const axis& xbins, const axis& ybins, const axis& zbins);
    void   Fill(const size_t& xbin, const size_t& ybin, const size_t& zbin,
                const double& value);
    double operator()(const double& x, const double& y, const double& z) const;
    double Value(const size_t& xbin, const size_t& ybin,
                 const size_t& zbin) const
    {
      return m_values[xbin][ybin][zbin];
    }
    void OutputToCSV(std::ofstream& outfile);
  };
}
#endif
