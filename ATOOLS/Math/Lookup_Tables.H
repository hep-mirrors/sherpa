#ifndef ATOOLS_Math_Lookup_Tables_H
#define ATOOLS_Math_Lookup_Tables_H

#include <cstddef>
#include <fstream>
#include <vector>
#include <memory> // For std::unique_ptr

namespace ATOOLS {
  struct axis_mode {
    enum code { linear = 1, log = 2 };
  };

  class axis {
  public:
    size_t          m_nbins;
    double          m_xmin, m_xmax, m_xstep;
    axis_mode::code m_mode;
    axis(size_t nbins, double xmin, double xmax,
         axis_mode::code mode = axis_mode::linear);
    axis(const axis& old) = default;
    ~axis() = default;
    double x(size_t bin) const;
    size_t bin(double x) const;
  };

  class OneDim_Table {
  private:
    axis                m_x;
    std::vector<double> m_values;

  public:
    explicit OneDim_Table(const axis& xbins);
    void          Fill(size_t xbin, double value);
    void          Rescale(double factor);
    double        operator()(double x) const;
    double        Inverse(double value) const;
    double        Integral() const;
    std::unique_ptr<OneDim_Table> Cumulative(double expo, double& integral) const;
    void          OutputToCSV(std::ofstream& outfile) const;

    inline double Value(size_t xbin) const { return m_values[xbin]; }
    inline size_t NBins() const { return m_x.m_nbins; }
    inline const axis& GetAxis() const { return m_x; }
  };

  class TwoDim_Table {
  private:
    axis                m_x, m_y;
    std::vector<double> m_values;

    inline size_t Index(size_t xbin, size_t ybin) const {
      return xbin * (m_y.m_nbins + 1) + ybin;
    }

  public:
    TwoDim_Table(const axis& xbins, const axis& ybins);
    void   Fill(size_t xbin, size_t ybin, double value);
    double operator()(double x, double y) const;
    double Integral(size_t naxis, size_t bin) const;
    std::unique_ptr<OneDim_Table> Cumulative(double expo, size_t naxis,
                                             size_t bin, double& integral) const;
    void          OutputToCSV(std::ofstream& outfile) const;

    double Value(size_t xbin, size_t ybin) const {
      return m_values[Index(xbin, ybin)];
    }
    inline const axis& GetAxis(size_t dim) const {
      return dim == 0 ? m_x : m_y;
    }
  };

  class ThreeDim_Table {
  private:
    axis                m_x, m_y, m_z;
    std::vector<double> m_values;

    inline size_t Index(size_t xbin, size_t ybin, size_t zbin) const {
      return (xbin * (m_y.m_nbins + 1) + ybin) * (m_z.m_nbins + 1) + zbin;
    }

  public:
    ThreeDim_Table(const axis& xbins, const axis& ybins, const axis& zbins);
    void   Fill(size_t xbin, size_t ybin, size_t zbin, double value);
    double operator()(double x, double y, double z) const;
    double Value(size_t xbin, size_t ybin, size_t zbin) const {
      return m_values[Index(xbin, ybin, zbin)];
    }
    void OutputToCSV(std::ofstream& outfile) const;
  };
}
#endif
