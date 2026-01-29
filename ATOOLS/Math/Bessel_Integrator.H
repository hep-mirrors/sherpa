#ifndef ATOOLS_Math_Bessel_Integrator_H
#define ATOOLS_Math_Bessel_Integrator_H

#include "ATOOLS/Math/Function_Base.H"
#include "ATOOLS/Math/Special_Functions.H"
#include <vector>

namespace ATOOLS {
  class Bessel_Integrator {
    class Kernel final: public ATOOLS::Function_Base {
    private:

      ATOOLS::Function_Base* p_func;
      size_t m_order;

    public:

      explicit Kernel(ATOOLS::Function_Base* func, const size_t& order = 1)
          : p_func(func), m_order(order)
      {
      }
      double operator()(double x) override
      {
        return (*p_func)(x)*SF.Jn(m_order, x);
      };
    };

  private:

    Kernel m_kernel;
    Gauss_Integrator m_gauss;
    size_t m_order, m_maxbins, m_depth, m_iterator;
    std::vector<double> m_zeroes, m_extrema, m_F, m_Psi;
    std::vector<std::vector<double>> m_M, m_N;
    void Test();
    void FixBins(const bool& output = false);
    double AdaptiveIntegrate(double a, double b, double tol);

  public:

    explicit Bessel_Integrator(Function_Base* f, const size_t& order = 1);

    double operator()(double tolerance = 1.e-4, bool output = false);
  };
} // namespace ATOOLS

#endif
