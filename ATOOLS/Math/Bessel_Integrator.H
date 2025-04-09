#ifndef ATOOLS_Math_Bessel_Integrator_H
#define ATOOLS_Math_Bessel_Integrator_H

#include "ATOOLS/Math/Function_Base.H"
#include <vector>

namespace ATOOLS {
  class Bessel_Integrator {
    class Kernel final : public ATOOLS::Function_Base {
    private:
      ATOOLS::Function_Base* p_func;
      size_t                 m_order;

    public:
      explicit Kernel(ATOOLS::Function_Base* func, const size_t& order = 1)
          : p_func(func), m_order(order)
      {}
      void   SetFunc(ATOOLS::Function_Base* func) { p_func = func; }
      void   SetOrder(const size_t& order) { m_order = order; }
      double operator()(double x) override;
    };
    class TestFunction final : public ATOOLS::Function_Base {
    private:
      size_t m_mode;

    public:
      explicit TestFunction(const size_t& mode = 0) : m_mode(mode) {}
      double operator()(double x) override;
    };

  private:
    Kernel              m_kernel;
    size_t              m_order, m_maxbins, m_depth, m_maxdepth, m_iterator;
    std::vector<double> m_zeroes, m_extrema, m_F, m_Psi;
    std::vector<std::vector<double>> m_M, m_N;
    void                             Test();
    bool                             FillBins(const bool& output = true);
    void                             FixBins(const bool& output = false);

  public:
    explicit Bessel_Integrator(ATOOLS::Function_Base* f,
                               const size_t&          order = 1);

    double operator()();

    inline void SetF(ATOOLS::Function_Base* f) { m_kernel.SetFunc(f); }
    inline void SetOrder(const size_t& order) { m_order = order; }
    inline void SetMaxBins(const size_t& maxbins) { m_maxbins = maxbins; }
    inline void SetDepth(const size_t& depth) { m_depth = depth; }
    inline void SetIterator(const size_t& iterator) { m_iterator = iterator; }
  };
}

#endif
