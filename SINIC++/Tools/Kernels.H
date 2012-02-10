#ifndef SINIC_Tools_Kernels_H
#define SINIC_Tools_Kernels_H

#include "SINIC++/Eikonals/Eikonal_Contributor.H"
#include "SINIC++/Tools/DEQ_Solver.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

namespace SINIC {
  /*!
    \class DEQ_Kernel_NoKT
    \brief The class for the \f$k_\perp\f$-independent realisation of 
    differential equations defining the single channel eikonals.  Due 
    to inheritance from SINIC::DEQ_Kernel_Base, again, the main method 
    here is the operator,  
    DEQ_Kernel_NoKT::operator()(const std::vector<double> &input,const double param=0.).
  */
  class DEQ_Kernel_NoKT : public DEQ_Kernel_Base {
  private:
    double              m_lambda, m_Delta, m_expfactor;
    absorption::code    m_absorp;
    std::vector<double> m_output;
  public:
    DEQ_Kernel_NoKT(const double & lambda,const double & Delta,
		    const absorption::code & absorp) :
      m_lambda(lambda), m_Delta(Delta), m_expfactor(1./2.),
      m_absorp(absorp)
    { 
      m_output = std::vector<double>(2,0.);
    }
    /*!
      \fn std::vector<double> & DEQ_Kernel_NoKT::operator()(const std::vector<double> &input,const double param=0.)
      \brief The representation of \f$\vec f(\vec x,\,y)\f$ for the case of 
      the \f$k_\perp\f$-independent definition of the single channel eikonal.

      The corresponding set of differential equations reads
      \f[
      \frac{\mbox{\rm d}\Omega_{i(k)}}{\mbox{\rm d}y} = 
      +\exp\left[-\lambda\cdot\frac{\Omega_{i(k)}+\Omega_{(i)k}}{2}\right]\cdot\Delta\cdot\Omega_{i(k)}
      \f]
      \f[
      \frac{\mbox{\rm d}\Omega_{(i)k}}{\mbox{\rm d}y} =
      -\exp\left[-\lambda\cdot\frac{\Omega_{i(k)}+\Omega_{(i)k}}{2}\right]\cdot\Delta\cdot\Omega_{(i)k}\,,
      \f]
      where the dependence of the \f$\Omega_{i(k)}\f$ and \f$\Omega_{i(k)}\f$ on impact parameters
      \f$\vec b_\perp^{(1,2)}\f$ and on the rapidity \f$y\f$ has been suppressed.  In the actual
      implementation, the following identification has been made:
      \f[
      \Omega_{i(k)} \to x_1\;,\;\;\;
      \Omega_{(i)k} \to x_2\;,\;\;\;\mbox{\rm and}\;\;\;
      y \to y\,.
      \f]
      Therefore
      \f[
      \left(\begin{array}{c}
      +\exp\left[-\lambda\cdot\frac{\Omega_{i(k)}+\Omega_{(i)k}}{2}\right]\cdot\Delta\\
      -\exp\left[-\lambda\cdot\frac{\Omega_{i(k)}+\Omega_{(i)k}}{2}\right]\cdot\Delta
      \end{array}\right) 
      \to \vec f(\vec x,\,y)\,.
      \f]
    */
    std::vector<double> & operator()(const std::vector<double> & input,
				     const double param=0.);
  };


  class Integration_Kernel_Theta : public ATOOLS::Function_Base {
  private:
    Eikonal_Contributor * p_Omega1, * p_Omega2;
    int                   m_test;
    double                m_errmax12, m_errmax21, m_maxvalue;
    double                m_b,m_b1,m_y;
  public:
    Integration_Kernel_Theta(Eikonal_Contributor * Omega1, 
			     Eikonal_Contributor * Omega2, 
			     const int & test=0); 
    inline void    SetYref(const double & y) { m_y  = y;  }
    inline void    SetB(const double & b)    { m_b  = b;  }
    inline void    Setb1(const double & b1)  { m_b1 = b1; }
    double         operator()(double theta1);

    inline void    ResetMax() { m_maxvalue = 0.; }
    const double & Max() const { return m_maxvalue; }
    void        PrintErrors();
  };

  class Integration_Kernel_B2 : public ATOOLS::Function_Base {
  private:
    Integration_Kernel_Theta m_kernel;
    ATOOLS::Gauss_Integrator m_integrator;
    double                   m_accu;
    int                      m_test;
    double                   m_b,m_y;
  public:
    Integration_Kernel_B2(Eikonal_Contributor * Omega1, 
			  Eikonal_Contributor * Omega2, 
			  const int & test=0);

    ~Integration_Kernel_B2() {}

    void SetYref(const double & y) { 
      m_y = y; 
      m_kernel.SetYref(m_y); 
    }
    void SetB(const double & b) { 
      m_b = b; 
      m_kernel.SetB(m_b); 
    }
    double operator()(double b1);

    inline void    ResetMax()  { m_kernel.ResetMax(); }
    const double & Max() const { return m_kernel.Max(); }

    void   PrintErrors() { m_kernel.PrintErrors(); }
  };
}

#endif

