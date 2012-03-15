#ifndef SHRIMPS_Tools_DEQ_Solver_H
#define SHRIMPS_Tools_DEQ_Solver_H

#include <vector>
#include <cstddef>

namespace SHRIMPS {
  /*!
    \class DEQ_Kernel_Base
    \brief The base class for possible kernels of systems of differential equations.
    The main method here is the purely virtual operator,  
    DEQ_Kernel_Base::operator()(const std::vector<double> &input,const double param=0.).
    
    Up to now this is only specified for differential equations of the form
    \f[
    \frac{\mbox{\rm d}\vec x}{\mbox{\rm d}y} = \vec f(\vec x,\,y)\,,
    \f]
    where the kernels actually represent the \f$\vec f(\vec x,\,y)\f$ of the equation above.
  */
  class DEQ_Kernel_Base {
  public:
    DEQ_Kernel_Base() {}
    virtual ~DEQ_Kernel_Base() {}
    /*!
      \fn std::vector<double> & DEQ_Kernel_Base::operator()(const std::vector<double> &input,const double param=0.)
      \brief The purely virtual representation of \f$\vec f(\vec x,\,y)\f$.
    */
    virtual std::vector<double> & operator()(const std::vector<double> &input,
					     const double param=0.)=0;
  };

  struct deqmode {
    /*!
      \struct deqmode
      \brief The mode for solving the equations defining the eikonal.
      Available or foreseen options up to now are Runge-Kutta methods of
      the second or fourth order, and a simple Euler step method.
    */
    enum code {
      RungeKutta4 = 4,
      RungeKutta2 = 2,
      SimpleEuler = 1
    };
  };

  /*!
    \class DEQ_Solver
    \brief A class to solve simple systems of linear differential equations
    with Runge-Kutta methods of order 2 or 4 or with a simple Euler
    step method.  
    
    Specifically, linear differential equations of the type
    \f[
    \frac{\mbox{\rm d}\vec x}{\mbox{\rm d}y} = \vec f(\vec x,\,y)
    \f]
    can be treated.  Note that so far only real values and parameters are
    foreseen.
    
    The key method is DEQ_Solver::SolveSystem(const std::vector<double> & x0,const int & ysteps)
    which will produce a (two-dimensional) grid \f$\vec x(y_i)\f$ at finite 
    values of \f$y\f$.  This grid can be recovered through DEQ_Solver::X().
  */
  class DEQ_Solver {
  private:
    DEQ_Kernel_Base                  * p_kernel;
    size_t                             m_dim;
    std::vector<std::vector<double> >  m_x;
    deqmode::code                      m_deqmode;
    double                             m_ystart, m_yend, m_stepsize;
    int                                m_test;

    /*!
      \fn void DEQ_Solver::SimpleEuler(const std::vector<double> & x0,const int & steps)
      \brief The method to solve the differential equation with a simple Euler
      step algorithm.

      Using \f$y_{i+1}=y_i+\Delta y\f$ with \f$\Delta y\f$ the stepsize 
      of the algorithm and \f$\vec x_i\equiv \vec x(y_i)\f$ this algorithm reads:
      \f[
      \vec x(y_{i+1}) = \vec x_i+\Delta y\cdot\vec f(\vec x_i,\,y_i)\,.
      \f]
    */
    void SimpleEuler(const int & steps);
    /*!
      \fn void DEQ_Solver::RungeKutta2(const int & steps)
      \brief The method to solve the differential equation with a Runge-Kutta
      algorithm of order 2.

      Using \f$y_{i+1}=y_i+\Delta y\f$ with \f$\Delta y\f$ the stepsize 
      of the algorithm and \f$\vec x_i\equiv \vec x(y_i)\f$ this algorithm makes use
      of some auxiliary point \f$\vec x^{(1)}\f$ and \f$y^{(1)}\f$ given by
      \f[
      \vec x^{(1)}_i = \vec x_i+\frac{\Delta y}{2}\cdot\vec f(\vec x_i,\,y_i)
      \;\;\;\mbox{\rm and}\;\;\;
      y^{(1)}_i = y_i+\frac{\Delta y}{2}\,.
      \f]
      With these two auxiliaries, the new point reads:
      \f[
      \vec x_{i+1} = \vec x_i+\Delta y\cdot \vec f(\vec x_i^{(1)},\,y_i^{(1)})\,. 
      \f]
    */
    void RungeKutta2(const int & steps);
    /*!
      \fn void DEQ_Solver::RungeKutta4(const int & steps)
      \brief The method to solve the differential equation with a Runge-Kutta
      algorithm of order 4.

      Using \f$y_{i+1}=y_i+\Delta y\f$ with \f$\Delta y\f$ the stepsize 
      of the algorithm and \f$\vec x_i\equiv \vec x(y_i)\f$ this algorithm makes use
      of some auxiliary points \f$\vec x^{(1,2,3,4)}\f$ and \f$y^{(1,2,3,4)}\f$ given by
      \f[
      \vec x^{(1)}_i = \vec x_i
      \;\;\;\mbox{\rm and}\;\;\;
      y^{(1)}_i = y_i
      \f]
      \f[
      \vec x^{(2)}_i = \vec x_i+\frac{\Delta y}{2}\cdot\vec f(\vec x_i^{(1)},\,y_i^{(1)})
      \;\;\;\mbox{\rm and}\;\;\;
      y^{(2)}_i = y_i+\frac{\Delta y}{2}\,.
      \f]
      \f[
      \vec x^{(3)}_i = \vec x_i+\frac{\Delta y}{2}\cdot\vec f(\vec x_i^{(2)},\,y_i^{(2)})
      \;\;\;\mbox{\rm and}\;\;\;
      y^{(3)}_i = y_i+\frac{\Delta y}{2}\,.
      \f]
      \f[
      \vec x^{(4)}_i = \vec x_i+\Delta y\cdot\vec f(\vec x_i^{(3)},\,y_i^{(3)})
      \;\;\;\mbox{\rm and}\;\;\;
      y^{(4)}_i = y_i+\Delta y\,.
      \f]
      With these two auxiliaries, the new point reads:
      \f[
      \vec x_{i+1} = \vec x_i+\frac{\Delta y}{6}\cdot 
      \left[\vec f(\vec x_i^{(1)},\,y_i^{(1)})+2\vec f(\vec x_i^{(2)},\,y_i^{(2)})+
            2\vec f(\vec x_i^{(3)},\,y_i^{(3)})+\vec f(\vec x_i^{(4)},\,y_i^{(4)})\right]\,. 
      \f]
    */
    void RungeKutta4(const int & steps);
  public:
    DEQ_Solver(DEQ_Kernel_Base * kernel,const size_t & dim=2,
	       const deqmode::code & deq=deqmode::RungeKutta4,
	       const int & m_test=0);
    inline ~DEQ_Solver() {}    

   /*!
      \fn void DEQ_Solver::SolveSystem(const int & steps)
      \brief The method to solve the differential equation with one of the methods
      implemented so far.  It will therefore call either
      DEQ_Solver::RungeKutta4(const int & steps), or
      DEQ_Solver::RungeKutta2(const int & steps), or
      DEQ_Solver::SimpleEuler(const int & steps).      
    */
    void SolveSystem(const std::vector<double> & x0,int & ysteps,
		     const double & accu=-1);

    inline void SetInterval(const double & ystart, const double & yend) {
      m_ystart = ystart; m_yend = yend;
    }

    const std::vector<std::vector<double> > & X() const { return m_x; }
  };
}

#endif
