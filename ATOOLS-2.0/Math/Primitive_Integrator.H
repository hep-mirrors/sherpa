#ifndef Primitive_Integrator_H
#define Primitive_Integrator_H

#include "MathTools.H"
#include <vector>

namespace ATOOLS {

  class Primitive_Integrand {
  public:

    // destructor
    virtual ~Primitive_Integrand();

    // member functions
    virtual double operator()(const std::vector<double> &point) const = 0;

  };// end of class Primitive_Integrand

  class Primitive_Channel {
  private:
    
    double m_alpha, m_weight;
    double m_sum, m_sum2, m_np;
    std::vector<double> m_this;
    std::vector<Primitive_Channel*> m_next;

    void SetWeight();

    friend std::ostream &
    operator<<(std::ostream &str,const Primitive_Channel &channel);

    Primitive_Channel();
    
  public:
    
    // constructor
    Primitive_Channel(const std::vector<Primitive_Channel*> &all,
		      Primitive_Channel *const prev,const size_t i);
    
    // destructor
    ~Primitive_Channel();

    // member functions
    double Point(const Primitive_Integrand *function,
		 std::vector<double> &point);

    void Reset();
    
    static void CreateRoot(const std::vector<double> &min,
			   const std::vector<double> &max,
			   std::vector<Primitive_Channel*> &channels);

    // inline functions
    inline void SetAlpha(const double &alpha) { m_alpha=alpha; }

    inline double Alpha() const  { return m_alpha;  }
    inline double Weight() const { return m_weight; }

    inline double Sum() const    { return m_sum;  }
    inline double Sum2() const   { return m_sum2; }
    inline double Points() const { return m_np;   }

    inline double Mean() const   { return m_sum/m_np; }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

    inline bool Boundary() const { return m_next[0]==NULL; }

  };// end of class Primitive_Channel

  std::ostream &
  operator<<(std::ostream &str,const Primitive_Channel &channel);

  class Primitive_Integrator {
  private:

    long unsigned int m_nopt, m_nmax, m_nmaxopt;

    double m_error;
    double m_sum, m_sum2, m_np;

    std::vector<double> m_min, m_max, m_point;

    std::vector<Primitive_Channel*> m_channels;

    const Primitive_Integrand *p_function;

    size_t m_lastdim, m_mode;

    void Point();
    void Split();
    void Update();

  public:

    // constructor
    Primitive_Integrator();

    // destructor
    ~Primitive_Integrator();

    // member functions
    double Integrate(const Primitive_Integrand *function);

    // inline functions
    inline void SetMin(const std::vector<double> &min) { m_min=min; }
    inline void SetMax(const std::vector<double> &max) { m_max=max; }

    inline void SetNOpt(const long unsigned int &nopt) { m_nopt=nopt; }
    inline void SetNMax(const long unsigned int &nmax) { m_nmax=nmax; }
    inline void SetNMaxOpt(const long unsigned int &nmaxopt) 
    { m_nmaxopt=nmaxopt; }

    inline void SetError(const double &error) { m_error=error; }
    inline void SetMode(const size_t mode)    { m_mode=mode;   }

    inline double Mean() const   { return m_sum/m_np; }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

  };// end of class Primitive_Integrator

}// end of namespace ATOOLS

#endif
