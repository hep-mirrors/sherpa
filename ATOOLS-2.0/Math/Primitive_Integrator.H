#ifndef Primitive_Integrator_H
#define Primitive_Integrator_H

#include "MathTools.H"
#include <vector>
#include <map>

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
    
    double m_alpha, m_oldalpha, m_weight;
    double m_sum, m_sum2, m_max, m_np;

    std::vector<double>             m_this;
    std::vector<Primitive_Channel*> m_next;

    size_t m_pos;

    friend std::ostream &
    operator<<(std::ostream &str,const Primitive_Channel &channel);

    Primitive_Channel();
    
  public:
    
    // constructor
    Primitive_Channel(Primitive_Channel *const prev,const size_t i,
		      const double &pos=std::numeric_limits<double>::max());
    
    // destructor
    ~Primitive_Channel();

    // member functions
    double Point(const Primitive_Integrand *function,
		 std::vector<double> &point);
    double Point(const Primitive_Integrand *function,
		 std::vector<double> &point,std::vector<double> &opt);

    bool Find(const std::vector<double> &point) const;
    
    void Reset();
    void SetWeight();
    void SetAlpha(const double &alpha);

    bool WriteOut(std::fstream *const file,
		  std::map<Primitive_Channel*,size_t> &pmap) const;
    bool ReadIn(std::fstream *const file,
		std::map<size_t,Primitive_Channel*> &pmap);

    static void CreateRoot(const std::vector<double> &min,
			   const std::vector<double> &max,
			   std::vector<Primitive_Channel*> &channels);

    // inline functions
    inline void SetPosition(const size_t pos) { m_pos=pos; }

    inline double Alpha() const    { return m_alpha;    }
    inline double OldAlpha() const { return m_oldalpha; }
    inline double Weight() const   { return m_weight;   }
    inline size_t Position() const { return m_pos;      }

    inline double Sum() const    { return m_sum;  }
    inline double Sum2() const   { return m_sum2; }
    inline double Max() const    { return m_max;  }
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
  public:

    typedef std::pair<size_t,std::pair<size_t,size_t> > Position_Pair;
    typedef std::map<std::string,Position_Pair>         Position_Map;

  private:

    long unsigned int m_nopt, m_nmax, m_nmaxopt, m_nfirst;

    double m_error, m_scale;
    double m_sum, m_sum2, m_max, m_np;

    std::vector<double>   m_rmin, m_rmax, m_point, m_opt;
    std::map<size_t,bool> m_nosplit;

    std::vector<Primitive_Channel*> m_channels;

    const Primitive_Integrand *p_function;

    size_t m_lastdim, m_mode;

    std::string  m_vname, m_uname;
    Position_Map m_reserved;

    void SelectDimension(const size_t pos);

    void Point();
    void Split();
    void Update();

  public:

    // constructor
    Primitive_Integrator();

    // destructor
    ~Primitive_Integrator();

    // member functions
    void SetDimension(const size_t dim);
    void SetMin(const std::vector<double> &min);
    void SetMax(const std::vector<double> &max);

    void   Initialize();
    double Integrate(const Primitive_Integrand *function);

    void   Point(std::vector<double> &x);
    double Weight(const std::vector<double> &x) const;

    bool WriteOut(const std::string &filename) const;
    bool ReadIn(const std::string &filename);

    void Reserve(const std::string &key,
		 const size_t n,const size_t nprev=0);
    void Split(const std::string &key,
	       const size_t nprev,const std::vector<double> &pos);
    const double *const Reserved(const std::string &key) const;

    // inline functions
    inline void SetNOpt(const long unsigned int &nopt) { m_nopt=nopt; }
    inline void SetNMax(const long unsigned int &nmax) { m_nmax=nmax; }
    inline void SetNMaxOpt(const long unsigned int &nmaxopt) 
    { m_nmaxopt=nmaxopt; }
    inline void SetNFirst(const long unsigned int &nfirst) 
    { m_nfirst=nfirst; }

    inline void SetError(const double &error) { m_error=error; }
    inline void SetScale(const double &scale) { m_scale=scale; }
    inline void SetMode(const size_t mode)    { m_mode=mode;   }

    inline void SetVariableName(const std::string &vname) { m_vname=vname; }
    inline void SetUnitName(const std::string &uname)     { m_uname=uname; } 

    inline void SetFunction(const Primitive_Integrand *function)
    { p_function=function; }

    inline double Mean() const     { return m_sum/m_np; }
    inline double Max() const      { return m_max;      }
    inline double Variance() const    
    { return (m_sum2-m_sum*m_sum/m_np)/(m_np-1.0); }
    inline double Sigma() const  
    { return sqrt(Variance()/m_np); }

  };// end of class Primitive_Integrator

}// end of namespace ATOOLS

#endif
