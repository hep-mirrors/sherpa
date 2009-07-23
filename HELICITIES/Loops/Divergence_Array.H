#ifndef HELICITIES_Loops_Divergence_Array_H
#define HELICITIES_Loops_Divergence_Array_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Exception.H"
#include <vector>

using namespace std;
using namespace ATOOLS;

namespace HELICITIES {

  template<typename T>
  class Divergence_Array {
  private:
    std::vector<T> m_result;

  public:
    Divergence_Array() : m_result(6) {
    }
    Divergence_Array(const T& t) : m_result(6, t) {
    }
    Divergence_Array(const T& uv, const T& ir, const T& ir2,
                     const T& fin, const T& eps, const T& eps2) {
      m_result.reserve(6);
      m_result.push_back(uv);
      m_result.push_back(ir);
      m_result.push_back(ir2);
      m_result.push_back(fin);
      m_result.push_back(eps);
      m_result.push_back(eps2);
    }
    Divergence_Array(const std::vector<T>& res) : m_result(res) {
      if (m_result.size()!=6) THROW(fatal_error, "Internal error.");
    }
    ~Divergence_Array() {}

    // extraction of entries
    inline const T&   GetUV()       const { return m_result[0]; }
    inline const T&   GetIR()       const { return m_result[1]; }
    inline const T&   GetIR2()      const { return m_result[2]; }
    inline const T&   GetFinite()   const { return m_result[3]; }
    inline const T&   GetEpsilon()  const { return m_result[4]; }
    inline const T&   GetEpsilon2() const { return m_result[5]; }
    // extraction of entries
    inline T&   UV()       { return m_result[0]; }
    inline T&   IR()       { return m_result[1]; }
    inline T&   IR2()      { return m_result[2]; }
    inline T&   Finite()   { return m_result[3]; }
    inline T&   Epsilon()  { return m_result[4]; }
    inline T&   Epsilon2() { return m_result[5]; }

    inline const std::vector<T>& GetResult() const { return m_result; }

    inline T&       operator[] (unsigned short int i)        {
      return m_result[i];
    }

    inline const T& operator[] (unsigned short int i) const  {
      return m_result[i];
    }

    // sign flip
    inline Divergence_Array<T> operator- () const {
      return Divergence_Array<T>(-m_result[0],-m_result[1],-m_result[2],
                                 -m_result[3],-m_result[4],-m_result[5]);
    }

    // addition/subtraction 
    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator+ (const Divergence_Array<T1>& da) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0]+da[0],
                                             m_result[1]+da[1],
                                             m_result[2]+da[2],
                                             m_result[3]+da[3],
                                             m_result[4]+da[4],
                                             m_result[5]+da[5]);
    }

    inline Divergence_Array<T>&
    operator+=(const Divergence_Array<T>& da) {
      m_result[0]+=da[0];
      m_result[1]+=da[1];
      m_result[2]+=da[2];
      m_result[3]+=da[3];
      m_result[4]+=da[4];
      m_result[5]+=da[5];
      return *this;
    }


    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator- (const Divergence_Array<T1>& da) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0]-da[0],
                                             m_result[1]-da[1],
                                             m_result[2]-da[2],
                                             m_result[3]-da[3],
                                             m_result[4]-da[4],
                                             m_result[5]-da[5]);
    }

    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator+ (const T1& c) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0],
                                             m_result[1],
                                             m_result[2],
                                             m_result[3]+c,
                                             m_result[4],
                                             m_result[5]);
    }

    inline Divergence_Array<T>&
    operator+=(const T& c) {
      m_result[3]+=c;
      return *this;
    }

    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator- (const T1& c) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0],
                                             m_result[1],
                                             m_result[2],
                                             m_result[3]-c,
                                             m_result[4],
                                             m_result[5]);
    }

    // multiplication
    // keep in mind:
    //   - terms of order 1/epsUV^2, 1/epsIR^3, 1/epsIR^4, eps^2 are not calculated
    //   - also overlapping IR/UV divergences are not calculated yet,
    //     -> their treatment needs to be thought about
    //   - terms 1/epsIR^2 * epsUV, 1/epsIR^2 * epsIR, etc. are not differentiated properly
    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator* (const Divergence_Array<T1>& da) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0]*da[3]+da[0]*m_result[3],
                                             m_result[1]*da[3]+da[1]*m_result[3]
                                            +m_result[2]*da[4]+da[2]*m_result[4],
                                             m_result[2]*da[3]+da[2]*m_result[3]
                                            +m_result[1]*da[1]+da[1]*m_result[1],
                                             m_result[3]*da[3]+da[3]*m_result[3]
                                            +m_result[0]*da[4]+da[0]*m_result[4]
                                            +m_result[1]*da[4]+da[1]*m_result[4]
                                            +m_result[2]*da[5]+da[2]*m_result[5],
                                             m_result[4]*da[3]+da[4]*m_result[3]
                                            +m_result[0]*da[5]+da[0]*m_result[5]
                                            +m_result[1]*da[5]+da[1]*m_result[5],
                                             m_result[5]*da[3]+da[5]*m_result[3]);
    }

    // multiplication with scalar
    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator* (const T1& s) const {
      return Divergence_Array<PROMOTE(T,T1)>(s*m_result[0],
                                             s*m_result[1],
                                             s*m_result[2],
                                             s*m_result[3],
                                             s*m_result[4],
                                             s*m_result[5]);
    }

    inline Divergence_Array<T>&
    operator*=(const T& s) {
      m_result[0]*=s;
      m_result[1]*=s;
      m_result[2]*=s;
      m_result[3]*=s;
      m_result[4]*=s;
      m_result[5]*=s;
      return *this;
    }

    template<typename T1>
    inline Divergence_Array<PROMOTE(T,T1)>
    operator/ (const T1& s) const {
      return Divergence_Array<PROMOTE(T,T1)>(m_result[0]/s,
                                             m_result[1]/s,
                                             m_result[2]/s,
                                             m_result[3]/s,
                                             m_result[4]/s,
                                             m_result[5]/s);
    }

  }; // end of class Divergence_Array

  // typedefs
  typedef Divergence_Array<double>  DivArrD;
  typedef Divergence_Array<Complex> DivArrC;

  template<typename T1, typename T2>
  inline Divergence_Array<PROMOTE(T1,T2)>
  operator+ (const T1& s, const Divergence_Array<T2>& a) {
    return a+s;
  }

  template<typename T1, typename T2>
  inline Divergence_Array<PROMOTE(T1,T2)>
  operator- (const T1& s, const Divergence_Array<T2>& a) {
    return -a+s;
  }

  template<typename T1, typename T2>
  inline Divergence_Array<PROMOTE(T1,T2)>
  operator* (const T1& s, const Divergence_Array<T2>& a) {
    return a*s;
  }

  template<typename T1, typename T2>
  inline Divergence_Array<PROMOTE(T1,T2)>
  operator/ (const T1& s, const Divergence_Array<T2>& a) {
    // only implemented for DivArrs with zero 1/eps parts
    // 1/(a+b*eps+c*eps^2)
    // = 1/a - b/a^2 * eps + (2b^2-ca)/(2a^3) * eps^2
    return s*Divergence_Array<PROMOTE(T1,T2)>(0.,
                                              0.,
                                              0.,
                                              1./a[3],
                                              -a[4]/sqr(a[3]),
                                              (2.*sqr(a[4])-a[5]*a[3])
                                                /(2.*a[3]*a[3]*a[3]));
  }

  template<typename T1, typename T2>
  inline Divergence_Array<PROMOTE(T1,T2)>
  operator/ (const Divergence_Array<T1>& a1, const Divergence_Array<T2>& a2) {
    return a1*(1./a2);
  }

  inline DivArrC conj(DivArrC divarr) {
    return DivArrC(std::conj(divarr[0]),
                   std::conj(divarr[1]),
                   std::conj(divarr[2]),
                   std::conj(divarr[3]),
                   std::conj(divarr[4]),
                   std::conj(divarr[5]));
  }

  inline DivArrD real(DivArrC divarr) {
    return DivArrD(std::real(divarr[0]),
                   std::real(divarr[1]),
                   std::real(divarr[2]),
                   std::real(divarr[3]),
                   std::real(divarr[4]),
                   std::real(divarr[5]));
  }

  inline DivArrD imag(DivArrC divarr) {
    return DivArrD(std::imag(divarr[0]),
                   std::imag(divarr[1]),
                   std::imag(divarr[2]),
                   std::imag(divarr[3]),
                   std::imag(divarr[4]),
                   std::imag(divarr[5]));
  }
  
} // end of namespace HELICITIES

#define RED(ARG) om::red<<ARG<<om::reset
#define GREEN(ARG) om::green<<ARG<<om::reset
#define BLUE(ARG) om::blue<<ARG<<om::reset
#define YELLOW(ARG) om::brown<<ARG<<om::reset
#define BLACK(ARG) ARG
#define BOLD(ARG) om::bold<<ARG<<om::reset

namespace ATOOLS {

  template<typename T>
  std::ostream& operator<<(std::ostream& s, const HELICITIES::Divergence_Array<T>& divarr) {
    return s<<BOLD('(')<<BLUE(divarr[0])<<BOLD(',')<<RED(divarr[1])<<BOLD(',')
		       <<YELLOW(divarr[2])<<BOLD(',')<<BOLD(divarr[3])<<BOLD(',')
		       <<GREEN(divarr[4])<<BOLD(',')<<GREEN(divarr[5])<<BOLD(')');
  }

  template<>
  inline bool IsNan<HELICITIES::DivArrD>(const HELICITIES::DivArrD& x) {
    for (size_t i=0; i<x.GetResult().size(); ++i)
      if (ATOOLS::IsNan(x[i])) return true;
    return false;
  }

  template<>
  inline bool IsNan<HELICITIES::DivArrC>(const HELICITIES::DivArrC& x) {
    for (size_t i=0; i<x.GetResult().size(); ++i)
      if (ATOOLS::IsNan(x[i])) return true;
    return false;
  }

  template<>
  inline bool IsZero<HELICITIES::DivArrD>(const HELICITIES::DivArrD& x) {
    for (size_t i=0; i<x.GetResult().size(); ++i)
      if (!ATOOLS::IsZero(x[i])) return false;
    return true;
  }

  template<>
  inline bool IsZero<HELICITIES::DivArrC>(const HELICITIES::DivArrC& x) {
    for (size_t i=0; i<x.GetResult().size(); ++i)
      if (!ATOOLS::IsZero(x[i])) return false;
    return true;
  }
}


#endif
