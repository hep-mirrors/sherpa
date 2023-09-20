#ifndef METOOLS_Explicit_C_Scalar_H
#define METOOLS_Explicit_C_Scalar_H

#include "METOOLS/Explicit/C_Object.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <vector>
#include <iostream>
 
namespace METOOLS {

  /// @brief Implementation of a scalar current.
  template <class Scalar>
  class CScalar: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  private:
    /// @brief The spin components of the current. In this case, a single complex number.
    SComplex m_x;

    /// @brief The accuracy for all scalars.
    /// @note Not used anywhere
    static double s_accu;

    /// @brief A vector of CScalars that cleans its memory upon destruction.
    ///
    /// If threading is enabled, all the CScalars created with the CScalar<Scalar>::New() method 
    /// will be allocated inside this vector behind the scenes.
    static ATOOLS::AutoDelete_Vector<CScalar> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &,const CScalar<_Scalar> &);

  public:
    /// @name New pointer creator methods
    /// @{    
    /// Create a pointer to a new CScalar object.
    /// If threading is enabled, existing slots inside CScalar<Scalar>::s_objects will be re-used.
    static CScalar *New();
    static CScalar *New(const CScalar &s);    
    /// @}

    CObject* Copy() const;

    /// @brief Delete this object.
    ///
    /// If threading is enabled, the memory it occupies is added to CScalar<Scalar>::s_objects
    /// in order to be re-used in the future.
    void Delete();

    /// @brief Return true if CScalar<Scalar>::m_x is zero.
    bool IsZero() const;

    /// @name Constructors
    /// @{

    /// @brief Default constructor.
    ///
    /// Initialize all data members to zero (m_x, m_c, m_h, m_s)
    inline CScalar(): m_x(Scalar(0.0))
    { 
      m_h=m_s=0; m_c[0]=m_c[1]=0;
    }
    /// @brief Copy constructor
    inline CScalar(const CScalar &s): m_x(s.m_x)
    { 
      m_h=s.m_h; m_s=s.m_s; m_c[0]=s.m_c[0]; m_c[1]=s.m_c[1];
    }
    /// @brief Constructor
    ///
    /// @param x  the spin components
    /// @param cr the color index
    /// @param ca the anti-color index
    /// @param h  the helicity
    /// @param s  the NLO subtraction flag (m_s)
    inline CScalar(const Scalar &x,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0): m_x(x)
    { 
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca;
    }
    /// @brief Constructor
    ///
    /// @param x  the spin components
    /// @param cr the color index
    /// @param ca the anti-color index
    /// @param h  the helicity
    /// @param s  the NLO subtraction flag (m_s)
    inline CScalar(const SComplex &x,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0): m_x(x)
    { 
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca;
    }
    /// @brief Scaled constructor
    ///
    /// @param s another CScalar
    /// @param c a complex number to scale by
    ///
    /// The created CScalar is the same as s, but its spin m_x is scaled by c.
    inline CScalar(const CScalar &s,const SComplex &c): m_x(s.m_x*c)
    {
      m_h=s.m_h; m_s=s.m_s; m_c[0]=s.m_c[0]; m_c[1]=s.m_c[1];
    }
    
    /// @}

    /// @name Mathematical operations
    /// @{
    /// These operations apply only to the CScalar::m_x member (spin components).

    /// @note The argument passed to this function should be a CScalar, otherwise the
    /// behavior is undefined
    void Add(const CObject *c);
    void Divide(const double &d);
    void Multiply(const Complex &c);
    void Invert();
    /// @}

    /// @name Access operators
    /// @{
    /// Return CScalar::m_x, no matter what the argument value is.
    inline SComplex &operator[](const int i) { return m_x;    }
    inline const SComplex &operator[](const int i) const { return m_x;    }
    /// @}

    /// @name Mathematical operators
    /// @{
    /// These operators apply only to the CScalar::m_x member (spin components).
    inline CScalar operator+(const CScalar &s) const  
    { 
      return CScalar(m_x+s.m_x,m_c[0],m_c[1],m_h,m_s); 
    }
    inline CScalar operator-(const CScalar &s) const
    { 
      return CScalar(m_x-s.m_x,m_c[0],m_c[1],m_h,m_s); 
    }
    inline CScalar operator-() const
    { 
      return CScalar(-m_x,m_c[0],m_c[1],m_h,m_s); 
    }

    inline CScalar& operator+=(const CScalar &s) 
    {
      m_x+=s.m_x;
      return *this;
    }
    inline CScalar& operator-=(const CScalar &s) 
    {
      m_x-=s.m_x;
      return *this;
    }
    inline CScalar& operator*=(const SComplex &c) 
    {
      m_x*=c;
      return *this;
    }
    /// @}

    /// @brief Return a CScalar with a complex conjugate m_x.
    inline CScalar Conj() const 
    {
      return CScalar(std::conj(m_x),m_c[0],m_c[1],m_h,m_s);
    }
    /// @brief Return the square of the magnitude of m_x
    inline SComplex Abs2() const 
    {
      return m_x*m_x;
    }
    /// @brief Return the magnitude of m_x
    inline SComplex Abs() const 
    { 
      return sqrt(Abs2()); 
    }

    /// @brief Return true if m_x is NAN
    inline bool Nan() const { return ATOOLS::IsNan<SComplex>(m_x); }

    
    /// @name Access or manipulate the CScaler::s_accu data member
    /// @{

    /// @todo This function is currently buggy. It does not reset s_accu.
    static void ResetAccu();
    inline static void   SetAccu(const Scalar &accu) { s_accu=accu;   }
    inline static Scalar Accu()                      { return s_accu; }
    /// @}

  };// end of class CScalar

  template <class Scalar> inline std::complex<Scalar>
  operator*(const Scalar &d,const CScalar<Scalar> &s)
  { return CScalar<Scalar>(s,d); }
  template <class Scalar> inline CScalar<Scalar> 
  operator*(const CScalar<Scalar> &s,const Scalar &d)
  { return CScalar<Scalar>(s,d); }
  template <class Scalar> inline CScalar<Scalar>
  operator*(const std::complex<Scalar> &c,const CScalar<Scalar> &s)
  { return CScalar<Scalar>(s,c); }
  template <class Scalar> inline CScalar<Scalar> 
  operator*(const CScalar<Scalar> &s,const std::complex<Scalar> &c)
  { return CScalar<Scalar>(s,c); }
  template <class Scalar> inline CScalar<Scalar>
  operator/(const CScalar<Scalar> &s,const std::complex<Scalar> &c)
  { return CScalar<Scalar>(s,1.0/c); }

  template <class Scalar> inline std::complex<Scalar> 
  operator*(const CScalar<Scalar> &s1,const CScalar<Scalar> &s2) 
  { return s1[0]*s2[0]; }

  template <class Scalar>
  std::ostream &operator<<(std::ostream &str,const CScalar<Scalar> &s);

}// end of namespace ATOOLS

#define DCScalar METOOLS::CScalar<double>
#define QCScalar METOOLS::CScalar<long double>

#endif
