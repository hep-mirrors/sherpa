#ifndef METOOLS_Explicit_C_Spinor_H
#define METOOLS_Explicit_C_Spinor_H

#include "METOOLS/Currents/C_Vector.H"
#include "ATOOLS/Phys/Spinor.H"
#include "ATOOLS/Org/STL_Tools.H"

#include <vector>

namespace METOOLS {

  template <class Scalar>
  class CSpinor: public CObject {
  public:

    typedef std::complex<Scalar> SComplex;

  protected:
    
    static double s_accu;
    static int    s_hsc;

    int m_r, m_b, m_on;

    SComplex m_u[4];

    static ATOOLS::AutoDelete_Vector<CSpinor> s_objects;

    template <class _Scalar> friend std::ostream &
    operator<<(std::ostream &ostr,const CSpinor<_Scalar> &s); 

  public:

    static CSpinor *New();
    static CSpinor *New(const CSpinor &s);
    static CSpinor *New(const int r,const int b,
			const int cr=0,const int ca=0,
			const size_t &h=0,const size_t &s=0,
			const int on=3);

    CObject* Copy() const;

    void Delete();

    bool IsZero() const;

    inline CSpinor(const int r=1,const int b=1,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0,
		   const int on=3): 
      m_r(r), m_b(b), m_on(on)
    { m_u[0]=m_u[1]=m_u[2]=m_u[3]=Scalar(0.0);
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca; }
    inline CSpinor(const int &r,const int &b,
		   const SComplex &u1,const SComplex &u2,
		   const SComplex &u3,const SComplex &u4,
		   const int cr=0,const int ca=0,
		   const size_t &h=0,const size_t &s=0,
		   const int on=3): 
      m_r(r), m_b(b), m_on(on)
    { m_u[0]=u1; m_u[1]=u2; m_u[2]=u3; m_u[3]=u4;
      m_h=h; m_s=s; m_c[0]=cr; m_c[1]=ca; }
    inline CSpinor(const int &r,const int &b,const int &h,
		   const ATOOLS::Vec4<Scalar> &p,
		   const int cr=0,const int ca=0,
		   const size_t &hh=0,const size_t &s=0,
		   const Scalar &m2=-1.0): 
      m_r(r), m_b(b)
    { m_h=hh; m_s=s; m_c[0]=cr; m_c[1]=ca; Construct(h,p,m2); }
    
    // member functions
    void Add(const CObject *c);
    void Divide(const double &d);
    void Invert();

    void Construct(const int h,const ATOOLS::Vec4<Scalar> &p,Scalar m2=-1.0);
    bool SetOn();

    SComplex operator*(const CSpinor &s) const;

    CSpinor operator*(const Scalar &d) const;
    CSpinor operator*(const SComplex &c) const;
    CSpinor operator/(const Scalar &d) const;
    CSpinor operator/(const SComplex &c) const;

    CSpinor operator*=(const Scalar &d); 
    CSpinor operator*=(const SComplex &c); 
    CSpinor operator/=(const Scalar &d);
    CSpinor operator/=(const SComplex &c); 

    CSpinor operator+(const CSpinor &s) const;
    CSpinor operator-(const CSpinor &s) const;

    CSpinor operator+=(const CSpinor &s); 
    CSpinor operator-=(const CSpinor &s);

    bool operator==(const CSpinor &s) const;

    // inline functions
    inline SComplex &operator[](const size_t &i) 
    { return m_u[i]; }
    inline const SComplex &operator[](const size_t &i) const 
    { return m_u[i]; }

    inline int R() const { return m_r; }
    inline int B() const { return m_b; }

    inline int On() const { return m_on; }

    inline CSpinor operator-() const
    { return CSpinor(m_r,m_b,-m_u[0],-m_u[1],-m_u[2],-m_u[3],
		     m_c[0],m_c[1],m_h,m_s,m_on); }

    inline CSpinor Bar() const
    { return CSpinor(m_r,-m_b,std::conj(m_u[2]),std::conj(m_u[3]),
		     std::conj(m_u[0]),std::conj(m_u[1]),m_c[0],m_c[1],
		     m_h,m_s,(m_on&1)<<1|(m_on&2)>>1); }
    
    inline static void SetAccuracy(const double &accu) 
    { s_accu=accu; }
    inline static void ResetAccuracy() 
    { s_accu=1.0e-12; }

    inline static double Accuracy() { return s_accu; }

    inline static void SetHSC(const int hsc) { s_hsc=hsc; }

  };// end of class CSpinor
  /*!
    \class CSpinor
    \brief Class representing a Dirac spinor.

    We employ \f$\gamma\f$-matrices in the Weyl representation, i.e.
    \f[
      \gamma_\mu=\left(\begin{array}{cc}0&\sigma_\mu\\
        \bar{\sigma}_\mu&0\end{array}\right)
    \f]
    where \f$\sigma_\mu=\left(1,-\vec{\sigma}\right)\f$ and 
    \f$\bar{\sigma}_\mu=\left(1,\vec{\sigma}\right)\f$.
    The pauli matrices \f$\sigma^i\f$ are given by
    \f[
      \sigma^1=\left(\begin{array}{cc}0&1\\1&0\end{array}\right)\;,\quad
      \sigma^2=\left(\begin{array}{cc}0&-i\\i&0\end{array}\right)\;,\quad
      \sigma^3=\left(\begin{array}{cc}1&0\\0&-1\end{array}\right)\;.
    \f]
    In this basis
    \f[
      \gamma^5=\left(\begin{array}{cc}-1&0\\0&1\end{array}\right)
    \f]
    The dirac equation thus yields the Eigenspinor problem
    \f[
      0=\left(p^\mu\gamma_\mu-m\right)u\\
       =\left(\begin{array}{cccc}
         -m&0&p^0-p^3&-p^1+ip^2\\0&-m&-p^1-ip^2&p^0+p^3\\
         p^0+p^3&p^1-ip^2&-m&0\\p^1+ip^2&p^0-p^3&0&-m\end{array}\right)u
    \f]
    Employing \f$p^\pm=p^0\pm p^3\f$ and \f$p_\perp=p^1+ip^2\f$ 
    this can be rewritten to give
    \f[
      0=\left(\begin{array}{cccc}
        -m&0&p^-&-p_\perp^*\\0&-m&-p_\perp&p^+\\
        p^+&p_\perp^*&-m&0\\p_\perp&p^-&0&-m\end{array}\right)u
    \f]
    The Eigenvalues are \f$\lambda=m\pm\sqrt{p^2}\f$. Eigenspinors
    can be found by firstly constructing solutions for \f$m=0\f$.
    One possible set of such solutions is
    \f[
      u_+(p,0)=v_-(p,0)=\left(\begin{array}{c}
        0\\\chi_+(p)\end{array}\right)\;,\quad
      u_-(p,0)=v_+(p,0)=\left(\begin{array}{c}
        \chi_-(p)\\0\end{array}\right)\;.
    \f]
    where
    \f[
      \chi_+(p)=\frac{1}{\sqrt{p^+}}\left(\begin{array}{c}
        p^+\\p_\perp\end{array}\right)=\left(\begin{array}{c}
        \sqrt{p^+}\\\sqrt{p^-}e^{i\phi_p}\end{array}\right)\;,\quad
      \chi_-(p)=\frac{e^{i\pi}}{\sqrt{p^+}}\left(\begin{array}{c}
        -p_\perp^*\\p^+\end{array}\right)=\left(\begin{array}{c}
        \sqrt{p^-}e^{-i\phi_p}\\-\sqrt{p^+}\end{array}\right)\;.
    \f]
    The spinors \f$\chi_\pm(p)\f$ are normalized to \f$2\,p_0\f$.
    Defining \f$\bar p=|\vec p|\f$ and 
    \f$\hat p=(\,{\rm sgn}(p_0)\,\bar p,\vec p\;)\f$, possible solutions 
    for \f$m\neq 0\f$ are [Nucl. Phys. B274 (1986) 1-32]
    \f[
      u_+(p,m)=\frac{1}{\sqrt{2\,\bar p}}\left(\begin{array}{r}
          \sqrt{p_0-\bar p}\;\chi_+(\hat p)\\
	  \sqrt{p_0+\bar p}\;\chi_+(\hat p)\end{array}\right)\;,\quad
      v_-(p,m)=\frac{1}{\sqrt{2\,\bar p}}\left(\begin{array}{r}
          -\sqrt{p_0-\bar p}\;\chi_+(\hat p)\\
	  \sqrt{p_0+\bar p}\;\chi_+(\hat p)\end{array}\right)\;,
    \f]
    \f[
      u_-(p,m)=\frac{1}{\sqrt{2\,\bar p}}\left(\begin{array}{r}
          \sqrt{p_0+\bar p}\;\chi_-(\hat p)\\
          \sqrt{p_0-\bar p}\;\chi_-(\hat p)\end{array}\right)\;,\quad
      v_+(p,m)=\frac{1}{\sqrt{2\,\bar p}}\left(\begin{array}{r}
          \sqrt{p_0+\bar p}\;\chi_-(\hat p)\\
	  -\sqrt{p_0-\bar p}\;\chi_-(\hat p)\end{array}\right)\;.
    \f]
    These Eigenspinors are orthogonal and satisfy the relations
    \f$\bar u_\lambda(p,m)u_\lambda(p,m)=2m\f$ 
    and \f$\bar v_\lambda(p,m)v_\lambda(p,m)=-2m\f$. 
    In the above spinor basis the color-ordered quark-quark-gluon vertex reads
    \f[
      \frac{i}{\sqrt{2}}\bar{u}\gamma^\mu v\\=
      \frac{i}{\sqrt{2}}\left(\bar{u}_0,\bar{u}_1,\bar{u}_2,\bar{u}_3\right)
      \left(\begin{array}{c}\left(\begin{array}{cc}0&1\\1&0\end{array}\right)\\
        \left(\begin{array}{cc}0&\sigma^1\\-\sigma^1&0\end{array}\right)\\
        \left(\begin{array}{cc}0&\sigma^2\\-\sigma^2&0\end{array}\right)\\
        \left(\begin{array}{cc}0&\sigma^3\\-\sigma^3&0\end{array}\right)
      \end{array}\right)
      \left(\begin{array}{c}v_0\\v_1\\v_2\\v_3\end{array}\right)\\=
      \frac{i}{\sqrt{2}}\left(\begin{array}{c}
        \bar{u}_0v_2+\bar{u}_1v_3+\bar{u}_2v_0+\bar{u}_3v_1\\
        \bar{u}_0v_3+\bar{u}_1v_2-\bar{u}_2v_1-\bar{u}_3v_0\\
        -i\left(\bar{u}_0v_3-\bar{u}_1v_2-\bar{u}_2v_1+\bar{u}_3v_0\right)\\
        \bar{u}_0v_2-\bar{u}_1v_3-\bar{u}_2v_0+\bar{u}_3v_1\end{array}\right)
    \f]
    The color-ordered quark-gluon-quark vertices read
    \f[
      \frac{i}{\sqrt{2}}j^\mu\gamma_\mu v\\=
      \frac{i}{\sqrt{2}}\left(\begin{array}{cccc}
        0&0&j^-&-j_\perp^*\\0&0&-j_\perp&j^+\\
        j^+&j_\perp^*&0&0\\j_\perp&j^-&0&0\end{array}\right)
      \left(\begin{array}{c}v_0\\v_1\\v_2\\v_3\end{array}\right)\\=
      \frac{i}{\sqrt{2}}\left(\begin{array}{c}
        j^-v_2-j_\perp^*v_3\\-j_\perp v_2+j^+v_3\\
        j^+v_0+j_\perp^*v_1\\j_\perp v_0+j^-v_1\end{array}\right)
    \f]
    and
    \f[
      \frac{i}{\sqrt{2}}\bar{u}j^\mu\gamma_\mu\\=
      \frac{i}{\sqrt{2}}\left(\bar{u}_0,\bar{u}_1,\bar{u}_2,\bar{u}_3\right)
      \left(\begin{array}{cccc}0&0&j^-&-j_\perp^*\\0&0&-j_\perp&j^+\\
        j^+&j_\perp^*&0&0\\j_\perp&j^-&0&0\end{array}\right)\\=
      \frac{i}{\sqrt{2}}\left(
        \bar{u}_2j^++\bar{u}_3j_\perp,\bar{u}_2j_\perp^*+\bar{u}_3j^-,
        \bar{u}_0j^--\bar{u}_1j_\perp,-\bar{u}_0j_\perp^*+\bar{u}_1j^+\right)
    \f]
  */
    
  template <class Scalar>
  std::ostream &operator<<(std::ostream &ostr,const CSpinor<Scalar> &s); 

}// end of namespace BCF

#define DDSpinor METOOLS::CSpinor<double>
#define QDSpinor METOOLS::CSpinor<long double>

#endif
