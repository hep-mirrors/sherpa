#ifndef Vector_H
#define Vector_H

/*
 * Vector.H: Declarations for vec.C
 * 3 dim. euklidean vectors and 4 dim. Minkowski vectors
 * Poincare transformations
 */


#include <iostream>
#include "MathTools.H"
 
namespace AMATOOLS {
  class Vec4D;
  class Vec3D;

  class Tag {
    static const struct Tsum {} sum;
    static const struct Tdiff {} diff;
    static const struct Tsmul {} smul;
    static const struct Tcross {} cross;

    friend class Vec4D;
    friend Vec4D operator+ (const Vec4D&, const Vec4D&);
    friend Vec4D operator- (const Vec4D&, const Vec4D&);
    friend Vec4D operator* (const double, const Vec4D&);

    friend class Vec3D;
    friend Vec3D operator+ (const Vec3D&, const Vec3D&);
    friend Vec3D operator- (const Vec3D&, const Vec3D&);
    friend Vec3D operator* (const double, const Vec3D&);
    friend Vec3D operator* (const Vec3D&, const double);
    friend Vec3D cross (const Vec3D&, const Vec3D&);
  };


  class Vec3D {
    double m_x[3];
  public:
    Vec3D(double x, double y, double z){
      m_x[0] = x;
      m_x[1] = y;
      m_x[2] = z;
    }
    Vec3D(const Vec4D& v);
    Vec3D(){
      m_x[0]=m_x[1]=m_x[2]=0.;
    }

    //!  Computational constructor
    inline Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tsum);
    //!  Computational constructor
    inline Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tdiff);
    //!  Computational constructor
    inline Vec3D(const Vec3D& v1 ,const double scal, const Tag::Tsmul); 
    //!  Computational constructor
    inline Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tcross); 


    inline const double Abs() const;
    inline const double Sqr() const;
    inline double& operator[] (int i);
    inline const double operator[] (int i) const;
  };

  class Vec4D {
    double m_x[4];
  public:
    Vec4D() {
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=0.;
    }
    Vec4D(const double x0, const double x1, const double x2, const double x3) { 
      m_x[0] = x0; 
      m_x[1] = x1; 
      m_x[2] = x2; 
      m_x[3] = x3; 
    }
    Vec4D(const double E, const Vec3D& xn) {
      m_x[0] = E;
      m_x[1] = xn[1];
      m_x[2] = xn[2];
      m_x[3] = xn[3];
    }

    //! Computational constructor
    Vec4D(const Vec4D& v1 ,const Vec4D& v2, const Tag::Tsum);
    //! Computational constructor
    Vec4D(const Vec4D& v1 ,const Vec4D& v2, const Tag::Tdiff);
    //! Computational constructor
    Vec4D(const Vec4D& v1 ,const double scal, const Tag::Tsmul); 
  
    inline double& operator[] (int i);
    inline const double operator[] (int i) const;
    const double Abs2() const {return m_x[0]*m_x[0]-m_x[1]*m_x[1]-m_x[2]*m_x[2]-m_x[3]*m_x[3];};

    Vec4D& operator+= (const Vec4D& v);
    Vec4D& operator-= (const Vec4D& v);
    Vec4D& operator*= (const double scal);

  };


  inline Vec3D::Vec3D(const Vec4D& v)
  {
    m_x[0]=v[1];m_x[1]=v[2];m_x[2]=v[3];
  }

  inline double& Vec3D::operator[] (int i) 
  {
#ifdef CHECK
    if(i<1 || i>3) {
      cerr<<"Vec3D: out of bound.\n";
      return m_x[0];
    }
#endif
    
    return m_x[--i];
  }

  inline const double Vec3D::operator[] (int i) const 
  {
#ifdef CHECK
    if(i<1 || i>3) {
      cerr<<"Vec3D: out of bound.\n";
      return m_x[0];
    }
#endif

    return m_x[--i];
  }

  inline const double Vec3D::Sqr() const {
    return m_x[0]*m_x[0]+m_x[1]*m_x[1]+m_x[2]*m_x[2];
  }

  inline const double Vec3D::Abs() const{
    return sqrt(Sqr());
  }

  inline Vec3D::Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tsum) 
  {
    m_x[0]= v1[1]+v2[1];
    m_x[1]= v1[2]+v2[2];
    m_x[2]= v1[3]+v2[3];
  }

  inline Vec3D::Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tdiff) 
  {
    m_x[0]= v1[1]-v2[1];
    m_x[1]= v1[2]-v2[2];
    m_x[2]= v1[3]-v2[3];

  }

  inline Vec3D::Vec3D(const Vec3D& v1 ,const double scal, const Tag::Tsmul) 
  {
    m_x[0]=scal*v1[1];
    m_x[1]=scal*v1[2];
    m_x[2]=scal*v1[3];
  }

  inline Vec3D::Vec3D(const Vec3D& v1 ,const Vec3D& v2, const Tag::Tcross) 
  {
    m_x[0]=v1[2]*v2[3]-v1[3]*v2[2];
    m_x[1]=v1[3]*v2[1]-v1[1]*v2[3];
    m_x[2]=v1[1]*v2[2]-v1[2]*v2[1];
  }

  // new operator definitions:
  inline Vec3D operator+ (const Vec3D& v1, const Vec3D& v2) {
    return Vec3D(v1,v2,Tag::sum);
  }

  inline Vec3D operator- (const Vec3D& v1, const Vec3D& v2) {
    return Vec3D(v1,v2,Tag::diff);
  }

  inline Vec3D operator* (const double scal,const Vec3D& v) {
    return Vec3D(v,scal,Tag::smul);
  }

  inline double operator* (const Vec3D& v1, const Vec3D& v2) {
    return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
  }

  inline Vec3D operator* (const Vec3D& v, const double scal) {
    return Vec3D(v,scal,Tag::smul);
  }

  inline Vec3D operator/ (const Vec3D& v, const double scal) {
    return (1./scal)*v;
  }


  inline Vec3D cross(const Vec3D& a, const  Vec3D& b)
  {
    return Vec3D(a,b,Tag::cross);
  }



  // **********************************************************************
  //   class Vec4D
  // **********************************************************************

  inline double& Vec4D::operator[] (int i) 
  {
#ifdef CHECK
    if(i<0 || i>3) {
      cerr<<"Vec4D: out of bound.\n";
      return m_x[0];
    }
#endif

    return m_x[i];
  }

  inline const double Vec4D::operator[] (int i) const 
  {
#ifdef CHECK
    if(i<0 || i>3) {
      cerr<<"Vec4D: out of bound.\n";
      return m_x[0];
    }
#endif

    return m_x[i];
  }

  inline Vec4D& Vec4D::operator+= (const Vec4D& v) 
  {
    m_x[0] += v[0];
    m_x[1] += v[1];
    m_x[2] += v[2];
    m_x[3] += v[3];
    return *this;
  }

  inline Vec4D& Vec4D::operator-= (const Vec4D& v) 
  {
    m_x[0] -= v[0];
    m_x[1] -= v[1];
    m_x[2] -= v[2];
    m_x[3] -= v[3];
    return *this;
  }

  inline Vec4D& Vec4D::operator*= (const double scal) 
  {
    m_x[0] *= scal;
    m_x[1] *= scal;
    m_x[2] *= scal;
    m_x[3] *= scal;
    return *this;
  }

  inline double operator* (const Vec4D& v1, const Vec4D& v2) 
  {
    return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];
  }

  inline bool operator==(const Vec4D& v1, const Vec4D& v2) 
  {
    double maxp=Max(v1[0],Max(v1[1],Max(v1[2],v1[3]))); 
    double q=1.;
    if (!IsZero(maxp)) q=1./maxp;
    for(short int i=0;i<4;i++) {
      if (!IsZero(q*(v1[i]-v2[i]))) return false;
    }
    return true;
  }

  inline bool operator!=(const Vec4D& v1, const Vec4D& v2) 
  {
    return !(v1==v2);
  }

  inline Vec4D::Vec4D(const Vec4D& v1 ,const Vec4D& v2, const Tag::Tsum) 
  {
    m_x[0] = v1[0]+v2[0];
    m_x[1] = v1[1]+v2[1];
    m_x[2] = v1[2]+v2[2];
    m_x[3] = v1[3]+v2[3];
  }

  inline Vec4D::Vec4D(const Vec4D& v1 ,const Vec4D& v2, const Tag::Tdiff) 
  {
    m_x[0] = v1[0]-v2[0];
    m_x[1] = v1[1]-v2[1];
    m_x[2] = v1[2]-v2[2];
    m_x[3] = v1[3]-v2[3];
  }

  inline Vec4D::Vec4D(const Vec4D& v1 ,const double scal, const Tag::Tsmul) 
  {
    m_x[0] = scal*v1[0];
    m_x[1] = scal*v1[1];
    m_x[2] = scal*v1[2];
    m_x[3] = scal*v1[3];
  }

  // new operator definitions:
  inline Vec4D operator* (const double scal, const Vec4D& v1) {
    return Vec4D(v1,scal,Tag::smul);
  }

  inline Vec4D operator+ (const Vec4D& v1, const Vec4D& v2) 
  {
    return Vec4D(v1,v2,Tag::sum);
  }

  inline Vec4D operator- (const Vec4D& v1, const Vec4D& v2) 
  {
    return Vec4D(v1,v2,Tag::diff);
  }


  // output streams:
  std::ostream& operator<<(std::ostream& s, const Vec4D& vec);
  std::ostream& operator<<(std::ostream& s, const Vec3D& vec);


  // standard vectors:
  const static Vec4D XVEC;
  const static Vec4D YVEC;
  const static Vec4D ZVEC;

  // end of namespace
}


// --------------------------------------------------
//           Doxygen part starts here
// --------------------------------------------------


/*!
 \file
 \brief   contains class AMATOOLS::Vec3D and class AMATOOLS::Vec4D 
*/

// --- class Vec3D ---

/*!
 \class AMATOOLS::Vec3D
 \brief implementation of a 3 dimensional vector and its algebraic operations

 This class can be used as a (real) 3 dimensional vector. All necessary
 operations, e.g. addition, scalar product, cross product, etc. are available.
*/

/*!
 \fn AMATOOLS::Vec3D::Vec3D()
 \brief Standard Constructor
*/

/*!
 \fn AMATOOLS::Vec3D::Vec3D(double x, double y, double z){
 \brief Special Constructor taking 3 single components
*/

/*!
 \fn AMATOOLS::Vec3D::Vec3D(const Vec4D& v)
 \brief Special Constructor extracting the space part of a minkowski vector
*/

/*!
 \fn inline const double AMATOOLS::Vec3D::Abs() const
 \brief returns \f$ \sqrt{x^2 + y^2 + z^2} \f$.
*/

/*!
 \fn inline const double AMATOOLS::Vec3D::Sqr() const
 \brief returns \f$ x^2 + y^2 + z^2 \f$.
*/

/*!
 \fn   inline double& AMATOOLS::Vec3D::operator[] (int i)
 \brief returns x,y,z for i=1,2,3 resp. (May be manipulated.)
*/

/*!
 \fn   inline const double AMATOOLS::Vec3D::operator[] (int i) const
 \brief returns x,y,z for i=1,2,3 resp.
*/

// --- class Vec4D ---

/*!
 \class AMATOOLS::Vec4D
 \brief implementation of a 4 dimensional (Minkowski) vector and its algebraic operations

 This class can be used as Minkowski vector. All necessary  operations, e.g. addition,
 subtraction, scalar product, etc. are available.
*/

/*!
 \fn AMATOOLS::Vec4D::Vec4D()
 \brief Standard Constructor
*/

/*!
 \fn AMATOOLS::Vec4D::Vec4D(double x0,double x1, double x2, double x3){
 \brief Special Constructor taking 4 single components
*/

/*!
 \fn AMATOOLS::Vec4D::Vec4D(const double E, const Vec3D& v)
 \brief Special Constructor taking the space part and a Energie
*/

/*!
 \fn inline const double AMATOOLS::Vec4D::Abs2() const
 \brief returns \f$ x_0^2 - (x_1^2 + x_2^2 + x_3^2) \f$.
*/

/*!
 \fn   inline double& AMATOOLS::Vec4D::operator[] (int i)
 \brief returns \f$x_i\f$. (May be manipulated.)
*/

/*!
 \fn   inline const double AMATOOLS::Vec4D::operator[] (int i) const
 \brief returns \f$x_i\f$.
*/


#endif
