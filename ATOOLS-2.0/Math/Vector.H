/*
 * Vector.H: Declarations for vec.C
 * 3 dim. euklidean vectors and 4 dim. Minkowski vectors
 * Poincare transformations
 */

#ifndef vector_h
#define vector_h

#include <iostream>
#include "MathTools.H"
 
namespace AMATOOLS {
  class vec4d;
  class vec3d;

  class Tag {
    static const struct Tsum {} sum;
    static const struct Tdiff {} diff;
    static const struct Tsmul {} smul;
    static const struct Tcross {} cross;

    friend class vec4d;
    friend vec4d operator+ (const vec4d&, const vec4d&);
    friend vec4d operator- (const vec4d&, const vec4d&);
    friend vec4d operator* (const double, const vec4d&);

    friend class vec3d;
    friend vec3d operator+ (const vec3d&, const vec3d&);
    friend vec3d operator- (const vec3d&, const vec3d&);
    friend vec3d operator* (const double, const vec3d&);
    friend vec3d operator* (const vec3d&, const double);
    friend vec3d cross (const vec3d&, const vec3d&);
  };


  class vec3d {
    double x[3];
  public:
    vec3d(double x0, double x1, double x2){
      x[0]=x0;x[1]=x1;x[2]=x2;};
    vec3d(const vec4d& v);
    vec3d(){
      x[0]=x[1]=x[2]=0;};
    ~vec3d() {}

    // computational constructors:
    vec3d(const vec3d& v1 ,const vec3d& v2, const Tag::Tsum);
    vec3d(const vec3d& v1 ,const vec3d& v2, const Tag::Tdiff);
    vec3d(const vec3d& v1 ,const double scal, const Tag::Tsmul); 
    vec3d(const vec3d& a ,const vec3d& b, const Tag::Tcross); 

    inline const double abs() const;
    inline const double sqr() const;
    inline double& operator[] (int idx);
    inline const double operator[] (int idx) const;
  };

  class vec4d {
    double x[4];
  public:
    vec4d(const double x0, const double x1, const double x2, const double x3) { 
      x[0]=x0; x[1]=x1; x[2]=x2; x[3]=x3; };
    vec4d(const double E, const vec3d& xn){
      x[0]=E;x[1]=xn[1];x[2]=xn[2];x[3]=xn[3];};
    vec4d(const vec4d& v1) {
      x[0]=v1[0]; x[1]=v1[1]; x[2]=v1[2]; x[3]=v1[3];};
    vec4d(){ x[0]=x[1]=x[2]=x[3]=0.0;};

    // computational constructors:
    vec4d(const vec4d& v1 ,const vec4d& v2, const Tag::Tsum);
    vec4d(const vec4d& v1 ,const vec4d& v2, const Tag::Tdiff);
    vec4d(const vec4d& v1 ,const double scal, const Tag::Tsmul); 
  
    ~vec4d() {}
    inline double& operator[] (int idx);
    inline const double operator[] (int idx) const;
    const double abs2() const {return x[0]*x[0]-x[1]*x[1]-x[2]*x[2]-x[3]*x[3];};
    vec4d& operator+= (const vec4d& v);
    vec4d& operator-= (const vec4d& v);
    vec4d& operator*= (const double scal);

  };


  inline vec3d::vec3d(const vec4d& v)
  {
    x[0]=v[1];x[1]=v[2];x[2]=v[3];
  }

  inline double& vec3d::operator[] (int idx) 
  {
#ifdef CHECK
    if(idx<1 || idx>3) {
      cerr<<"vec3d: out of bound.\n";
      return x[0];
    }
#endif
    
    return x[--idx];
  }

  inline const double vec3d::operator[] (int idx) const 
  {
#ifdef CHECK
    if(idx<1 || idx>3) {
      cerr<<"vec3d: out of bound.\n";
      return x[0];
    }
#endif

    return x[--idx];
  }

  inline const double vec3d::sqr() const {
    return x[0]*x[0]+x[1]*x[1]+x[2]*x[2];
  }

  inline const double vec3d::abs() const{
    return sqrt(sqr());
  }

  // computational constructors:
  inline vec3d::vec3d(const vec3d& v1 ,const vec3d& v2, const Tag::Tsum) 
  {
    x[0]= v1[1]+v2[1];
    x[1]= v1[2]+v2[2];
    x[2]= v1[3]+v2[3];
  }

  inline vec3d::vec3d(const vec3d& v1 ,const vec3d& v2, const Tag::Tdiff) 
  {
    x[0]= v1[1]-v2[1];
    x[1]= v1[2]-v2[2];
    x[2]= v1[3]-v2[3];

  }

  inline vec3d::vec3d(const vec3d& v1 ,const double scal, const Tag::Tsmul) 
  {
    x[0]=scal*v1[1];
    x[1]=scal*v1[2];
    x[2]=scal*v1[3];
  }

  inline vec3d::vec3d(const vec3d& a ,const vec3d& b, const Tag::Tcross) 
  {
    x[0]=a[2]*b[3]-a[3]*b[2];
    x[1]=a[3]*b[1]-a[1]*b[3];
    x[2]=a[1]*b[2]-a[2]*b[1];
  }


  // new operator definitions:
  inline vec3d operator+ (const vec3d& v1, const vec3d& v2) {
    return vec3d(v1,v2,Tag::sum);
  }

  inline vec3d operator- (const vec3d& v1, const vec3d& v2) {
    return vec3d(v1,v2,Tag::diff);
  }

  inline vec3d operator* (const double scal,const vec3d& v) {
    return vec3d(v,scal,Tag::smul);
  }

  inline double operator* (const vec3d& v1, const vec3d& v2) {
    return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3];
  }

  inline vec3d operator* (const vec3d& v, const double scal) {
    return vec3d(v,scal,Tag::smul);
  }

  inline vec3d operator/ (const vec3d& v, const double scal) {
    return (1./scal)*v;
  }


  inline vec3d cross(const vec3d& a, const  vec3d& b)
  {
    return vec3d(a,b,Tag::cross);
  }



  // **********************************************************************
  //   class vec4d
  // **********************************************************************

  inline double& vec4d::operator[] (int idx) 
  {
#ifdef CHECK
    if(idx<0 || idx>3) {
      cerr<<"vec4d: out of bound.\n";
      return x[0];
    }
#endif

    return x[idx];
  }

  inline const double vec4d::operator[] (int idx) const 
  {
#ifdef CHECK
    if(idx<0 || idx>3) {
      cerr<<"vec4d: out of bound.\n";
      return x[0];
    }
#endif

    return x[idx];
  }

  inline vec4d& vec4d::operator+= (const vec4d& v) 
  {
    x[0] += v[0];
    x[1] += v[1];
    x[2] += v[2];
    x[3] += v[3];
    return *this;
  }

  inline vec4d& vec4d::operator-= (const vec4d& v) 
  {
    x[0] -= v[0];
    x[1] -= v[1];
    x[2] -= v[2];
    x[3] -= v[3];
    return *this;
  }

  inline vec4d& vec4d::operator*= (const double scal) 
  {
    x[0] *=scal;
    x[1] *=scal;
    x[2] *=scal;
    x[3] *=scal;
    return *this;
  }

  inline double operator* (const vec4d& v1, const vec4d& v2) 
  {
    return v1[0]*v2[0]-v1[1]*v2[1]-v1[2]*v2[2]-v1[3]*v2[3];
  }

  inline bool operator==(const vec4d& v1, const vec4d& v2) 
  {
    double maxp=Max(v1[0],Max(v1[1],Max(v1[2],v1[3]))); 
    double q=1.;
    if (!IsZero(maxp)) q=1./maxp;
    for(short int i=0;i<4;i++) {
      if (!IsZero(q*(v1[i]-v2[i]))) return false;
    }
    return true;
  }

  inline bool operator!=(const vec4d& v1, const vec4d& v2) 
  {
    return !(v1==v2);
  }

  // computational constructors
  inline vec4d::vec4d(const vec4d& v1 ,const vec4d& v2, const Tag::Tsum) 
  {
    x[0]=v1[0]+v2[0];
    x[1]=v1[1]+v2[1];
    x[2]=v1[2]+v2[2];
    x[3]=v1[3]+v2[3];
  }

  inline vec4d::vec4d(const vec4d& v1 ,const vec4d& v2, const Tag::Tdiff) 
  {
    x[0]=v1[0]-v2[0];
    x[1]=v1[1]-v2[1];
    x[2]=v1[2]-v2[2];
    x[3]=v1[3]-v2[3];
  }

  inline vec4d::vec4d(const vec4d& v1 ,const double scal, const Tag::Tsmul) 
  {
    x[0]=scal*v1[0];
    x[1]=scal*v1[1];
    x[2]=scal*v1[2];
    x[3]=scal*v1[3];
  }


  // new operator definitions:
  inline vec4d operator* (const double scal, const vec4d& v1) {
    return vec4d(v1,scal,Tag::smul);
  }

  inline vec4d operator+ (const vec4d& v1, const vec4d& v2) 
  {
    return vec4d(v1,v2,Tag::sum);
  }

  inline vec4d operator- (const vec4d& v1, const vec4d& v2) 
  {
    return vec4d(v1,v2,Tag::diff);
  }


  // output streams:
  std::ostream& operator<<(std::ostream& s, const vec4d& vec);
  std::ostream& operator<<(std::ostream& s, const vec3d& vec);


  // standard vectors:
  const static vec4d xvec;
  const static vec4d yvec;
  const static vec4d zvec;

  // end of namespace
};

#endif
