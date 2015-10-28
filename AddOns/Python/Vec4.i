//%module Vec4
%{
#include <ATOOLS/Math/MathTools.H>
#include <ATOOLS/Math/Vec4.H>
#include <ATOOLS/Math/Vector.H>
#include <ATOOLS/Org/MyStrStream.H>
#include <iostream>

using namespace ATOOLS;
 
%}

namespace ATOOLS {
  template<typename Scalar> class Vec4;
  template<typename Scalar> class Vec3;

  typedef Vec4<double> Vec4D;
  typedef std::vector<Vec4D> Vec4D_Vector;

  template<typename Scalar> class Vec4 {
    Scalar m_x[4];
    //template <typename Scalar2> friend class Vec4;
  public:
    Vec4() {
      m_x[0]=m_x[1]=m_x[2]=m_x[3]=Scalar(0.0);
    }

    inline Vec4(const Scalar& x0, const Scalar& x1,
                const Scalar& x2, const Scalar& x3) {
      m_x[0] = x0; m_x[1] = x1; m_x[2] = x2; m_x[3] = x3;
    }

    // This extension is called when Vec4[i] is evaluated
    // NOTE: a proper typemap might defined if "Scalar" is not
    // of one of the standard C++ types (double/float/int)
    %extend{
        Scalar __getitem__(unsigned int i){
	  return (*self)[i];
        };
    };

    %rename(__sub__) operator-;

    inline Vec4<Scalar> operator-() const {
      return Vec4<Scalar>(-m_x[0],-m_x[1],-m_x[2],-m_x[3]);
    }

    double Mass() const { 
      return sqrt(ATOOLS::Abs<Scalar>(Abs2()));
    }

    %extend {
      std::string __str__() {
	MyStrStream conv;
	conv<<*self;
	return conv.str();
      };
    };

  };


  
}

%include <std_vector.i>

// Instantiate a "double" version of the Vec4-template that will be available as a Class Vec4D in python
%template(Vec4D) ATOOLS::Vec4<double>;
// Instantiate a "Vec4D" version of the std::vector-template that will be available as a Class Vec4D_Vector in python
%template(Vec4D_Vector) std::vector<ATOOLS::Vec4D>;
