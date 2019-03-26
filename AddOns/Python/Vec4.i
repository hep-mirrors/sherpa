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

  public:

    Vec4();

    Vec4(const Scalar& x0, const Scalar& x1,
	 const Scalar& x2, const Scalar& x3);

    Scalar Mass() const;

    %extend{

      Scalar __getitem__(unsigned int i){
	return (*self)[i];
      };

      Scalar __mul__(const Vec4<Scalar>& v){
	return (*self)*v;
      };

      Vec4<Scalar> __mul__(const Scalar& s){
	return (*self)*s;
      };

      Vec4<Scalar> __rmul__(const Scalar& s){
	return (*self)*s;
      };

      Vec4<Scalar> __add__(const Vec4<Scalar>& v){
	return (*self)+v;
      };

      Vec4<Scalar> __sub__(const Vec4<Scalar>& v){
	return (*self)-v;
      }

      Vec4<Scalar> __neg__(){
	return -(*self);
      }

    };

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
