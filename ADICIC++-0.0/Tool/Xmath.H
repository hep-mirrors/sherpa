//bof
//Version: 1 ADICIC++-0.0/2004/07/12

//A collection of additional mathematical tools.



#ifndef _Xmath_H_
#define _Xmath_H_ _Xmath_H_


//#include





//namespace ... {



//==================================================================
//  Functions evaluating positive integer powers of double numbers
//==================================================================


inline double power(const double& x, unsigned i) {
  double r=1.0;
  for( ;i>0;--i) r*=x;
  return r;
}


template<unsigned i> inline double power(const double& x) {
  return power<i-1>(x)*x;
}
template<> double inline power<9>(const double& x) {
  return x*x*x*x*x*x*x*x*x;
}
template<> double inline power<8>(const double& x) {
  return x*x*x*x*x*x*x*x;
}
template<> double inline power<7>(const double& x) {
  return x*x*x*x*x*x*x;
}
template<> double inline power<6>(const double& x) {
  return x*x*x*x*x*x;
}
template<> double inline power<5>(const double& x) {
  return x*x*x*x*x;
}
template<> double inline power<4>(const double& x) {
  return x*x*x*x;
}
template<> double inline power<3>(const double& x) {
  return x*x*x;
}
template<> double inline power<2>(const double& x) {
  return x*x;
}
template<> double inline power<1>(const double& x) {
  return x;
}
template<> double inline power<0>(const double& x) {
  return 1.0;
}



//==========================================
//  Functions evaluating factorial numbers
//==========================================


inline double factorial(int j) {
  double r=1.0;
  for( ;j>1;--j) r*=j;
  return r;
}

//OR:

//inline double factorial(unsigned j) {
//  return (j==0) ? 1 : j*factorial(j-1);
//}



template<unsigned j> struct Factorial {
  static const double RET;
};
template<> struct Factorial<0> {
  static const double RET;
};
template<unsigned j>
const double Factorial<j>::RET=j*Factorial<j-1>::RET;

//OR:

//template<unsigned j> struct Factorial {
//  enum { RET=j*Factorial<j-1>::RET};
//};
//template<> struct Factorial<0> {
//  enum { RET=1};
//};



//}//eo namespace ...





#endif    //eo _Xmath_H_



//eof
