//bof
//Version: 1 ADICIC++-0.0/2004/05/04

//A collection of additional mathematical tools.



#ifndef _Xmath_H_
#define _Xmath_H_ _Xmath_H_


//#include





//namespace ... {



//===================================================================
//  A function evaluating positive integer powers of double numbers
//===================================================================


template<unsigned short i> inline double power(double x) {
  return x*power<i-1>(x);
}
template<> inline double power<0>(double x) {
  return 1.0;
}
template<> inline double power<1>(double x) {
  return x;
}



//}//eo namespace ...





#endif    //eo _Xmath_H_



//eof
