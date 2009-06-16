#ifndef ATOOLS_Math_Tensor_Build_H
#define ATOOLS_Math_Tensor_Build_H

/********************************************************************************************
******                                                                                 ******
******       tensor builds                                                             ******
******                                                                                 ******
********************************************************************************************/

template<typename Scal1, typename Scal2> 
Lorentz_Ten2<PROMOTE(Scal1,Scal2)> 
ATOOLS::BuildTensor(const Vec4<Scal1>& p, const Vec4<Scal2>& q) {
  // T^{mu,nu} = p^{mu} q^{nu}
  PROMOTE(Scal1,Scal2) x[4][4];
  for (unsigned short int i=0; i<4; i++)
    for (unsigned short int j=0; j<4; j++)
      x[i][j] = p[i]*q[j];
  return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>(x);
//   return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>(p[0]*q[0],p[1]*q[0],p[2]*q[0],p[3]*q[0],
//                                             p[0]*q[1],p[1]*q[1],p[2]*q[1],p[3]*q[1],
//                                             p[0]*q[2],p[1]*q[2],p[2]*q[2],p[3]*q[2],
//                                             p[0]*q[3],p[1]*q[3],p[2]*q[3],p[3]*q[3]);
}

template<typename Scal1, typename Scal2, typename Scal3> 
Lorentz_Ten3<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))> 
ATOOLS::BuildTensor(const Vec4<Scal1>& p1, const Vec4<Scal2>& p2, const Vec4<Scal3>& p3) {
  // T^{mu,nu,rho} = p1^{mu} p2^{nu} p3^{rho}
  PROMOTE(Scal1,PROMOTE(Scal2,Scal3)) x[4][4][4];
  for (unsigned short int i=0; i<4; ++i)
    for (unsigned short int j=0; j<4; ++j)
      for (unsigned short int k=0; k<4; ++k)
        x[i][j][k] = p1[i]*p2[j]*p3[k];
  return Lorentz_Ten3<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))>(x);
//   return BuildTensor(BuildTensor(p1,p2),p3);
//   return Lorentz_Ten3<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))>(
//             p1[0]*p2[0]*p3[0],p1[1]*p2[0]*p3[0],p1[2]*p2[0]*p3[0],p1[3]*p2[0]*p3[0],
//             p1[0]*p2[1]*p3[0],p1[1]*p2[1]*p3[0],p1[2]*p2[1]*p3[0],p1[3]*p2[1]*p3[0],
//             p1[0]*p2[2]*p3[0],p1[1]*p2[2]*p3[0],p1[2]*p2[2]*p3[0],p1[3]*p2[2]*p3[0],
//             p1[0]*p2[3]*p3[0],p1[1]*p2[3]*p3[0],p1[2]*p2[3]*p3[0],p1[3]*p2[3]*p3[0],
//             p1[0]*p2[0]*p3[1],p1[1]*p2[0]*p3[1],p1[2]*p2[0]*p3[1],p1[3]*p2[0]*p3[1],
//             p1[0]*p2[1]*p3[1],p1[1]*p2[1]*p3[1],p1[2]*p2[1]*p3[1],p1[3]*p2[1]*p3[1],
//             p1[0]*p2[2]*p3[1],p1[1]*p2[2]*p3[1],p1[2]*p2[2]*p3[1],p1[3]*p2[2]*p3[1],
//             p1[0]*p2[3]*p3[1],p1[1]*p2[3]*p3[1],p1[2]*p2[3]*p3[1],p1[3]*p2[3]*p3[1],
//             p1[0]*p2[0]*p3[2],p1[1]*p2[0]*p3[2],p1[2]*p2[0]*p3[2],p1[3]*p2[0]*p3[2],
//             p1[0]*p2[1]*p3[2],p1[1]*p2[1]*p3[2],p1[2]*p2[1]*p3[2],p1[3]*p2[1]*p3[2],
//             p1[0]*p2[2]*p3[2],p1[1]*p2[2]*p3[2],p1[2]*p2[2]*p3[2],p1[3]*p2[2]*p3[2],
//             p1[0]*p2[3]*p3[2],p1[1]*p2[3]*p3[2],p1[2]*p2[3]*p3[2],p1[3]*p2[3]*p3[2],
//             p1[0]*p2[0]*p3[3],p1[1]*p2[0]*p3[3],p1[2]*p2[0]*p3[3],p1[3]*p2[0]*p3[3],
//             p1[0]*p2[1]*p3[3],p1[1]*p2[1]*p3[3],p1[2]*p2[1]*p3[3],p1[3]*p2[1]*p3[3],
//             p1[0]*p2[2]*p3[3],p1[1]*p2[2]*p3[3],p1[2]*p2[2]*p3[3],p1[3]*p2[2]*p3[3],
//             p1[0]*p2[3]*p3[3],p1[1]*p2[3]*p3[3],p1[2]*p2[3]*p3[3],p1[3]*p2[3]*p3[3]);
}

template<typename Scal1, typename Scal2> 
Lorentz_Ten3<PROMOTE(Scal1,Scal2)> 
ATOOLS::BuildTensor(const Lorentz_Ten2<Scal1>& t, const Vec4<Scal2>& p) {
  // T^{mu,nu,rho} = t^{mu,nu} p^{rho}
  PROMOTE(Scal1,Scal2) x[4][4][4];
  for (unsigned short int i=0; i<4; ++i)
    for (unsigned short int j=0; j<4; ++j)
      for (unsigned short int k=0; k<4; ++k)
        x[i][j][k] = t.at(i,j)*p[k];
  return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>(x);
//   return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>(
//             t.at(0,0)*p[0],t.at(1,0)*p[0],t.at(2,0)*p[0],t.at(3,0)*p[0],
//             t.at(0,1)*p[0],t.at(1,1)*p[0],t.at(2,1)*p[0],t.at(3,1)*p[0],
//             t.at(0,2)*p[0],t.at(1,2)*p[0],t.at(2,2)*p[0],t.at(3,2)*p[0],
//             t.at(0,3)*p[0],t.at(1,3)*p[0],t.at(2,3)*p[0],t.at(3,3)*p[0],
//             t.at(0,0)*p[1],t.at(1,0)*p[1],t.at(2,0)*p[1],t.at(3,0)*p[1],
//             t.at(0,1)*p[1],t.at(1,1)*p[1],t.at(2,1)*p[1],t.at(3,1)*p[1],
//             t.at(0,2)*p[1],t.at(1,2)*p[1],t.at(2,2)*p[1],t.at(3,2)*p[1],
//             t.at(0,3)*p[1],t.at(1,3)*p[1],t.at(2,3)*p[1],t.at(3,3)*p[1],
//             t.at(0,0)*p[2],t.at(1,0)*p[2],t.at(2,0)*p[2],t.at(3,0)*p[2],
//             t.at(0,1)*p[2],t.at(1,1)*p[2],t.at(2,1)*p[2],t.at(3,1)*p[2],
//             t.at(0,2)*p[2],t.at(1,2)*p[2],t.at(2,2)*p[2],t.at(3,2)*p[2],
//             t.at(0,3)*p[2],t.at(1,3)*p[2],t.at(2,3)*p[2],t.at(3,3)*p[2],
//             t.at(0,0)*p[3],t.at(1,0)*p[3],t.at(2,0)*p[3],t.at(3,0)*p[3],
//             t.at(0,1)*p[3],t.at(1,1)*p[3],t.at(2,1)*p[3],t.at(3,1)*p[3],
//             t.at(0,2)*p[3],t.at(1,2)*p[3],t.at(2,2)*p[3],t.at(3,2)*p[3],
//             t.at(0,3)*p[3],t.at(1,3)*p[3],t.at(2,3)*p[3],t.at(3,3)*p[3]);
}

template<typename Scal1, typename Scal2, typename Scal3, typename Scal4> 
Lorentz_Ten4<PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4)))> 
ATOOLS::BuildTensor(const Vec4<Scal1>& p1, const Vec4<Scal2>& p2,
                    const Vec4<Scal3>& p3, const Vec4<Scal4>& p4) {
  // T^{mu,nu,rho,sigma} = p1^{mu} p2^{nu} p3^{rho} p4^{sigma}
  PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4))) x[4][4][4][4];
  for (unsigned short int i=0; i<4; ++i)
    for (unsigned short int j=0; j<4; ++j)
      for (unsigned short int k=0; k<4; ++k)
        for (unsigned short int l=0; l<4; ++l)
          x[i][j][k][l] = p1[i]*p2[j]*p3[k]*p4[l];
  return Lorentz_Ten4<PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4)))>(x);
//   return BuildTensor(BuildTensor(p1,p2),BuildTensor(p3,p4));
//   return BuildTensor(BuildTensor(p1,p2,p3),p4);
}

template<typename Scal1, typename Scal2>
Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
ATOOLS::BuildTensor(const Lorentz_Ten2<Scal1>& t, const Lorentz_Ten2<Scal2>& s) {
  // T^{mu,nu,rho,sigma} = t^{mu,nu} s^{rho,sigma}
  PROMOTE(Scal1,Scal2) x[4][4][4][4];
  for (unsigned short int i=0; i<4; ++i)
    for (unsigned short int j=0; j<4; ++j)
      for (unsigned short int k=0; k<4; ++k)
        for (unsigned short int l=0; l<4; ++l)
          x[i][j][k][l] = t.at(i,j)*s.at(k,l);
  return Lorentz_Ten4<PROMOTE(Scal1,Scal2)>(x);
}

template<typename Scal1, typename Scal2>
Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
ATOOLS::BuildTensor(const Lorentz_Ten2<Scal1>& t, int i, int j,
                    const Lorentz_Ten2<Scal2>& s, int k, int l) {
  msg_Error()<<METHOD<<" not implemented."<<std::endl;
  abort();
  return Lorentz_Ten4<PROMOTE(Scal1,Scal2)>();
}

template<typename Scal1, typename Scal2> 
Lorentz_Ten4<PROMOTE(Scal1,Scal2)> 
ATOOLS::BuildTensor(const Lorentz_Ten3<Scal1>& t, const Vec4<Scal2>& p) {
  // T^{mu,nu,rho,sigma} = t^{mu,nu,rho} p^{sigma}
  PROMOTE(Scal1,Scal2) x[4][4][4][4];
  for (unsigned short int i=0; i<4; ++i)
    for (unsigned short int j=0; j<4; ++j)
      for (unsigned short int k=0; k<4; ++k)
        for (unsigned short int l=0; l<4; ++l)
          x[i][j][k][l] = t.at(i,j,k)*p[l];
  return Lorentz_Ten4<PROMOTE(Scal1,Scal2)>(x);
}

#endif
