#ifndef ATOOLS_Math_Tensor_Contractions_Epsilon_H
#define ATOOLS_Math_Tensor_Contractions_Epsilon_H

/********************************************************************************************
******                                                                                 ******
******       epsilon tensor contractions explicitely                                   ******
******                                                                                 ******
********************************************************************************************/

template<typename Scalar>
Lorentz_Ten3<Scalar>
ATOOLS::EpsilonTensorContraction(const Vec4<Scalar>& v) {
  // T^{mu,nu,rho} = eps^{mu,nu,rho,alpha} g_{alpha,alpha'} v^{alpha}
  return Lorentz_Ten3<Scalar>(
             0.  , 0.  , 0.  , 0.  ,
             0.  , 0.  ,-v[3], v[2],
             0.  , v[3], 0.  ,-v[1],
             0.  ,-v[2], v[1], 0.  ,
             0.  , 0.  , v[3],-v[2],
             0.  , 0.  , 0.  , 0.  ,
            -v[3], 0.  , 0.  , v[0],
             v[2], 0.  ,-v[0], 0.  ,
             0.  ,-v[3], 0.  , v[1],
             v[3], 0.  , 0.  ,-v[0],
             0.  , 0.  , 0.  , 0.  ,
            -v[1], v[0], 0.  , 0.  ,
             0.  , v[2],-v[1], 0.  ,
            -v[2], 0.  , v[0], 0.  ,
             v[1],-v[0], 0.  , 0.  ,
             0.  , 0.  , 0.  , 0.  );
}

template<typename Scal1, typename Scal2>
Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
ATOOLS::EpsilonTensorContraction(const Vec4<Scal1>& p, const Vec4<Scal2>& q) {
  // returns epsilon tensor contractions 
  //  T^{mu,nu} = eps^{mu,nu,alpha,beta}g_{alpha,alpha'}g_{beta,beta'}p^alpha' q^beta' 
  return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>( 
                               0.   , -p[2]*q[3]+p[3]*q[2] ,  p[1]*q[3]-p[3]*q[1] , -p[1]*q[2]+p[2]*q[1] ,
                p[2]*q[3]-p[3]*q[2] ,                 0.   , -p[0]*q[3]+p[3]*q[0] ,  p[0]*q[2]-p[2]*q[0] ,
               -p[1]*q[3]+p[3]*q[1] ,  p[0]*q[3]-p[3]*q[0] ,                 0.   , -p[0]*q[1]+p[1]*q[0] ,
                p[1]*q[2]-p[2]*q[1] , -p[0]*q[2]+p[2]*q[0] ,  p[0]*q[1]-p[1]*q[0] ,                 0.   );
}

template<typename Scal1, typename Scal2, typename Scal3>
Vec4<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))>
ATOOLS::EpsilonTensorContraction(const Vec4<Scal1>& p1, const Vec4<Scal2>& p2,
                                 const Vec4<Scal3>& p3) {
  // v^{mu} = eps^{mu,alpha,beta,gamma}g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  //              p1^{alpha'} p2^{beta'} p3^{gamma'} and permutations thereof
// maybe uncomment, but doubtful if any computational improvement
//   if (IsEqual(p1,p2) || IsEqual(p1,p3)  || IsEqual(p2,p3))
//     return Vec4<PROMOTE(Scal1,PROMOTE(Scal2,Scal3))>(0.,0.,0.,0.);
  return EpsilonTensorContraction(BuildTensor(p1,p2,p3),2,3,4);
}

template<typename Scal1, typename Scal2, typename Scal3, typename Scal4>
PROMOTE(Scal1,PROMOTE(Scal2,PROMOTE(Scal3,Scal4)))
ATOOLS::EpsilonTensorContraction(const Vec4<Scal1>& p1, const Vec4<Scal2>& p2,
                                 const Vec4<Scal3>& p3, const Vec4<Scal4>& p4) {
  // c = eps^{alpha,beta,gamma,delta} p1^alpha p2^beta p3^gamma p4^delta
  //                          gamma_{alpha,alpha'} gamma_{beta,beta'}
  //                          gamma_{gamma,gamma'} gamma_{delta,delta'}
  return   p1[0]*p2[1]*p3[2]*p4[3]
          -p1[0]*p2[1]*p3[3]*p4[2]
          -p1[0]*p2[2]*p3[1]*p4[3]
          +p1[0]*p2[2]*p3[3]*p4[1]
          +p1[0]*p2[3]*p3[1]*p4[2]
          -p1[0]*p2[3]*p3[2]*p4[1]

          -p1[1]*p2[0]*p3[2]*p4[3]
          +p1[1]*p2[0]*p3[3]*p4[2]
          +p1[1]*p2[2]*p3[0]*p4[3]
          -p1[1]*p2[2]*p3[3]*p4[0]
          -p1[1]*p2[3]*p3[0]*p4[2]
          +p1[1]*p2[3]*p3[2]*p4[0]

          +p1[2]*p2[0]*p3[1]*p4[3]
          -p1[2]*p2[0]*p3[3]*p4[1]
          -p1[2]*p2[1]*p3[0]*p4[3]
          +p1[2]*p2[1]*p3[3]*p4[0]
          +p1[2]*p2[3]*p3[0]*p4[1]
          -p1[2]*p2[3]*p3[1]*p4[0]

          -p1[3]*p2[0]*p3[1]*p4[2]
          +p1[3]*p2[0]*p3[2]*p4[1]
          +p1[3]*p2[1]*p3[0]*p4[2]
          -p1[3]*p2[1]*p3[2]*p4[0]
          -p1[3]*p2[2]*p3[0]*p4[1]
          +p1[3]*p2[2]*p3[1]*p4[0];
}

template<typename Scalar>
Lorentz_Ten2<Scalar>
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten2<Scalar>& t, int i, int j) {
// returns the contraction eps^{mu,nu,alpha,beta}g_{alpha,alpha'}g_{beta,beta'}t^{alpha'beta'}
  if ((i==3) && (j==4))
    return Lorentz_Ten2<Scalar>(
      0. , -t.at(2,3)+t.at(3,2) , t.at(1,3)-t.at(3,1) , -t.at(1,2)+t.at(2,1) ,
      t.at(2,3)-t.at(3,2) , 0. , -t.at(0,3)+t.at(3,0) ,  t.at(0,2)-t.at(2,0) ,
     -t.at(1,3)+t.at(3,1) ,  t.at(0,3)-t.at(3,0) , 0. , -t.at(0,1)+t.at(1,0) ,
      t.at(1,2)-t.at(2,1) , -t.at(0,2)+t.at(2,0) , t.at(0,1)-t.at(1,0) , 0. );
  else if ((j==3) && (i==4))
    return -Lorentz_Ten2<Scalar>(
      0. , -t.at(2,3)+t.at(3,2) , t.at(1,3)-t.at(3,1) , -t.at(1,2)+t.at(2,1) ,
      t.at(2,3)-t.at(3,2) , 0. , -t.at(0,3)+t.at(3,0) ,  t.at(0,2)-t.at(2,0) ,
     -t.at(1,3)+t.at(3,1) ,  t.at(0,3)-t.at(3,0) , 0. , -t.at(0,1)+t.at(1,0) ,
      t.at(1,2)-t.at(2,1) , -t.at(0,2)+t.at(2,0) , t.at(0,1)-t.at(1,0) , 0. );
  return Lorentz_Ten2<Scalar>();
}

template<typename Scal1, typename Scal2>
Vec4<PROMOTE(Scal1,Scal2)>
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten2<Scal1>& t, int i, int j,
                                 const Vec4<Scal2>& p, int k) {
  // v^{mu} = eps^{mu,alpha,beta,gamma}g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  //              t^{alpha',beta'} p^{gamma'} and permutations thereof
  return EpsilonTensorContraction(BuildTensor(t,p),i,j,k);
}

template<typename Scalar>
Vec4<Scalar>
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten3<Scalar>& ten, int i, int j, int k) {
  Lorentz_Ten3<Scalar> t;
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{alpha',beta',gamma'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  if      ((i==2) && (j==3) && (k==4)) 
     t = ten;
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{alpha',gamma',beta'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  else if ((i==2) && (j==4) && (k==3))
     t = ten.Transpose(2,3);
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{beta',alpha',gamma'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  else if ((i==3) && (j==2) && (k==4))
     t = ten.Transpose(1,2);
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{beta',gamma',alpha'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  else if ((i==3) && (j==4) && (k==2))
    {t = ten.Transpose(2,3); t = t.Transpose(1,2);}
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{gamma',alpha',beta'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  else if ((i==4) && (j==2) && (k==3))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3);}
  // v^{mu} = eps^{mu,alpha,beta,gamma} t^{gamma',beta',alpha'} g_{alpha,alpha'}g_{beta,beta'}g_{gamma,gamma'}
  else if ((i==4) && (j==3) && (k==2))
     t = ten.Transpose(1,3);
  else
    return Vec4<Scalar>(0.,0.,0.,0.);
  return Vec4<Scalar>(
                            t.at(3,2,1)-t.at(2,3,1)+t.at(1,3,2)-t.at(3,1,2)+t.at(2,1,3)-t.at(1,2,3),
    t.at(3,2,0)-t.at(2,3,0)+                        t.at(0,3,2)-t.at(3,0,2)+t.at(2,0,3)-t.at(0,2,3),
    t.at(1,3,0)-t.at(3,1,0)+t.at(3,0,1)-t.at(0,3,1)+                        t.at(0,1,3)-t.at(1,0,3),
    t.at(2,1,0)-t.at(1,2,0)+t.at(0,2,1)-t.at(2,1,0)+t.at(1,0,2)-t.at(0,1,2)                        );
}

template<typename Scal1, typename Scal2>
PROMOTE(Scal1,Scal2)
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten2<Scal1>& ten, int i, int j,
                                 const Lorentz_Ten2<Scal2>& sen, int k, int l) {
  // c = eps^{alpha,beta,gamma,delta} t^{alpha',beta'} s^{gamma',delta'}
  //                                  gamma_{alpha,alpha'} gamma_{beta,beta'}
  //                                  gamma_{gamma,gamma'} gamma_{delta,delta'}
  return EpsilonTensorContraction(BuildTensor(ten,sen),i,j,k,l);
}

template<typename Scal1, typename Scal2>
PROMOTE(Scal1,Scal2)
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten3<Scal1>& ten, int i, int j, int k,
                                 const Vec4<Scal2>& v, int l) {
  // c = eps^{alpha,beta,gamma,delta} t^{alpha',beta',gamma'} v^delta'
  //                                  gamma_{alpha,alpha'} gamma_{beta,beta'}
  //                                  gamma_{gamma,gamma'} gamma_{delta,delta'}
  return EpsilonTensorContraction(BuildTensor(ten,v),i,j,k,l);
}

template<typename Scalar>
Scalar
ATOOLS::EpsilonTensorContraction(const Lorentz_Ten4<Scalar>& ten, int i, int j, int k, int l) {
  Lorentz_Ten4<Scalar> t;
  // c = eps^{alpha,beta,gamma,delta} t^{alpha',beta',gamma',delta'}
  //                                  gamma_{alpha,alpha'} gamma_{beta,beta'}
  //                                  gamma_{gamma,gamma'} gamma_{delta,delta'}
  if      ((i == 1) && (j == 2) && (k == 3) && (l == 4))
    {t = ten;}
  else if ((i == 1) && (j == 2) && (k == 4) && (l == 3))
    {t = ten.Transpose(3,4);}
  else if ((i == 1) && (j == 3) && (k == 2) && (l == 4))
    {t = ten.Transpose(2,3);}
  else if ((i == 1) && (j == 3) && (k == 4) && (l == 2))
    {t = ten.Transpose(3,4); t = t.Transpose(2,3);}
  else if ((i == 1) && (j == 4) && (k == 2) && (l == 3))
    {t = ten.Transpose(2,3); t = t.Transpose(3,4);}
  else if ((i == 1) && (j == 4) && (k == 3) && (l == 2))
    {t = ten.Transpose(2,4);}
  else if ((i == 2) && (j == 1) && (k == 3) && (l == 4))
    {t = ten.Transpose(1,2);}
  else if ((i == 2) && (j == 1) && (k == 4) && (l == 3))
    {t = ten.Transpose(1,2); t = t.Transpose(3,4);}
  else if ((i == 2) && (j == 3) && (k == 1) && (l == 4))
    {t = ten.Transpose(2,3); t = t.Transpose(3,4);}
  else if ((i == 2) && (j == 3) && (k == 4) && (l == 1))
    {t = ten.Transpose(3,4); t = t.Transpose(2,3); t = t.Transpose(1,2);}
  else if ((i == 2) && (j == 4) && (k == 1) && (l == 3))
    {t = ten.Transpose(1,3); t = t.Transpose(2,3); t = t.Transpose(3,4);}
  else if ((i == 2) && (j == 4) && (k == 3) && (l == 1))
    {t = ten.Transpose(2,4); t = t.Transpose(1,2);}
  else if ((i == 3) && (j == 1) && (k == 2) && (l == 4))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3);}
  else if ((i == 3) && (j == 1) && (k == 4) && (l == 2))
    {t = ten.Transpose(1,2); t = t.Transpose(2,4); t = t.Transpose(3,4);}
  else if ((i == 3) && (j == 2) && (k == 1) && (l == 4))
    {t = ten.Transpose(1,3);}
  else if ((i == 3) && (j == 2) && (k == 4) && (l == 1))
    {t = ten.Transpose(1,4); t = t.Transpose(3,4);}
  else if ((i == 3) && (j == 4) && (k == 1) && (l == 2))
    {t = ten.Transpose(1,3); t = t.Transpose(2,4);}
  else if ((i == 3) && (j == 4) && (k == 2) && (l == 1))
    {t = ten.Transpose(1,4); t = t.Transpose(2,3); t = t.Transpose(3,4);}
  else if ((i == 4) && (j == 1) && (k == 2) && (l == 3))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3); t = t.Transpose(3,4);}
  else if ((i == 4) && (j == 1) && (k == 3) && (l == 2))
    {t = ten.Transpose(1,2); t = t.Transpose(2,4);}
  else if ((i == 4) && (j == 2) && (k == 1) && (l == 3))
    {t = ten.Transpose(1,3); t = t.Transpose(3,4);}
  else if ((i == 4) && (j == 2) && (k == 3) && (l == 1))
    {t = ten.Transpose(1,4);}
  else if ((i == 4) && (j == 3) && (k == 1) && (l == 2))
    {t = ten.Transpose(1,3); t = t.Transpose(2,4); t = t.Transpose(3,4);}
  else if ((i == 4) && (j == 3) && (k == 2) && (l == 1))
    {t = ten.Transpose(1,4); t = t.Transpose(2,3);}
  else
    return 0.;
  return   t.at(0,1,2,3)
          -t.at(0,1,3,2)
          -t.at(0,2,1,3)
          +t.at(0,2,3,1)
          +t.at(0,3,1,2)
          -t.at(0,3,2,1)

          -t.at(1,0,2,3)
          +t.at(1,0,3,2)
          +t.at(1,2,0,3)
          -t.at(1,2,3,0)
          -t.at(1,3,0,2)
          +t.at(1,3,2,0)

          +t.at(2,0,1,3)
          -t.at(2,0,3,1)
          -t.at(2,1,0,3)
          +t.at(2,1,3,0)
          +t.at(2,3,0,1)
          -t.at(2,3,1,0)

          -t.at(3,0,1,2)
          +t.at(3,0,2,1)
          +t.at(3,1,0,2)
          -t.at(3,1,2,0)
          -t.at(3,2,0,1)
          +t.at(3,2,1,0);
}


#endif
