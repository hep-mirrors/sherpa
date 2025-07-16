#ifndef ATOOLS_Math_Tensor_Contractions_H
#define ATOOLS_Math_Tensor_Contractions_H

/********************************************************************************************
******                                                                                 ******
******       tensor with itself                                                        ******
******                                                                                 ******
********************************************************************************************/

template<typename Scalar>
Scalar
ATOOLS::Contraction(const Lorentz_Ten2<Scalar>& t) {
  // c = t^{alpha,beta} g_{alpha,beta}
  return t.at(0,0)-t.at(1,1)-t.at(2,2)-t.at(3,3);
}

template<typename Scalar>
Vec4<Scalar>
ATOOLS::Contraction(const Lorentz_Ten3<Scalar>& ten, int i, int j) {
  Lorentz_Ten3<Scalar> t;
  // v^{mu} = t{mu,alpha,alpha'} g_{alpha,alpha'}
  if      (((i==2) && (j==3)) || ((j==2) && (i==3)))
    t = ten;
  // v^{mu} = t{alpha,mu,alpha'} g_{alpha,alpha'}
  else if (((i==1) && (j==3)) || ((j==1) && (i==3)))
    t = ten.Transpose(1,2);
  // v^{mu} = t{alpha,alpha',mu} g_{alpha,alpha'}
  else if (((i==1) && (j==2)) || ((j==1) && (i==2)))
    t = ten.Transpose(1,3);
  else
    return Vec4<Scalar>(0.,0.,0.,0.);
  return Vec4<Scalar>(t.at(0,0,0)-t.at(0,1,1)-t.at(0,2,2)-t.at(0,3,3),
                      t.at(1,0,0)-t.at(1,1,1)-t.at(1,2,2)-t.at(1,3,3),
                      t.at(2,0,0)-t.at(2,1,1)-t.at(2,2,2)-t.at(2,3,3),
                      t.at(3,0,0)-t.at(3,1,1)-t.at(3,2,2)-t.at(3,3,3));
}

template<typename Scalar>
Lorentz_Ten2<Scalar>
ATOOLS::Contraction(const Lorentz_Ten4<Scalar>& ten, int i, int j) {
  Lorentz_Ten4<Scalar> t;
  // T^{mu,nu} = t^{mu,nu,alpha,alpha'} g_{alpha,alpha'}
  if      (((i==3) && (j==4)) || ((i==4) && (j==3)))
     t = ten;
  // T^{mu,nu} = t^{mu,alpha,nu,alpha'} g_{alpha,alpha'}
  else if (((i==2) && (j==4)) || ((i==4) && (j==2)))
     t = ten.Transpose(2,3);
  // T^{mu,nu} = t^{alpha,mu,nu,alpha'} g_{alpha,alpha'}
  else if (((i==1) && (j==4)) || ((i==4) && (j==1)))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3);}
  // T^{mu,nu} = t^{mu,alpha,alpha',nu} g_{alpha,alpha'}
  else if (((i==2) && (j==3)) || ((i==3) && (j==2)))
     t = ten.Transpose(2,4);
  // T^{mu,nu} = t^{alpha,mu,alpha',nu} g_{alpha,alpha'}
  else if (((i==1) && (j==3)) || ((i==3) && (j==1)))
    {t = ten.Transpose(1,2); t = t.Transpose(2,4);}
  // T^{mu,nu} = t^{alpha,alpha',mu,nu} g_{alpha,alpha'}
  else if (((i==1) && (j==2)) || ((i==2) && (j==1)))
    {t = ten.Transpose(1,3); t = t.Transpose(2,4);}
  else
    return Lorentz_Ten2<Scalar>();
  return Lorentz_Ten2<Scalar>(t.at(0,0,0,0)-t.at(0,0,1,1)-t.at(0,0,2,2)-t.at(0,0,3,3),
                              t.at(1,0,0,0)-t.at(1,0,1,1)-t.at(1,0,2,2)-t.at(1,0,3,3),
                              t.at(2,0,0,0)-t.at(2,0,1,1)-t.at(2,0,2,2)-t.at(2,0,3,3),
                              t.at(3,0,0,0)-t.at(3,0,1,1)-t.at(3,0,2,2)-t.at(3,0,3,3),
                              t.at(0,1,0,0)-t.at(0,1,1,1)-t.at(0,1,2,2)-t.at(0,1,3,3),
                              t.at(1,1,0,0)-t.at(1,1,1,1)-t.at(1,1,2,2)-t.at(1,1,3,3),
                              t.at(2,1,0,0)-t.at(2,1,1,1)-t.at(2,1,2,2)-t.at(2,1,3,3),
                              t.at(3,1,0,0)-t.at(3,1,1,1)-t.at(3,1,2,2)-t.at(3,1,3,3),
                              t.at(0,2,0,0)-t.at(0,2,1,1)-t.at(0,2,2,2)-t.at(0,2,3,3),
                              t.at(1,2,0,0)-t.at(1,2,1,1)-t.at(1,2,2,2)-t.at(1,2,3,3),
                              t.at(2,2,0,0)-t.at(2,2,1,1)-t.at(2,2,2,2)-t.at(2,2,3,3),
                              t.at(3,2,0,0)-t.at(3,2,1,1)-t.at(3,2,2,2)-t.at(3,2,3,3),
                              t.at(0,3,0,0)-t.at(0,3,1,1)-t.at(0,3,2,2)-t.at(0,3,3,3),
                              t.at(1,3,0,0)-t.at(1,3,1,1)-t.at(1,3,2,2)-t.at(1,3,3,3),
                              t.at(2,3,0,0)-t.at(2,3,1,1)-t.at(2,3,2,2)-t.at(2,3,3,3),
                              t.at(3,3,0,0)-t.at(3,3,1,1)-t.at(3,3,2,2)-t.at(3,3,3,3));
}

template<typename Scalar>
Scalar
ATOOLS::Contraction(const Lorentz_Ten4<Scalar>& ten, int i, int j, int k, int l) {
  Lorentz_Ten4<Scalar> t;
  // c = t^{alpha,alpha',beta,beta'} g_{alpha,alpha'} g_{beta,beta'}
  if      (((((i==1) && (j==2)) || ((i==2) && (j==1)))  && 
            (((k==3) && (l==4)) || ((k==4) && (l==3)))) ||
           ((((i==3) && (j==4)) || ((i==4) && (j==3)))  && 
            (((k==1) && (l==2)) || ((k==2) && (l==1)))))
    t = ten;
  // c = t^{alpha,beta,alpha',beta'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((((i==1) && (j==3)) || ((i==3) && (j==1)))  && 
            (((k==2) && (l==4)) || ((k==4) && (l==2)))) ||
           ((((i==2) && (j==4)) || ((i==4) && (j==2)))  && 
            (((k==1) && (l==3)) || ((k==3) && (l==1)))))
    t = ten.Transpose(2,3);
  // c = t^{alpha,beta,beta',alpha'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((((i==1) && (j==4)) || ((i==4) && (j==1)))  && 
            (((k==2) && (l==3)) || ((k==3) && (l==2)))) ||
           ((((i==2) && (j==3)) || ((i==3) && (j==2)))  && 
            (((k==1) && (l==4)) || ((k==4) && (l==1)))))
    t = ten.Transpose(2,4);
  else
    return Scalar(0.);
  return Scalar( t.at(0,0,0,0)-t.at(0,0,1,1)-t.at(0,0,2,2)-t.at(0,0,3,3)
                -t.at(1,1,0,0)+t.at(1,1,1,1)+t.at(1,1,2,2)+t.at(1,1,3,3)
                -t.at(2,2,0,0)+t.at(2,2,1,1)+t.at(2,2,2,2)+t.at(2,2,3,3)
                -t.at(3,3,0,0)+t.at(3,3,1,1)+t.at(3,3,2,2)+t.at(3,3,3,3));
}

/********************************************************************************************
******                                                                                 ******
******       tensor with vector                                                        ******
******                                                                                 ******
********************************************************************************************/

template<typename Scal1, typename Scal2> 
Vec4<PROMOTE(Scal1,Scal2)> 
ATOOLS::Contraction(const Lorentz_Ten2<Scal1>& ten, int i, const Vec4<Scal2>& p) {
  Lorentz_Ten2<Scal1> t;
  if ((i==1) || (i==2)) {
    // v^{mu} = t^{alpha,mu} p^{beta} g_{alpha,beta}
    if      (i==1)
      t = ten.Transpose();
    // v^{mu} = t^{mu,alpha} p^{beta} g_{alpha,beta}
    else if (i==2)
      t = ten;
    return Vec4<PROMOTE(Scal1,Scal2)>(
                        t.at(0,0)*p[0]-t.at(0,1)*p[1]-t.at(0,2)*p[2]-t.at(0,3)*p[3],
                        t.at(1,0)*p[0]-t.at(1,1)*p[1]-t.at(1,2)*p[2]-t.at(1,3)*p[3],
                        t.at(2,0)*p[0]-t.at(2,1)*p[1]-t.at(2,2)*p[2]-t.at(2,3)*p[3],
                        t.at(3,0)*p[0]-t.at(3,1)*p[1]-t.at(3,2)*p[2]-t.at(3,3)*p[3]);
  }
  return Vec4<PROMOTE(Scal1,Scal2)>();
}

template<typename Scal1, typename Scal2, typename Scal3> 
PROMOTE(PROMOTE(Scal1,Scal2),Scal3)
ATOOLS::Contraction(const Lorentz_Ten2<Scal1>& t,
                    const Vec4<Scal2>& p, const Vec4<Scal3>& q) {
  // c = t^{alpha,beta} p^{alpha'} q^{beta'} g_{alpha',alpha} g_{beta,beta'}
  return  q[0]*(t.at(0,0)*p[0]-t.at(1,0)*p[1]-t.at(2,0)*p[2]-t.at(3,0)*p[3])
         -q[1]*(t.at(0,1)*p[0]-t.at(1,1)*p[1]-t.at(2,1)*p[2]-t.at(3,1)*p[3])
         -q[2]*(t.at(0,2)*p[0]-t.at(1,2)*p[1]-t.at(2,2)*p[2]-t.at(3,2)*p[3])
         -q[3]*(t.at(0,3)*p[0]-t.at(1,3)*p[1]-t.at(2,3)*p[2]-t.at(3,3)*p[3]);
}

template<typename Scal1, typename Scal2>
Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten3<Scal1>& ten, int i, const Vec4<Scal2>& v) {
  Lorentz_Ten3<Scal1> t;
  // c^{mu,nu} = t^{alpha,mu,nu} g_{alpha,alpha'} v^{alpha'}
  if      (i==1)
    {t = ten.Transpose(1,2); t = t.Transpose(2,3);}
  // c^{mu,nu} = t^{mu,alpha,nu} g_{alpha,alpha'} v^{alpha'}
  else if (i==2)
     t = ten.Transpose(2,3);
  // c^{mu,nu} = t^{mu,nu,alpha} g_{alpha,alpha'} v^{alpha'}
  else if (i==3)
     t = ten;
  else
    return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>();
  return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>(
            t.at(0,0,0)*v[0]-t.at(0,0,1)*v[1]-t.at(0,0,2)*v[2]-t.at(0,0,3)*v[3],
            t.at(1,0,0)*v[0]-t.at(1,0,1)*v[1]-t.at(1,0,2)*v[2]-t.at(1,0,3)*v[3],
            t.at(2,0,0)*v[0]-t.at(2,0,1)*v[1]-t.at(2,0,2)*v[2]-t.at(2,0,3)*v[3],
            t.at(3,0,0)*v[0]-t.at(3,0,1)*v[1]-t.at(3,0,2)*v[2]-t.at(3,0,3)*v[3],
            t.at(0,1,0)*v[0]-t.at(0,1,1)*v[1]-t.at(0,1,2)*v[2]-t.at(0,1,3)*v[3],
            t.at(1,1,0)*v[0]-t.at(1,1,1)*v[1]-t.at(1,1,2)*v[2]-t.at(1,1,3)*v[3],
            t.at(2,1,0)*v[0]-t.at(2,1,1)*v[1]-t.at(2,1,2)*v[2]-t.at(2,1,3)*v[3],
            t.at(3,1,0)*v[0]-t.at(3,1,1)*v[1]-t.at(3,1,2)*v[2]-t.at(3,1,3)*v[3],
            t.at(0,2,0)*v[0]-t.at(0,2,1)*v[1]-t.at(0,2,2)*v[2]-t.at(0,2,3)*v[3],
            t.at(1,2,0)*v[0]-t.at(1,2,1)*v[1]-t.at(1,2,2)*v[2]-t.at(1,2,3)*v[3],
            t.at(2,2,0)*v[0]-t.at(2,2,1)*v[1]-t.at(2,2,2)*v[2]-t.at(2,2,3)*v[3],
            t.at(3,2,0)*v[0]-t.at(3,2,1)*v[1]-t.at(3,2,2)*v[2]-t.at(3,2,3)*v[3],
            t.at(0,3,0)*v[0]-t.at(0,3,1)*v[1]-t.at(0,3,2)*v[2]-t.at(0,3,3)*v[3],
            t.at(1,3,0)*v[0]-t.at(1,3,1)*v[1]-t.at(1,3,2)*v[2]-t.at(1,3,3)*v[3],
            t.at(2,3,0)*v[0]-t.at(2,3,1)*v[1]-t.at(2,3,2)*v[2]-t.at(2,3,3)*v[3],
            t.at(3,3,0)*v[0]-t.at(3,3,1)*v[1]-t.at(3,3,2)*v[2]-t.at(3,3,3)*v[3]);
}

template<typename Scal1, typename Scal2>
Lorentz_Ten3<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, const Vec4<Scal2>& v) {
  Lorentz_Ten4<Scal1> t;
  // c^{mu,nu,rho} = t^{alpha,mu,nu,rho} g_{alpha,alpha'} v^{alpha'}
  if      (i==1)
    {t = ten.Transpose(1,2); t = t.Transpose(2,3); t = t.Transpose(3,4);}
  // c^{mu,nu,rho} = t^{mu,alpha,nu,rho} g_{alpha,alpha'} v^{alpha'}
  else if (i==2)
    {t = ten.Transpose(2,3); t = t.Transpose(3,4);}
  // c^{mu,nu,rho} = t^{mu,nu,alpha,rho} g_{alpha,alpha'} v^{alpha'}
  else if (i==3)
     t = ten.Transpose(3,4);
  // c^{mu,nu,rho} = t^{mu,nu,rho,alpha} g_{alpha,alpha'} v^{alpha'}
  else if (i==4)
     t = ten;
  else
    return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>();
  return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>(
    t.at(0,0,0,0)*v[0]-t.at(0,0,0,1)*v[1]-t.at(0,0,0,2)*v[2]-t.at(0,0,0,3)*v[3],
    t.at(1,0,0,0)*v[0]-t.at(1,0,0,1)*v[1]-t.at(1,0,0,2)*v[2]-t.at(1,0,0,3)*v[3],
    t.at(2,0,0,0)*v[0]-t.at(2,0,0,1)*v[1]-t.at(2,0,0,2)*v[2]-t.at(2,0,0,3)*v[3],
    t.at(3,0,0,0)*v[0]-t.at(3,0,0,1)*v[1]-t.at(3,0,0,2)*v[2]-t.at(3,0,0,3)*v[3],
    t.at(0,1,0,0)*v[0]-t.at(0,1,0,1)*v[1]-t.at(0,1,0,2)*v[2]-t.at(0,1,0,3)*v[3],
    t.at(1,1,0,0)*v[0]-t.at(1,1,0,1)*v[1]-t.at(1,1,0,2)*v[2]-t.at(1,1,0,3)*v[3],
    t.at(2,1,0,0)*v[0]-t.at(2,1,0,1)*v[1]-t.at(2,1,0,2)*v[2]-t.at(2,1,0,3)*v[3],
    t.at(3,1,0,0)*v[0]-t.at(3,1,0,1)*v[1]-t.at(3,1,0,2)*v[2]-t.at(3,1,0,3)*v[3],
    t.at(0,2,0,0)*v[0]-t.at(0,2,0,1)*v[1]-t.at(0,2,0,2)*v[2]-t.at(0,2,0,3)*v[3],
    t.at(1,2,0,0)*v[0]-t.at(1,2,0,1)*v[1]-t.at(1,2,0,2)*v[2]-t.at(1,2,0,3)*v[3],
    t.at(2,2,0,0)*v[0]-t.at(2,2,0,1)*v[1]-t.at(2,2,0,2)*v[2]-t.at(2,2,0,3)*v[3],
    t.at(3,2,0,0)*v[0]-t.at(3,2,0,1)*v[1]-t.at(3,2,0,2)*v[2]-t.at(3,2,0,3)*v[3],
    t.at(0,3,0,0)*v[0]-t.at(0,3,0,1)*v[1]-t.at(0,3,0,2)*v[2]-t.at(0,3,0,3)*v[3],
    t.at(1,3,0,0)*v[0]-t.at(1,3,0,1)*v[1]-t.at(1,3,0,2)*v[2]-t.at(1,3,0,3)*v[3],
    t.at(2,3,0,0)*v[0]-t.at(2,3,0,1)*v[1]-t.at(2,3,0,2)*v[2]-t.at(2,3,0,3)*v[3],
    t.at(3,3,0,0)*v[0]-t.at(3,3,0,1)*v[1]-t.at(3,3,0,2)*v[2]-t.at(3,3,0,3)*v[3],

    t.at(0,0,1,0)*v[0]-t.at(0,0,1,1)*v[1]-t.at(0,0,1,2)*v[2]-t.at(0,0,1,3)*v[3],
    t.at(1,0,1,0)*v[0]-t.at(1,0,1,1)*v[1]-t.at(1,0,1,2)*v[2]-t.at(1,0,1,3)*v[3],
    t.at(2,0,1,0)*v[0]-t.at(2,0,1,1)*v[1]-t.at(2,0,1,2)*v[2]-t.at(2,0,1,3)*v[3],
    t.at(3,0,1,0)*v[0]-t.at(3,0,1,1)*v[1]-t.at(3,0,1,2)*v[2]-t.at(3,0,1,3)*v[3],
    t.at(0,1,1,0)*v[0]-t.at(0,1,1,1)*v[1]-t.at(0,1,1,2)*v[2]-t.at(0,1,1,3)*v[3],
    t.at(1,1,1,0)*v[0]-t.at(1,1,1,1)*v[1]-t.at(1,1,1,2)*v[2]-t.at(1,1,1,3)*v[3],
    t.at(2,1,1,0)*v[0]-t.at(2,1,1,1)*v[1]-t.at(2,1,1,2)*v[2]-t.at(2,1,1,3)*v[3],
    t.at(3,1,1,0)*v[0]-t.at(3,1,1,1)*v[1]-t.at(3,1,1,2)*v[2]-t.at(3,1,1,3)*v[3],
    t.at(0,2,1,0)*v[0]-t.at(0,2,1,1)*v[1]-t.at(0,2,1,2)*v[2]-t.at(0,2,1,3)*v[3],
    t.at(1,2,1,0)*v[0]-t.at(1,2,1,1)*v[1]-t.at(1,2,1,2)*v[2]-t.at(1,2,1,3)*v[3],
    t.at(2,2,1,0)*v[0]-t.at(2,2,1,1)*v[1]-t.at(2,2,1,2)*v[2]-t.at(2,2,1,3)*v[3],
    t.at(3,2,1,0)*v[0]-t.at(3,2,1,1)*v[1]-t.at(3,2,1,2)*v[2]-t.at(3,2,1,3)*v[3],
    t.at(0,3,1,0)*v[0]-t.at(0,3,1,1)*v[1]-t.at(0,3,1,2)*v[2]-t.at(0,3,1,3)*v[3],
    t.at(1,3,1,0)*v[0]-t.at(1,3,1,1)*v[1]-t.at(1,3,1,2)*v[2]-t.at(1,3,1,3)*v[3],
    t.at(2,3,1,0)*v[0]-t.at(2,3,1,1)*v[1]-t.at(2,3,1,2)*v[2]-t.at(2,3,1,3)*v[3],
    t.at(3,3,1,0)*v[0]-t.at(3,3,1,1)*v[1]-t.at(3,3,1,2)*v[2]-t.at(3,3,1,3)*v[3],

    t.at(0,0,2,0)*v[0]-t.at(0,0,2,1)*v[1]-t.at(0,0,2,2)*v[2]-t.at(0,0,2,3)*v[3],
    t.at(1,0,2,0)*v[0]-t.at(1,0,2,1)*v[1]-t.at(1,0,2,2)*v[2]-t.at(1,0,2,3)*v[3],
    t.at(2,0,2,0)*v[0]-t.at(2,0,2,1)*v[1]-t.at(2,0,2,2)*v[2]-t.at(2,0,2,3)*v[3],
    t.at(3,0,2,0)*v[0]-t.at(3,0,2,1)*v[1]-t.at(3,0,2,2)*v[2]-t.at(3,0,2,3)*v[3],
    t.at(0,1,2,0)*v[0]-t.at(0,1,2,1)*v[1]-t.at(0,1,2,2)*v[2]-t.at(0,1,2,3)*v[3],
    t.at(1,1,2,0)*v[0]-t.at(1,1,2,1)*v[1]-t.at(1,1,2,2)*v[2]-t.at(1,1,2,3)*v[3],
    t.at(2,1,2,0)*v[0]-t.at(2,1,2,1)*v[1]-t.at(2,1,2,2)*v[2]-t.at(2,1,2,3)*v[3],
    t.at(3,1,2,0)*v[0]-t.at(3,1,2,1)*v[1]-t.at(3,1,2,2)*v[2]-t.at(3,1,2,3)*v[3],
    t.at(0,2,2,0)*v[0]-t.at(0,2,2,1)*v[1]-t.at(0,2,2,2)*v[2]-t.at(0,2,2,3)*v[3],
    t.at(1,2,2,0)*v[0]-t.at(1,2,2,1)*v[1]-t.at(1,2,2,2)*v[2]-t.at(1,2,2,3)*v[3],
    t.at(2,2,2,0)*v[0]-t.at(2,2,2,1)*v[1]-t.at(2,2,2,2)*v[2]-t.at(2,2,2,3)*v[3],
    t.at(3,2,2,0)*v[0]-t.at(3,2,2,1)*v[1]-t.at(3,2,2,2)*v[2]-t.at(3,2,2,3)*v[3],
    t.at(0,3,2,0)*v[0]-t.at(0,3,2,1)*v[1]-t.at(0,3,2,2)*v[2]-t.at(0,3,2,3)*v[3],
    t.at(1,3,2,0)*v[0]-t.at(1,3,2,1)*v[1]-t.at(1,3,2,2)*v[2]-t.at(1,3,2,3)*v[3],
    t.at(2,3,2,0)*v[0]-t.at(2,3,2,1)*v[1]-t.at(2,3,2,2)*v[2]-t.at(2,3,2,3)*v[3],
    t.at(3,3,2,0)*v[0]-t.at(3,3,2,1)*v[1]-t.at(3,3,2,2)*v[2]-t.at(3,3,2,3)*v[3],

    t.at(0,0,3,0)*v[0]-t.at(0,0,3,1)*v[1]-t.at(0,0,3,2)*v[2]-t.at(0,0,3,3)*v[3],
    t.at(1,0,3,0)*v[0]-t.at(1,0,3,1)*v[1]-t.at(1,0,3,2)*v[2]-t.at(1,0,3,3)*v[3],
    t.at(2,0,3,0)*v[0]-t.at(2,0,3,1)*v[1]-t.at(2,0,3,2)*v[2]-t.at(2,0,3,3)*v[3],
    t.at(3,0,3,0)*v[0]-t.at(3,0,3,1)*v[1]-t.at(3,0,3,2)*v[2]-t.at(3,0,3,3)*v[3],
    t.at(0,1,3,0)*v[0]-t.at(0,1,3,1)*v[1]-t.at(0,1,3,2)*v[2]-t.at(0,1,3,3)*v[3],
    t.at(1,1,3,0)*v[0]-t.at(1,1,3,1)*v[1]-t.at(1,1,3,2)*v[2]-t.at(1,1,3,3)*v[3],
    t.at(2,1,3,0)*v[0]-t.at(2,1,3,1)*v[1]-t.at(2,1,3,2)*v[2]-t.at(2,1,3,3)*v[3],
    t.at(3,1,3,0)*v[0]-t.at(3,1,3,1)*v[1]-t.at(3,1,3,2)*v[2]-t.at(3,1,3,3)*v[3],
    t.at(0,2,3,0)*v[0]-t.at(0,2,3,1)*v[1]-t.at(0,2,3,2)*v[2]-t.at(0,2,3,3)*v[3],
    t.at(1,2,3,0)*v[0]-t.at(1,2,3,1)*v[1]-t.at(1,2,3,2)*v[2]-t.at(1,2,3,3)*v[3],
    t.at(2,2,3,0)*v[0]-t.at(2,2,3,1)*v[1]-t.at(2,2,3,2)*v[2]-t.at(2,2,3,3)*v[3],
    t.at(3,2,3,0)*v[0]-t.at(3,2,3,1)*v[1]-t.at(3,2,3,2)*v[2]-t.at(3,2,3,3)*v[3],
    t.at(0,3,3,0)*v[0]-t.at(0,3,3,1)*v[1]-t.at(0,3,3,2)*v[2]-t.at(0,3,3,3)*v[3],
    t.at(1,3,3,0)*v[0]-t.at(1,3,3,1)*v[1]-t.at(1,3,3,2)*v[2]-t.at(1,3,3,3)*v[3],
    t.at(2,3,3,0)*v[0]-t.at(2,3,3,1)*v[1]-t.at(2,3,3,2)*v[2]-t.at(2,3,3,3)*v[3],
    t.at(3,3,3,0)*v[0]-t.at(3,3,3,1)*v[1]-t.at(3,3,3,2)*v[2]-t.at(3,3,3,3)*v[3]);
}

/********************************************************************************************
******                                                                                 ******
******       tensor with tensor                                                        ******
******                                                                                 ******
********************************************************************************************/

template<typename Scal1, typename Scal2> 
PROMOTE(Scal1,Scal2) 
ATOOLS::Contraction(const Lorentz_Ten2<Scal1>& ten, int i, int j,
                    const Lorentz_Ten2<Scal2>& sen, int k, int l) {
  Lorentz_Ten2<Scal1> t; Lorentz_Ten2<Scal2> s;
  // c = t^{alpha,beta} s^{alpha',beta'} g_{alpha',alpha} g_{beta,beta'}
  if      (((i==1) && (j==2) && (k==1) && (l==2)) ||
           ((i==2) && (j==1) && (k==2) && (l==1)))
    {t = ten.Transpose(); s = sen;}
  // c = t^{alpha,beta} s^{beta',alpha'} g_{alpha',alpha} g_{beta,beta'}
  else if (((i==1) && (j==2) && (k==2) && (l==1)) ||
           ((i==2) && (j==1) && (k==1) && (l==2)))
    {t = ten;             s = sen;}
  else 
    return PROMOTE(Scal1,Scal2)();
  return PROMOTE(Scal1,Scal2)
    (t.at(0,0)*s.at(0,0)-t.at(1,0)*s.at(0,1)-t.at(2,0)*s.at(0,2)-t.at(3,0)*s.at(0,3))
   -(t.at(0,1)*s.at(1,0)-t.at(1,1)*s.at(1,1)-t.at(2,1)*s.at(1,2)-t.at(3,1)*s.at(1,3))
   -(t.at(0,2)*s.at(2,0)-t.at(1,2)*s.at(2,1)-t.at(2,2)*s.at(2,2)-t.at(3,2)*s.at(2,3))
   -(t.at(0,3)*s.at(3,0)-t.at(1,3)*s.at(3,1)-t.at(2,3)*s.at(3,2)-t.at(3,3)*s.at(3,3));
}

template<typename Scal1, typename Scal2> 
Lorentz_Ten2<PROMOTE(Scal1,Scal2)> 
ATOOLS::Contraction(const Lorentz_Ten2<Scal1>& ten,int i,
                    const Lorentz_Ten2<Scal2>& sen, int j) {
  Lorentz_Ten2<Scal1> t; Lorentz_Ten2<Scal2> s;
  if (((i==1) || (i==2)) && ((j==1) || (j==2))) {
    // T^{mu,nu} = t^{alpha,mu} s^{beta,nu} g_{alpha,beta}
    if      ((i==1) && (j==1))
      {t = ten.Transpose(); s = sen;}
    // T^{mu,nu} = t^{alpha,mu} s^{nu,beta} g_{alpha,beta}
    else if ((i==1) && (j==2))
      {t = ten.Transpose(); s = sen.Transpose();}
    // T^{mu,nu} = t^{mu,alpha} s^{beta,nu} g_{alpha,beta}
    else if ((i==2) && (j==1))
      {t = ten; s = sen;}
    // T^{mu,nu} = t^{mu,alpha} s^{nu,beta} g_{alpha,beta}
    else if ((i==2) && (j==2))
      {t = ten; s = sen.Transpose();}
    return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>(
      s.at(0,0)*t.at(0,0)-s.at(1,0)*t.at(0,1)-s.at(2,0)*t.at(0,2)-s.at(3,0)*t.at(0,3),
      s.at(0,0)*t.at(1,0)-s.at(1,0)*t.at(1,1)-s.at(2,0)*t.at(1,2)-s.at(3,0)*t.at(1,3),
      s.at(0,0)*t.at(2,0)-s.at(1,0)*t.at(2,1)-s.at(2,0)*t.at(2,2)-s.at(3,0)*t.at(2,3),
      s.at(0,0)*t.at(3,0)-s.at(1,0)*t.at(3,1)-s.at(2,0)*t.at(3,2)-s.at(3,0)*t.at(3,3),
      s.at(0,1)*t.at(0,0)-s.at(1,1)*t.at(0,1)-s.at(2,1)*t.at(0,2)-s.at(3,1)*t.at(0,3),
      s.at(0,1)*t.at(1,0)-s.at(1,1)*t.at(1,1)-s.at(2,1)*t.at(1,2)-s.at(3,1)*t.at(1,3),
      s.at(0,1)*t.at(2,0)-s.at(1,1)*t.at(2,1)-s.at(2,1)*t.at(2,2)-s.at(3,1)*t.at(2,3),
      s.at(0,1)*t.at(3,0)-s.at(1,1)*t.at(3,1)-s.at(2,1)*t.at(3,2)-s.at(3,1)*t.at(3,3),
      s.at(0,2)*t.at(0,0)-s.at(1,2)*t.at(0,1)-s.at(2,2)*t.at(0,2)-s.at(3,2)*t.at(0,3),
      s.at(0,2)*t.at(1,0)-s.at(1,2)*t.at(1,1)-s.at(2,2)*t.at(1,2)-s.at(3,2)*t.at(1,3),
      s.at(0,2)*t.at(2,0)-s.at(1,2)*t.at(2,1)-s.at(2,2)*t.at(2,2)-s.at(3,2)*t.at(2,3),
      s.at(0,2)*t.at(3,0)-s.at(1,2)*t.at(3,1)-s.at(2,2)*t.at(3,2)-s.at(3,2)*t.at(3,3),
      s.at(0,3)*t.at(0,0)-s.at(1,3)*t.at(0,1)-s.at(2,3)*t.at(0,2)-s.at(3,3)*t.at(0,3),
      s.at(0,3)*t.at(1,0)-s.at(1,3)*t.at(1,1)-s.at(2,3)*t.at(1,2)-s.at(3,3)*t.at(1,3),
      s.at(0,3)*t.at(2,0)-s.at(1,3)*t.at(2,1)-s.at(2,3)*t.at(2,2)-s.at(3,3)*t.at(2,3),
      s.at(0,3)*t.at(3,0)-s.at(1,3)*t.at(3,1)-s.at(2,3)*t.at(3,2)-s.at(3,3)*t.at(3,3));
  }
  return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>();
}

template<typename Scal1, typename Scal2> 
Vec4<PROMOTE(Scal1,Scal2)> 
ATOOLS::Contraction(const Lorentz_Ten3<Scal1>& ten, int i, int j,
                    const Lorentz_Ten2<Scal2>& sen, int k, int l) {
  Lorentz_Ten3<Scal1> t; Lorentz_Ten2<Scal2> s;
  // v^{mu} = t^{mu,alpha,beta} s^{alpha',beta'} g_{alpha,alpha'} g_{beta,beta'}
  if      (((i==2) && (j==3) && (k==1) && (l==2)) || ((i==3) && (j==2) && (k==2) && (l==1)))
  // v^{mu} = t^{alpha,mu,beta} s^{alpha',beta'} g_{alpha,alpha'} g_{beta,beta'}
    {t = ten; s = sen;}
  else if (((i==1) && (j==3) && (k==1) && (l==2)) || ((i==3) && (j==1) && (k==2) && (l==1)))
    {t = ten.Transpose(1,2); s = sen;}
  // v^{mu} = t^{alpha,beta,mu} s^{alpha',beta'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((i==1) && (j==2) && (k==1) && (l==2)) || ((i==2) && (j==1) && (k==2) && (l==1)))
    {t = ten.Transpose(2,3); t = t.Transpose(1,2); s = sen;}
  // v^{mu} = t^{mu,alpha,beta} s^{beta',alpha'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((i==2) && (j==3) && (k==2) && (l==1)) || ((i==3) && (j==2) && (k==1) && (l==2)))
    {t = ten; s = sen.Transpose();}
  // v^{mu} = t^{alpha,mu,beta} s^{beta',alpha'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((i==1) && (j==3) && (k==2) && (l==1)) || ((i==3) && (j==1) && (k==1) && (l==2)))
    {t = ten.Transpose(1,2); s = sen.Transpose();}
  // v^{mu} = t^{alpha,beta,mu} s^{beta',alpha'} g_{alpha,alpha'} g_{beta,beta'}
  else if (((i==1) && (j==2) && (k==2) && (l==1)) || ((i==2) && (j==1) && (k==1) && (l==2)))
    {t = ten.Transpose(2,3); t = t.Transpose(1,2); s = sen.Transpose();}
  else
    return Vec4<PROMOTE(Scal1,Scal2)>(0.,0.,0.,0.);
  return Vec4<PROMOTE(Scal1,Scal2)>(
    t.at(0,0,0)*s.at(0,0)-t.at(0,1,0)*s.at(1,0)-t.at(0,2,0)*s.at(2,0)-t.at(0,3,0)*s.at(3,0)
   -t.at(0,0,1)*s.at(0,1)+t.at(0,1,1)*s.at(1,1)+t.at(0,2,1)*s.at(2,1)+t.at(0,3,1)*s.at(3,1)
   -t.at(0,0,2)*s.at(0,2)+t.at(0,1,2)*s.at(1,2)+t.at(0,2,2)*s.at(2,2)+t.at(0,3,2)*s.at(3,2)
   -t.at(0,0,3)*s.at(0,3)+t.at(0,1,3)*s.at(1,3)+t.at(0,2,3)*s.at(2,3)+t.at(0,3,3)*s.at(3,3),
    t.at(1,0,0)*s.at(0,0)-t.at(1,1,0)*s.at(1,0)-t.at(1,2,0)*s.at(2,0)-t.at(1,3,0)*s.at(3,0)
   -t.at(1,0,1)*s.at(0,1)+t.at(1,1,1)*s.at(1,1)+t.at(1,2,1)*s.at(2,1)+t.at(1,3,1)*s.at(3,1)
   -t.at(1,0,2)*s.at(0,2)+t.at(1,1,2)*s.at(1,2)+t.at(1,2,2)*s.at(2,2)+t.at(1,3,2)*s.at(3,2)
   -t.at(1,0,3)*s.at(0,3)+t.at(1,1,3)*s.at(1,3)+t.at(1,2,3)*s.at(2,3)+t.at(1,3,3)*s.at(3,3),
    t.at(2,0,0)*s.at(0,0)-t.at(2,1,0)*s.at(1,0)-t.at(2,2,0)*s.at(2,0)-t.at(2,3,0)*s.at(3,0)
   -t.at(2,0,1)*s.at(0,1)+t.at(2,1,1)*s.at(1,1)+t.at(2,2,1)*s.at(2,1)+t.at(2,3,1)*s.at(3,1)
   -t.at(2,0,2)*s.at(0,2)+t.at(2,1,2)*s.at(1,2)+t.at(2,2,2)*s.at(2,2)+t.at(2,3,2)*s.at(3,2)
   -t.at(2,0,3)*s.at(0,3)+t.at(2,1,3)*s.at(1,3)+t.at(2,2,3)*s.at(2,3)+t.at(2,3,3)*s.at(3,3),
    t.at(3,0,0)*s.at(0,0)-t.at(3,1,0)*s.at(1,0)-t.at(3,2,0)*s.at(2,0)-t.at(3,3,0)*s.at(3,0)
   -t.at(3,0,1)*s.at(0,1)+t.at(3,1,1)*s.at(1,1)+t.at(3,2,1)*s.at(2,1)+t.at(3,3,1)*s.at(3,1)
   -t.at(3,0,2)*s.at(0,2)+t.at(3,1,2)*s.at(1,2)+t.at(3,2,2)*s.at(2,2)+t.at(3,3,2)*s.at(3,2)
   -t.at(3,0,3)*s.at(0,3)+t.at(3,1,3)*s.at(1,3)+t.at(3,2,3)*s.at(2,3)+t.at(3,3,3)*s.at(3,3));
}

template<typename Scal1, typename Scal2> 
Lorentz_Ten3<PROMOTE(Scal1,Scal2)> 
ATOOLS::Contraction(const Lorentz_Ten3<Scal1>& ten, int i,
                    const Lorentz_Ten2<Scal2>& sen, int j) {
  Lorentz_Ten3<Scal1> t; Lorentz_Ten2<Scal2> s;
  // T^{mu,nu,rho} = t^{mu,nu,alpha} s^{rho,alpha'} g_{alpha,alpha'}
  if      ((i==3) && (j==2))
    {t = ten; s = sen;}
  else if ((i==2) && (j==2))
    {t = ten.Transpose(2,3); s = sen;}
  else if ((i==1) && (j==2))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3); s = sen;}
  else if ((i==3) && (j==1))
    {t = ten; s = sen.Transpose();}
  else if ((i==2) && (j==1))
    {t = ten.Transpose(2,3); s = sen.Transpose();}
  else if ((i==1) && (j==1))
    {t = ten.Transpose(1,2); t = t.Transpose(2,3); s = sen.Transpose();}
  else
    return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>();
  return Lorentz_Ten3<PROMOTE(Scal1,Scal2)>( 
              t.at(0,0,0)*s.at(0,0)-t.at(0,0,1)*s.at(0,1)
                            -t.at(0,0,2)*s.at(0,2)-t.at(0,0,3)*s.at(0,3),
              t.at(1,0,0)*s.at(0,0)-t.at(1,0,1)*s.at(0,1)
                            -t.at(1,0,2)*s.at(0,2)-t.at(1,0,3)*s.at(0,3),
              t.at(2,0,0)*s.at(0,0)-t.at(2,0,1)*s.at(0,1)
                            -t.at(2,0,2)*s.at(0,2)-t.at(2,0,3)*s.at(0,3),
              t.at(3,0,0)*s.at(0,0)-t.at(3,0,1)*s.at(0,1)
                            -t.at(3,0,2)*s.at(0,2)-t.at(3,0,3)*s.at(0,3),
              t.at(0,1,0)*s.at(0,0)-t.at(0,1,1)*s.at(0,1)
                            -t.at(0,1,2)*s.at(0,2)-t.at(0,1,3)*s.at(0,3),
              t.at(1,1,0)*s.at(0,0)-t.at(1,1,1)*s.at(0,1)
                            -t.at(1,1,2)*s.at(0,2)-t.at(1,1,3)*s.at(0,3),
              t.at(2,1,0)*s.at(0,0)-t.at(2,1,1)*s.at(0,1)
                            -t.at(2,1,2)*s.at(0,2)-t.at(2,1,3)*s.at(0,3),
              t.at(3,1,0)*s.at(0,0)-t.at(3,1,1)*s.at(0,1)
                            -t.at(3,1,2)*s.at(0,2)-t.at(3,1,3)*s.at(0,3),
              t.at(0,2,0)*s.at(0,0)-t.at(0,2,1)*s.at(0,1)
                            -t.at(0,2,2)*s.at(0,2)-t.at(0,2,3)*s.at(0,3),
              t.at(1,2,0)*s.at(0,0)-t.at(1,2,1)*s.at(0,1)
                            -t.at(1,2,2)*s.at(0,2)-t.at(1,2,3)*s.at(0,3),
              t.at(2,2,0)*s.at(0,0)-t.at(2,2,1)*s.at(0,1)
                            -t.at(2,2,2)*s.at(0,2)-t.at(2,2,3)*s.at(0,3),
              t.at(3,2,0)*s.at(0,0)-t.at(3,2,1)*s.at(0,1)
                            -t.at(3,2,2)*s.at(0,2)-t.at(3,2,3)*s.at(0,3),
              t.at(0,3,0)*s.at(0,0)-t.at(0,3,1)*s.at(0,1)
                            -t.at(0,3,2)*s.at(0,2)-t.at(0,3,3)*s.at(0,3),
              t.at(1,3,0)*s.at(0,0)-t.at(1,3,1)*s.at(0,1)
                            -t.at(1,3,2)*s.at(0,2)-t.at(1,3,3)*s.at(0,3),
              t.at(2,3,0)*s.at(0,0)-t.at(2,3,1)*s.at(0,1)
                            -t.at(2,3,2)*s.at(0,2)-t.at(2,3,3)*s.at(0,3),
              t.at(3,3,0)*s.at(0,0)-t.at(3,3,1)*s.at(0,1)
                            -t.at(3,3,2)*s.at(0,2)-t.at(3,3,3)*s.at(0,3),

              t.at(0,0,0)*s.at(1,0)-t.at(0,0,1)*s.at(1,1)
                            -t.at(0,0,2)*s.at(1,2)-t.at(0,0,3)*s.at(1,3),
              t.at(1,0,0)*s.at(1,0)-t.at(1,0,1)*s.at(1,1)
                            -t.at(1,0,2)*s.at(1,2)-t.at(1,0,3)*s.at(1,3),
              t.at(2,0,0)*s.at(1,0)-t.at(2,0,1)*s.at(1,1)
                            -t.at(2,0,2)*s.at(1,2)-t.at(2,0,3)*s.at(1,3),
              t.at(3,0,0)*s.at(1,0)-t.at(3,0,1)*s.at(1,1)
                            -t.at(3,0,2)*s.at(1,2)-t.at(3,0,3)*s.at(1,3),
              t.at(0,1,0)*s.at(1,0)-t.at(0,1,1)*s.at(1,1)
                            -t.at(0,1,2)*s.at(1,2)-t.at(0,1,3)*s.at(1,3),
              t.at(1,1,0)*s.at(1,0)-t.at(1,1,1)*s.at(1,1)
                            -t.at(1,1,2)*s.at(1,2)-t.at(1,1,3)*s.at(1,3),
              t.at(2,1,0)*s.at(1,0)-t.at(2,1,1)*s.at(1,1)
                            -t.at(2,1,2)*s.at(1,2)-t.at(2,1,3)*s.at(1,3),
              t.at(3,1,0)*s.at(1,0)-t.at(3,1,1)*s.at(1,1)
                            -t.at(3,1,2)*s.at(1,2)-t.at(3,1,3)*s.at(1,3),
              t.at(0,2,0)*s.at(1,0)-t.at(0,2,1)*s.at(1,1)
                            -t.at(0,2,2)*s.at(1,2)-t.at(0,2,3)*s.at(1,3),
              t.at(1,2,0)*s.at(1,0)-t.at(1,2,1)*s.at(1,1)
                            -t.at(1,2,2)*s.at(1,2)-t.at(1,2,3)*s.at(1,3),
              t.at(2,2,0)*s.at(1,0)-t.at(2,2,1)*s.at(1,1)
                            -t.at(2,2,2)*s.at(1,2)-t.at(2,2,3)*s.at(1,3),
              t.at(3,2,0)*s.at(1,0)-t.at(3,2,1)*s.at(1,1)
                            -t.at(3,2,2)*s.at(1,2)-t.at(3,2,3)*s.at(1,3),
              t.at(0,3,0)*s.at(1,0)-t.at(0,3,1)*s.at(1,1)
                            -t.at(0,3,2)*s.at(1,2)-t.at(0,3,3)*s.at(1,3),
              t.at(1,3,0)*s.at(1,0)-t.at(1,3,1)*s.at(1,1)
                            -t.at(1,3,2)*s.at(1,2)-t.at(1,3,3)*s.at(1,3),
              t.at(2,3,0)*s.at(1,0)-t.at(2,3,1)*s.at(1,1)
                            -t.at(2,3,2)*s.at(1,2)-t.at(2,3,3)*s.at(1,3),
              t.at(3,3,0)*s.at(1,0)-t.at(3,3,1)*s.at(1,1)
                            -t.at(3,3,2)*s.at(1,2)-t.at(3,3,3)*s.at(1,3),

              t.at(0,0,0)*s.at(2,0)-t.at(0,0,1)*s.at(2,1)
                            -t.at(0,0,2)*s.at(2,2)-t.at(0,0,3)*s.at(2,3),
              t.at(1,0,0)*s.at(2,0)-t.at(1,0,1)*s.at(2,1)
                            -t.at(1,0,2)*s.at(2,2)-t.at(1,0,3)*s.at(2,3),
              t.at(2,0,0)*s.at(2,0)-t.at(2,0,1)*s.at(2,1)
                            -t.at(2,0,2)*s.at(2,2)-t.at(2,0,3)*s.at(2,3),
              t.at(3,0,0)*s.at(2,0)-t.at(3,0,1)*s.at(2,1)
                            -t.at(3,0,2)*s.at(2,2)-t.at(3,0,3)*s.at(2,3),
              t.at(0,1,0)*s.at(2,0)-t.at(0,1,1)*s.at(2,1)
                            -t.at(0,1,2)*s.at(2,2)-t.at(0,1,3)*s.at(2,3),
              t.at(1,1,0)*s.at(2,0)-t.at(1,1,1)*s.at(2,1)
                            -t.at(1,1,2)*s.at(2,2)-t.at(1,1,3)*s.at(2,3),
              t.at(2,1,0)*s.at(2,0)-t.at(2,1,1)*s.at(2,1)
                            -t.at(2,1,2)*s.at(2,2)-t.at(2,1,3)*s.at(2,3),
              t.at(3,1,0)*s.at(2,0)-t.at(3,1,1)*s.at(2,1)
                            -t.at(3,1,2)*s.at(2,2)-t.at(3,1,3)*s.at(2,3),
              t.at(0,2,0)*s.at(2,0)-t.at(0,2,1)*s.at(2,1)
                            -t.at(0,2,2)*s.at(2,2)-t.at(0,2,3)*s.at(2,3),
              t.at(1,2,0)*s.at(2,0)-t.at(1,2,1)*s.at(2,1)
                            -t.at(1,2,2)*s.at(2,2)-t.at(1,2,3)*s.at(2,3),
              t.at(2,2,0)*s.at(2,0)-t.at(2,2,1)*s.at(2,1)
                            -t.at(2,2,2)*s.at(2,2)-t.at(2,2,3)*s.at(2,3),
              t.at(3,2,0)*s.at(2,0)-t.at(3,2,1)*s.at(2,1)
                            -t.at(3,2,2)*s.at(2,2)-t.at(3,2,3)*s.at(2,3),
              t.at(0,3,0)*s.at(2,0)-t.at(0,3,1)*s.at(2,1)
                            -t.at(0,3,2)*s.at(2,2)-t.at(0,3,3)*s.at(2,3),
              t.at(1,3,0)*s.at(2,0)-t.at(1,3,1)*s.at(2,1)
                            -t.at(1,3,2)*s.at(2,2)-t.at(1,3,3)*s.at(2,3),
              t.at(2,3,0)*s.at(2,0)-t.at(2,3,1)*s.at(2,1)
                            -t.at(2,3,2)*s.at(2,2)-t.at(2,3,3)*s.at(2,3),
              t.at(3,3,0)*s.at(2,0)-t.at(3,3,1)*s.at(2,1)
                            -t.at(3,3,2)*s.at(2,2)-t.at(3,3,3)*s.at(2,3),

              t.at(0,0,0)*s.at(3,0)-t.at(0,0,1)*s.at(3,1)
                            -t.at(0,0,2)*s.at(3,2)-t.at(0,0,3)*s.at(3,3),
              t.at(1,0,0)*s.at(3,0)-t.at(1,0,1)*s.at(3,1)
                            -t.at(1,0,2)*s.at(3,2)-t.at(1,0,3)*s.at(3,3),
              t.at(2,0,0)*s.at(3,0)-t.at(2,0,1)*s.at(3,1)
                            -t.at(2,0,2)*s.at(3,2)-t.at(2,0,3)*s.at(3,3),
              t.at(3,0,0)*s.at(3,0)-t.at(3,0,1)*s.at(3,1)
                            -t.at(3,0,2)*s.at(3,2)-t.at(3,0,3)*s.at(3,3),
              t.at(0,1,0)*s.at(3,0)-t.at(0,1,1)*s.at(3,1)
                            -t.at(0,1,2)*s.at(3,2)-t.at(0,1,3)*s.at(3,3),
              t.at(1,1,0)*s.at(3,0)-t.at(1,1,1)*s.at(3,1)
                            -t.at(1,1,2)*s.at(3,2)-t.at(1,1,3)*s.at(3,3),
              t.at(2,1,0)*s.at(3,0)-t.at(2,1,1)*s.at(3,1)
                            -t.at(2,1,2)*s.at(3,2)-t.at(2,1,3)*s.at(3,3),
              t.at(3,1,0)*s.at(3,0)-t.at(3,1,1)*s.at(3,1)
                            -t.at(3,1,2)*s.at(3,2)-t.at(3,1,3)*s.at(3,3),
              t.at(0,2,0)*s.at(3,0)-t.at(0,2,1)*s.at(3,1)
                            -t.at(0,2,2)*s.at(3,2)-t.at(0,2,3)*s.at(3,3),
              t.at(1,2,0)*s.at(3,0)-t.at(1,2,1)*s.at(3,1)
                            -t.at(1,2,2)*s.at(3,2)-t.at(1,2,3)*s.at(3,3),
              t.at(2,2,0)*s.at(3,0)-t.at(2,2,1)*s.at(3,1)
                            -t.at(2,2,2)*s.at(3,2)-t.at(2,2,3)*s.at(3,3),
              t.at(3,2,0)*s.at(3,0)-t.at(3,2,1)*s.at(3,1)
                            -t.at(3,2,2)*s.at(3,2)-t.at(3,2,3)*s.at(3,3),
              t.at(0,3,0)*s.at(3,0)-t.at(0,3,1)*s.at(3,1)
                            -t.at(0,3,2)*s.at(3,2)-t.at(0,3,3)*s.at(3,3),
              t.at(1,3,0)*s.at(3,0)-t.at(1,3,1)*s.at(3,1)
                            -t.at(1,3,2)*s.at(3,2)-t.at(1,3,3)*s.at(3,3),
              t.at(2,3,0)*s.at(3,0)-t.at(2,3,1)*s.at(3,1)
                            -t.at(2,3,2)*s.at(3,2)-t.at(2,3,3)*s.at(3,3),
              t.at(3,3,0)*s.at(3,0)-t.at(3,3,1)*s.at(3,1)
                            -t.at(3,3,2)*s.at(3,2)-t.at(3,3,3)*s.at(3,3));

}

template<typename Scal1, typename Scal2>
Vec4<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, int j, int k,
                    const Lorentz_Ten3<Scal2>& sen, int l, int m, int n) {
  // v^{mu} = t^{mu,alpha,beta,gamma} s^{alpha',beta',gamma'}
  //                              g_{alpha,alpha'} g_{beta,beta'} g_{gamma,gamma'}
  if ((i!=2) && (j!=3) && (k!=4) && (l!=2) && (m!=3) && (n!=4)) {
    THROW(fatal_error,"tensor contraction not specified for these indizes ...");
  }
  Vec4<PROMOTE(Scal1,Scal2)> end;
  for (unsigned short int I=0; I<4; ++I) {
    PROMOTE(Scal1,Scal2) res(ten.at(0,0,0,0)-ten.at(0,0,0,0));
    for (unsigned short int a=0; a<4; ++a)
      for (unsigned short int b=0; b<4; ++b)
        for (unsigned short int c=0; c<4; ++c) {
          double ga(-1.),gb(-1.),gc(-1.);
          if (a==0) ga=1.; if (b==0) gb=1.; if (c==0) gc=1.;
          res += ga*gb*gc*ten.at(I,a,b,c)*sen.at(a,b,c);
        }
    end[I] = res;
  }
  return end;
}

template<typename Scal1, typename Scal2>
Lorentz_Ten3<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, int j,
                    const Lorentz_Ten3<Scal2>& sen, int k, int l) {
  // T^{mu,nu,rho} = t^{mu,nu,alpha,beta} s^{rho,alpha',beta'}
  //                                g_{alpha,alpha'} g_{beta,beta'}
  if ((i!=3) && (j!=4) && (k!=3) && (l!=4)) {
    THROW(fatal_error,"tensor contraction not specified for these indizes ...");
  }
  Vec4<PROMOTE(Scal1,Scal2)> end;
  for (unsigned short int I=0; I<4; ++I)
    for (unsigned short int J=0; J<4; ++J)
      for (unsigned short int K=0; K<4; ++K) {
        PROMOTE(Scal1,Scal2) res(ten.at(0,0,0,0)-ten.at(0,0,0,0));
        for (unsigned short int a=0; a<4; ++a)
          for (unsigned short int b=0; b<4; ++b) {
            double ga(-1.),gb(-1.);
            if (a==0) ga=1.; if (b==0) gb=1.;
            res += ga*gb*ten.at(I,J,a,b)*sen.at(K,a,b);
          }
        end.at(I,J,K) = res;
      }
  return end;
}

template<typename Scal1, typename Scal2>
PROMOTE(Scal1,Scal2)
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, int j, int k, int l,
                    const Lorentz_Ten4<Scal2>& sen, int m, int n, int o, int p) {
  Lorentz_Ten4<Scal1> t; Lorentz_Ten4<Scal2> s;
  // c = t^{alpha,beta,gamma,delta} s^{alpha',beta',gamma',delta'}
  //            g_{alpha,alpha'} g_{beta,beta'} g_{gamma.gamma'} g_{delta,delta'}
  if       ((i==m) && (j==n) && (k==o) && (l==p))
    {t=ten; s=sen;}
  else if  ((i==m) && (j==n) && (k==p) && (l==o))
    {t=ten; s=sen.Transpose(3,4);}
  else if  ((i==m) && (j==o) && (k==n) && (l==p))
    {t=ten; s=sen.Transpose(2,3);}
  else if  ((i==m) && (j==o) && (k==p) && (l==n))
    {t=ten; s=sen.Transpose(3,4); s=s.Transpose(2,3);}
  else if  ((i==n) && (j==m) && (k==o) && (l==p))
    {t=ten; s=sen.Transpose(1,2);}
  else if  ((i==n) && (j==m) && (k==p) && (l==o))
    {t=ten; s=sen.Transpose(1,2); s=s.Transpose(3,4);}
  else if  ((i==n) && (j==o) && (k==m) && (l==p))
    {t=ten; s=sen.Transpose(1,3); s=s.Transpose(2,3);}
  else if  ((i==n) && (j==o) && (k==p) && (l==m))
    {t=ten; s=sen.Transpose(3,4); s=s.Transpose(2,3); s=s.Transpose(1,2);}
  else if  ((i==n) && (j==p) && (k==m) && (l==o))
    {t=ten; s=sen.Transpose(1,3); s=s.Transpose(2,4); s=s.Transpose(2,3);}
  else if  ((i==n) && (j==p) && (k==o) && (l==m))
    {t=ten; s=sen.Transpose(2,4); s=s.Transpose(1,2);}
  else if  ((i==o) && (j==m) && (k==n) && (l==p))
    {t=ten; s=sen.Transpose(1,2); s=s.Transpose(2,3);}
  else if  ((i==o) && (j==m) && (k==p) && (l==n))
    {t=ten; s=sen.Transpose(1,2); s=s.Transpose(2,3); s=s.Transpose(2,4);}
  else if  ((i==o) && (j==n) && (k==m) && (l==p))
    {t=ten; s=sen.Transpose(1,3);}
  else if  ((i==o) && (j==n) && (k==p) && (l==m))
    {t=ten; s=sen.Transpose(1,4); s=s.Transpose(3,4);}
  else if  ((i==o) && (j==p) && (k==m) && (l==n))
    {t=ten; s=sen.Transpose(1,3); s=s.Transpose(2,4);}
  else if  ((i==o) && (j==p) && (k==n) && (l==m))
    {t=ten; s=sen.Transpose(1,4); s=s.Transpose(2,3); s=s.Transpose(3,4);}
  else if  ((i==p) && (j==m) && (k==n) && (l==o))
    {t=ten; s=sen.Transpose(1,2); s=s.Transpose(2,3); s=s.Transpose(3,4);}
  else if  ((i==p) && (j==m) && (k==o) && (l==n))
    {t=ten; s=sen.Transpose(1,2); s=s.Transpose(2,4);}
  else if  ((i==p) && (j==n) && (k==m) && (l==o))
    {t=ten; s=sen.Transpose(1,3); s=s.Transpose(3,4);}
  else if  ((i==p) && (j==n) && (k==o) && (l==m))
    {t=ten; s=sen.Transpose(1,4);}
  else if  ((i==p) && (j==o) && (k==m) && (l==n))
    {t=ten; s=sen.Transpose(1,3); s=s.Transpose(2,4); s=s.Transpose(3,4);}
  else if  ((i==p) && (j==o) && (k==n) && (l==m))
    {t=ten; s=sen.Transpose(1,4); s=s.Transpose(2,3);}
  else
    return PROMOTE(Scal1,Scal2)();
  PROMOTE(Scal1,Scal2) res(ten.at(0,0,0,0)-ten.at(0,0,0,0));
  for (unsigned short int a=0; a<4; ++a)
    for (unsigned short int b=0; b<4; ++b)
      for (unsigned short int c=0; c<4; ++c)
        for (unsigned short int d=0; d<4; ++d) {
          double ga(-1.),gb(-1.),gc(-1.),gd(-1.);
          if (a==0) ga=1.; if (b==0) gb=1.; if (c==0) gc=1.; if (d==0) gd=1.;
          res += ga*gb*gc*gd*t.at(a,b,c,d)*s.at(a,b,c,d);
        }
  return res;
}

template<typename Scal1, typename Scal2>
Lorentz_Ten2<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, int j, int k,
                    const Lorentz_Ten4<Scal2>& sen, int l, int m, int n) {
  Lorentz_Ten4<Scal1> t; Lorentz_Ten4<Scal2> s;
  // T^{mu,nu} = t^{mu,alpha,beta,gamma} s^{nu,alpha',beta',gamma'}
  //                    g_{alpha,alpha'} g_{beta,beta'} g_{gamma,gamma'}
  if      (((i==l) && (j==m) && (k==n)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen;}
  else if (((i==l) && (j==m) && (k==n)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2);}
  else if (((i==l) && (j==m) && (k==n)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3);}
  else if (((i==l) && (j==m) && (k==n)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4);}

  else if (((i==l) && (j==n) && (k==m)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen.Transpose(3,4);}
  else if (((i==l) && (j==n) && (k==m)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2); s=s.Transpose(3,4);}
  else if (((i==l) && (j==n) && (k==m)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3); s=s.Transpose(2,4);}
  else if (((i==l) && (j==n) && (k==m)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4); s=s.Transpose(2,3);}

  else if (((i==m) && (j==l) && (k==n)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen.Transpose(2,3);}
  else if (((i==m) && (j==l) && (k==n)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2); s=s.Transpose(2,3);}
  else if (((i==m) && (j==l) && (k==n)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3); s=s.Transpose(2,3);}
  else if (((i==m) && (j==l) && (k==n)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4); s=s.Transpose(2,4);}

  else if (((i==m) && (j==n) && (k==l)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen.Transpose(2,4); s=s.Transpose(3,4);}
  else if (((i==m) && (j==n) && (k==l)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2); s=s.Transpose(2,4); s=s.Transpose(3,4);}
  else if (((i==m) && (j==n) && (k==l)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3); s=s.Transpose(2,4); s=s.Transpose(2,3);}
  else if (((i==m) && (j==n) && (k==l)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4); s=s.Transpose(2,4); s=s.Transpose(3,4);}

  else if (((i==n) && (j==l) && (k==m)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen.Transpose(2,3); s=s.Transpose(3,4);}
  else if (((i==n) && (j==l) && (k==m)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2); s=s.Transpose(2,3); s=s.Transpose(3,4);}
  else if (((i==n) && (j==l) && (k==m)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3); s=s.Transpose(3,4); s=s.Transpose(2,3);}
  else if (((i==n) && (j==l) && (k==m)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4); s=s.Transpose(2,3); s=s.Transpose(3,4);}

  else if (((i==n) && (j==m) && (k==l)) && ((i==2) && (j==3) && (k==4)))
    {t=ten; s=sen.Transpose(2,4);}
  else if (((i==n) && (j==m) && (k==l)) && ((i==1) && (j==3) && (k==4)))
    {t=ten.Transpose(1,2); s=sen.Transpose(1,2); s=s.Transpose(2,4);}
  else if (((i==n) && (j==m) && (k==l)) && ((i==1) && (j==2) && (k==4)))
    {t=ten.Transpose(1,3); s=sen.Transpose(1,3); s=s.Transpose(3,4);}
  else if (((i==n) && (j==m) && (k==l)) && ((i==1) && (j==2) && (k==3)))
    {t=ten.Transpose(1,4); s=sen.Transpose(1,4); s=s.Transpose(3,4);}
  else {
    THROW(fatal_error,"tensor contraction not specified for these indizes ...");
    return Lorentz_Ten2<PROMOTE(Scal1,Scal2)>();
  }
  Lorentz_Ten2<PROMOTE(Scal1,Scal2)> end;
  for (unsigned short int I=0; I<4; ++I)
    for (unsigned short int J=0; J<4; ++J) {
      PROMOTE(Scal1,Scal2) res(ten.at(0,0,0,0)-ten.at(0,0,0,0));
      for (unsigned short int a=0; a<4; ++a)
        for (unsigned short int b=0; b<4; ++b)
          for (unsigned short int c=0; c<4; ++c) {
            double ga(-1.),gb(-1.),gc(-1.);
            if (a==0) ga=1.; if (b==0) gb=1.; if (c==0) gc=1.;
            res += ga*gb*gc*t.at(I,a,b,c)*s.at(J,a,b,c);
          }
      end.at(I,J) = res;
    }
  return end;
}

template<typename Scal1, typename Scal2>
Lorentz_Ten4<PROMOTE(Scal1,Scal2)>
ATOOLS::Contraction(const Lorentz_Ten4<Scal1>& ten, int i, int j,
                    const Lorentz_Ten4<Scal2>& sen, int k, int l) {
  // T^{mu,nu,rho,sigma} = t^{mu,nu,alpha,beta} s^{rho,sigma,alpha',beta'}
  //                                        g_{alpha,alpha'} g_{beta,beta'}
  if      ((i!=3) && (j!=4) && (k!=3) && (l!=4)) {
    THROW(fatal_error,"tensor contraction not specified for these indizes ...");
    return Lorentz_Ten4<PROMOTE(Scal1,Scal2)>();
  }
  Lorentz_Ten4<PROMOTE(Scal1,Scal2)> end;
  for (unsigned short int I=0; I<4; ++I)
    for (unsigned short int J=0; J<4; ++J)
      for (unsigned short int K=0; K<4; ++K)
        for (unsigned short int L=0; L<4; ++L) {
          PROMOTE(Scal1,Scal2) res(ten.at(0,0,0,0)-ten.at(0,0,0,0));
          for (unsigned short int a=0; a<4; ++a)
            for (unsigned short int b=0; b<4; ++b) {
              double ga(-1.),gb(-1.);
              if (a==0) ga=1.; if (b==0) gb=1.;
              res += ga*gb*ten.at(I,J,a,b)*sen.at(K,L,a,b);
            }
          end.at(I,J) = res;
        }
  return end;
}
#endif
