#ifndef gggg_H
#define gggg_H

double ggggosmec(double q,double s,double t,double u,double s1,double s2,
		double z1,double z2,double a4,double b4) {
return (2*q*s1*s2*(s1*s2 - t*u)*z1*z2)/(t*u*(-s1 - s2 + t + u)) + 
  (pow(q,2)*(4*pow(s1,3)*s2*
        ((2 - 3*b4 + pow(b4,2))*pow(t,2) + t*u + 
          b4*(1 + b4)*pow(u,2)) + 
       pow(s1,2)*(t*u*(-((-1 + pow(b4,2))*pow(t,2)) + 
             2*(1 - b4 + pow(b4,2))*t*u - (-2 + b4)*b4*pow(u,2)) + 
          2*pow(s2,2)*((5 + 4*a4 - 4*b4)*pow(t,2) + 
             (5 + 2*a4 + 2*b4 - 4*a4*b4)*t*u + 
             (5 - 4*a4 + 4*b4)*pow(u,2)) + 
          2*s2*((-1 + a4)*pow(t,3) + (-11 - 2*a4 + 4*b4)*pow(t,2)*u + 
             (-9 + 2*a4 - 4*b4)*t*pow(u,2) - a4*pow(u,3))) + 
       t*u*(-(pow(a4,2)*pow(s2,2)*pow(t - u,2)) + 
          2*a4*s2*(t - u)*(s2*t - 3*pow(t,2) - 4*t*u - 3*pow(u,2)) + 
          u*(pow(s2,2)*(2*t + u) + 
             8*t*(pow(t,2) + t*u + pow(u,2)) - 
             2*s2*(pow(t,2) + 2*t*u + 3*pow(u,2)))) + 
       2*s1*(2*pow(a4,2)*pow(s2,3)*(pow(t,2) + pow(u,2)) - 
          u*(-2*pow(s2,3)*(t + 2*u) + 
             pow(t,2)*(3*pow(t,2) + 2*t*u + pow(u,2)) - 
             2*s2*t*(pow(t,2) + 4*t*u + pow(u,2)) + 
             pow(s2,2)*(9*pow(t,2) + 11*t*u + pow(u,2))) + 
          a4*s2*(4*s2*t*u*(-t + u) + 
             2*pow(s2,2)*(pow(t,2) - 3*pow(u,2)) - 
             b4*pow(pow(t,2) + pow(u,2),2) + 
             u*(pow(t,3) + pow(t,2)*u - t*pow(u,2) + pow(u,3))) + 
          b4*(t*u*(3*pow(t,3) + pow(t,2)*u - t*pow(u,2) - 
                3*pow(u,3)) + 
             pow(s2,2)*(-pow(t,3) + 2*pow(t,2)*u - 2*t*pow(u,2) + 
                pow(u,3)) + s2*t*
              (pow(t,3) - pow(t,2)*u + t*pow(u,2) + pow(u,3)))))*
     pow(z1,2)*pow(z2,2))/
   (pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) - 
  (4*pow(q,3)*(-2*pow(s1,2)*pow(t,3) - pow(s1,2)*pow(t,2)*u + 
       3*s1*s2*pow(t,2)*u + s1*pow(t,3)*u + 3*s1*s2*t*pow(u,2) - 
       pow(s2,2)*t*pow(u,2) - s1*pow(t,2)*pow(u,2) - 
       s2*pow(t,2)*pow(u,2) - 2*pow(s2,2)*pow(u,3) + 
       s2*t*pow(u,3) + pow(b4,2)*s1*
        (-2*pow(t,4) + 2*t*pow(u,3) - 
          2*pow(s2,2)*(pow(t,2) + pow(u,2)) + 
          s2*(2*pow(t,3) + 3*pow(t,2)*u + 5*t*pow(u,2) + 
             2*pow(u,3)) - s1*
           (t*(t - u)*u + s2*(6*pow(t,2) + 4*t*u + 6*pow(u,2)))) + 
       b4*(pow(s1,2)*(pow(t,2)*(2*t + u) + 
             2*s2*(pow(t,2) - t*u - 3*pow(u,2))) + 
          u*(-8*pow(t,2)*(pow(t,2) + t*u + pow(u,2)) + 
             pow(s2,2)*(-pow(t,2) + 2*t*u + 2*pow(u,2)) + 
             s2*t*(3*pow(t,2) + 4*t*u + 5*pow(u,2))) + 
          2*s1*(-(pow(s2,2)*(pow(t,2) + 3*t*u + 5*pow(u,2))) + 
             s2*(3*pow(t,3) + pow(t,2)*u + 4*t*pow(u,2) + 
                pow(u,3)) + 
             t*(pow(t,3) + pow(t,2)*u + 2*t*pow(u,2) + 2*pow(u,3))\
)) + a4*(pow(s1,2)*(-2*(-1 + b4)*pow(t,3) + 
             (2 - 3*b4 + 2*pow(b4,2))*pow(t,2)*u - 
             (1 - b4 + 2*pow(b4,2))*t*pow(u,2) - 2*b4*pow(u,3) + 
             2*s2*((-5 + 8*b4)*pow(t,2) + 3*(-1 + 2*b4)*t*u + 
                (-1 + 8*b4)*pow(u,2))) + 
          u*(pow(s2,2)*u*(t + 2*u) - 
             8*t*u*(pow(t,2) + t*u + pow(u,2)) + 
             2*s2*(2*pow(t,3) + 2*pow(t,2)*u + t*pow(u,2) + 
                pow(u,3))) + 
          b4*(pow(s2,2)*(-2*pow(t,3) + pow(t,2)*u - 
                3*t*pow(u,2) - 2*pow(u,3)) + 
             8*t*u*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                pow(u,3)) + 
             2*s2*(pow(t,4) - 3*pow(t,3)*u - 
                4*pow(t,2)*pow(u,2) - 5*t*pow(u,3) - pow(u,4))) + 
          s1*(t*u*(5*pow(t,2) + 4*t*u + 3*pow(u,2)) + 
             2*pow(s2,2)*((-3 + 8*b4)*pow(t,2) + (-1 + 6*b4)*t*u + 
                (1 + 8*b4)*pow(u,2)) - 
             2*s2*((-1 + 6*b4)*pow(t,3) - 
                (4 - 9*b4 + pow(b4,2))*pow(t,2)*u + 
                (-1 + 9*b4 + pow(b4,2))*t*pow(u,2) + 
                3*(-1 + 2*b4)*pow(u,3)) + 
             2*pow(b4,2)*(pow(t,4) + pow(t,3)*u - t*pow(u,3) - 
                pow(u,4)) - 
             2*b4*(pow(t,4) + 5*pow(t,3)*u + 
                4*pow(t,2)*pow(u,2) + 3*t*pow(u,3) - pow(u,4)))) + 
       pow(a4,2)*s2*(-2*pow(s1,2)*(pow(t,2) + pow(u,2)) + 
          s1*(2*pow(t,3) + (5 - 2*b4)*pow(t,2)*u + 
             (3 + 2*b4)*t*pow(u,2) + 2*pow(u,3) - 
             2*s2*(3*pow(t,2) + 2*t*u + 3*pow(u,2))) - 
          (t - u)*(2*b4*(pow(t,3) + 2*pow(t,2)*u + pow(u,3) + 
                t*u*(s2 + 2*u)) - 
             u*(s2*t + 2*(pow(t,2) + t*u + pow(u,2))))))*pow(z1,3)*
     pow(z2,3))/(pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) + 
  (4*pow(q,4)*(-2*a4*(pow(s1,2)*(2*t + u)*
           ((-1 - b4 + 2*pow(b4,2))*t - 2*b4*u) + 
          s2*u*(t*(t - u) + s2*(t + 2*u)) - 
          b4*(8*t*u*(pow(t,2) + t*u + pow(u,2)) + 
             pow(s2,2)*(2*pow(t,2) + 5*t*u + 2*pow(u,2)) - 
             2*s2*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                pow(u,3))) + 
          s1*(t*(t - u)*u + 2*s2*
              ((-1 + 2*b4 + 4*pow(b4,2))*pow(t,2) + 
                b4*(-1 + 2*b4)*t*u + (1 + 2*b4)*pow(u,2)) + 
             2*b4*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                pow(u,3)) - 
             2*pow(b4,2)*(3*pow(t,3) + 4*pow(t,2)*u + 
                4*t*pow(u,2) + pow(u,3))) + 
          2*pow(b4,2)*(pow(s2,2)*t*(2*t + u) + 
             4*t*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                pow(u,3)) - s2*
              (3*pow(t,3) + 4*pow(t,2)*u + 4*t*pow(u,2) + pow(u,3))\
)) + pow(a4,2)*(pow(s1,2)*
           (4*(-1 + pow(b4,2))*pow(t,2) + 
             (-5 - 4*b4 + 4*pow(b4,2))*t*u + 4*(-2 + b4)*b4*pow(u,2)) \
+ 4*pow(b4,2)*(pow(t,2) + t*u + pow(u,2))*
           (pow(s2,2) - 2*s2*(t + u) + 2*pow(t + u,2)) + 
          u*(-(pow(s2,2)*(t - 4*u)) - 
             4*s2*(pow(t,2) + t*u + pow(u,2)) + 
             8*u*(pow(t,2) + t*u + pow(u,2))) - 
          4*b4*(pow(s2,2)*u*(t + 2*u) + 
             4*u*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                pow(u,3)) - 
             s2*(pow(t,3) + 4*pow(t,2)*u + 4*t*pow(u,2) + 
                3*pow(u,3))) + 
          2*s1*(s2*((2 + 4*pow(b4,2))*pow(t,2) + 
                pow(1 - 2*b4,2)*t*u + 
                2*(3 - 4*b4 + 2*pow(b4,2))*pow(u,2)) - 
             2*(-1 + 2*b4)*(pow(t,2) + t*u + pow(u,2))*
              (-u + b4*(t + u)))) + 
       b4*(pow(s1,2)*t*(4*(-1 + b4)*t - (2 + b4)*u) + 
          2*s2*u*(t*(t - u) + s2*(t + 2*u)) - 
          b4*(pow(s2,2)*u*(5*t + 4*u) + 
             4*s2*t*(pow(t,2) + t*u + pow(u,2)) - 
             8*pow(t,2)*(pow(t,2) + t*u + pow(u,2))) + 
          2*s1*(s2*((-2 + 6*b4)*pow(t,2) + b4*t*u + 
                2*(1 + b4)*pow(u,2)) - 
             t*(u*(-t + u) + 2*b4*(pow(t,2) + t*u + pow(u,2))))))*
   pow(z1,4)*pow(z2,4))/(pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2));
}

double ggggosmet(double q,double s,double t,double u,double s1,double s2,
		double z1,double z2,double a4,double b4) {
  return (q*s2*(-(t*u*(6*pow(t,2) + s2*((-2 + a4)*t - 3*(-1 + a4)*u))) + 
       pow(s1,2)*(s2*(2*t - 3*u) + u*((-3 + b4)*t + 3*(-1 + b4)*u)) + 
       s1*(4*(1 + a4)*pow(s2,2)*t + 
          u*((11 - 3*b4)*pow(t,2) + 2*(5 + 2*b4)*t*u - 
             3*(-1 + b4)*pow(u,2)) + 
          s2*(u*(-((10 + a4)*t) + (-4 + 3*a4)*u) + 
             b4*(-2*pow(t,2) + 7*t*u + pow(u,2)))))*z1*z2)/
   (2.*(-1 + b4)*t*pow(u,2)*(-s1 - s2 + t + u)) + 
  (pow(q,2)*(pow(s1,3)*s2*
        (2*(-5 + a4 + 10*b4 - 8*pow(b4,2) + 2*pow(b4,3))*pow(t,2) + 
          (-5 - 2*a4*(-3 + b4) + b4 - 4*pow(b4,2))*t*u - 
          4*a4*(-1 + b4)*pow(u,2)) + 
       pow(s1,2)*(2*pow(t,2)*u*
           ((-2 + b4 + 2*pow(b4,2) - pow(b4,3))*t + 
             (-4 + 10*b4 - 7*pow(b4,2) + 3*pow(b4,3))*u) + 
          pow(s2,2)*(2*(-11 - 5*a4 + 24*b4 + 8*a4*b4 - 
                14*pow(b4,2) + 2*pow(b4,3))*pow(t,2) + 
             (-13 + 3*a4 + 12*b4 + 2*a4*b4 - 8*pow(b4,2))*t*u + 
             4*(1 - 3*a4 + pow(a4,2))*(-1 + b4)*pow(u,2)) + 
          s2*(-2*(-1 + b4)*(7 + a4 - 8*b4 + 2*pow(b4,2))*pow(t,3) + 
             (32 - 45*b4 + 27*pow(b4,2) + 
                a4*(-19 + 9*b4 - 2*pow(b4,2)))*pow(t,2)*u + 
             (20 - 27*b4 + 25*pow(b4,2) - 
                3*a4*(9 - 7*b4 + 2*pow(b4,2)))*t*pow(u,2) + 
             4*a4*(-1 + b4)*pow(u,3))) + 
       s1*(pow(s2,3)*(4*(-3 + 7*b4 - 3*pow(b4,2) + a4*(-3 + 4*b4))*
              pow(t,2) + (-8 + 11*b4 - 4*pow(b4,2) + 
                a4*(-3 + 4*b4))*t*u + 
             4*pow(-1 + a4,2)*(-1 + b4)*pow(u,2)) + 
          2*pow(t,2)*u*(2*(2 - 3*b4 + pow(b4,2))*pow(t,2) - 
             (-6 + 7*b4 + pow(b4,2))*t*u + 
             (8 - 14*b4 + 5*pow(b4,2))*pow(u,2)) + 
          pow(s2,2)*((16 + a4*(16 - 18*b4) - 42*b4 + 16*pow(b4,2))*
              pow(t,3) + (41 + 2*pow(a4,2) - 72*b4 + 
                35*pow(b4,2) - 2*a4*b4*(4 + 3*b4))*pow(t,2)*u + 
             (26 + 2*pow(a4,2) - 32*b4 + 17*pow(b4,2) + 
                a4*(-19 + 8*b4 - 2*pow(b4,2)))*t*pow(u,2) - 
             4*pow(-1 + a4,2)*(-1 + b4)*pow(u,3)) - 
          s2*t*(2*(2 - 7*b4 + 2*pow(b4,2) + a4*(2 + b4))*pow(t,3) + 
             (25 - 48*b4 + 27*pow(b4,2) + 
                a4*(-17 + 17*b4 - 4*pow(b4,2)))*pow(t,2)*u + 
             2*(17 - 24*b4 + 11*pow(b4,2) + 
                a4*(-12 + 5*b4 + 4*pow(b4,2)))*t*pow(u,2) + 
             (21 - 32*b4 + 17*pow(b4,2) + 
                a4*(-23 + 21*b4 - 4*pow(b4,2)))*pow(u,3))) + 
       t*(b4*(8*t*pow(u,2)*(pow(t,2) + t*u + pow(u,2)) + 
             pow(s2,3)*(-4*pow(t,2) + (3 - 5*a4)*t*u + 
                3*(-1 + a4)*pow(u,2)) + 
             pow(s2,2)*(-4*(-2 + a4)*pow(t,3) + 
                (-9 + a4)*pow(t,2)*u + 
                2*(4 - 6*a4 + 3*pow(a4,2))*t*pow(u,2) + 
                (5 - 3*a4 - 2*pow(a4,2))*pow(u,3)) + 
             2*s2*(2*(-1 + a4)*pow(t,4) + (3 + 4*a4)*pow(t,3)*u + 
                (-3 + 8*a4)*pow(t,2)*pow(u,2) + 
                (-4 + a4)*t*pow(u,3) + 3*(-1 + a4)*pow(u,4))) - 
          u*(-3*pow(s2,3)*u + 2*pow(a4,2)*pow(s2,2)*(t - u)*u + 
             8*t*u*(pow(t,2) + t*u + pow(u,2)) + 
             pow(s2,2)*(-8*pow(t,2) + 9*t*u + 5*pow(u,2)) + 
             2*s2*(4*pow(t,3) - 3*pow(t,2)*u - 4*t*pow(u,2) - 
                3*pow(u,3)) + 
             a4*s2*(pow(s2,2)*(-2*t + 3*u) - 
                s2*(2*pow(t,2) + 9*t*u + 3*pow(u,2)) + 
                2*(2*pow(t,3) + 6*pow(t,2)*u + t*pow(u,2) + 
                   3*pow(u,3))))))*pow(z1,2)*pow(z2,2))/
   ((-1 + b4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) - 
  (2*pow(q,3)*(pow(b4,2)*(16*t*u*((-1 + a4)*t + a4*u)*
           (pow(t,2) + t*u + pow(u,2)) + 
          pow(s2,3)*(4*(-2 + a4)*pow(t,2) + (-7 + 4*a4)*t*u + 
             4*(-1 + a4)*pow(u,2)) + 
          pow(s2,2)*(-4*(-4 + 3*a4)*pow(t,3) + 
             (-3 + 10*a4 - 10*pow(a4,2))*pow(t,2)*u + 
             (7 - 6*a4 + 2*pow(a4,2))*t*pow(u,2) - 
             4*(-1 + pow(a4,2))*pow(u,3)) + 
          2*s2*(4*(-1 + a4)*pow(t,4) + 
             (13 - 12*a4 + 5*pow(a4,2))*pow(t,3)*u + 
             pow(3 - 2*a4,2)*pow(t,2)*pow(u,2) + 
             (6 - 14*a4 + pow(a4,2))*t*pow(u,3) + 
             2*(-1 + a4)*a4*pow(u,4))) + 
       b4*(-16*t*u*((-1 + a4)*t + 2*a4*u)*
           (pow(t,2) + t*u + pow(u,2)) + 
          pow(s2,3)*(4*pow(t,2) - 3*(-3 + a4)*t*u - 
             8*(-1 + a4)*pow(u,2)) + 
          pow(s2,2)*(4*(-2 + a4)*pow(t,3) + 
             (1 - a4 + 6*pow(a4,2))*pow(t,2)*u + 
             3*(-5 + 3*a4)*t*pow(u,2) + 8*(-1 + pow(a4,2))*pow(u,3)\
) - 2*s2*(2*(-1 + a4)*pow(t,4) + 
             (13 - 8*a4 + 4*pow(a4,2))*pow(t,3)*u + 
             (8 - 18*a4 + 9*pow(a4,2))*pow(t,2)*pow(u,2) + 
             (5 - 20*a4 + pow(a4,2))*t*pow(u,3) + 
             4*(-1 + a4)*a4*pow(u,4))) + 
       u*(-(pow(s2,3)*((2 + a4)*t - 4*(-1 + a4)*u)) + 
          16*a4*t*u*(pow(t,2) + t*u + pow(u,2)) + 
          2*s2*((1 + 2*a4)*pow(t,3) + 2*a4*(-3 + 2*a4)*pow(t,2)*u - 
             (1 + 6*a4)*t*pow(u,2) + 2*(-1 + a4)*a4*pow(u,3)) - 
          pow(s2,2)*(3*a4*t*(t + u) - 4*u*(2*t + u) + 
             2*pow(a4,2)*u*(t + 2*u))) + 
       s1*(-2*pow(s2,2)*(2*
              (-1 + 2*b4 + pow(a4,2)*b4 - 4*pow(b4,2) + 
                4*pow(b4,3) - 
                a4*(1 - 3*b4 + 2*pow(b4,2) + 2*pow(b4,3)))*
              pow(t,2) + (pow(a4,2)*(-3 + 4*b4 - 2*pow(b4,2)) + 
                b4*(-3 - 2*b4 + 2*pow(b4,2)) + 
                a4*(1 + 5*b4 - 2*pow(b4,3)))*t*u - 
             2*(-1 + b4)*(1 - b4 + pow(a4,2)*(-3 + 2*b4) + 
                a4*(-1 + 2*b4))*pow(u,2)) + 
          s2*(4*(-1 + 4*b4 + pow(a4,2)*b4 - 7*pow(b4,2) + 
                5*pow(b4,3) + a4*(-1 + 2*b4 - 3*pow(b4,3)))*
              pow(t,3) + (-8 + 25*b4 - 47*pow(b4,2) + 
                10*pow(b4,3) - 
                2*pow(a4,2)*(6 - 11*b4 + 7*pow(b4,2)) + 
                a4*(9 - 11*b4 + 36*pow(b4,2) - 10*pow(b4,3)))*
              pow(t,2)*u + (-2*pow(a4,2)*
                 (7 - 10*b4 + 5*pow(b4,2)) + 
                b4*(13 - 31*b4 + 8*pow(b4,2)) + 
                a4*(9 - 25*b4 + 44*pow(b4,2) - 14*pow(b4,3)))*t*
              pow(u,2) - 4*a4*(-1 + b4)*(1 + a4*(-4 + 3*b4))*pow(u,3)\
) + 2*t*(-4*a4*u*(2*pow(t,2) + 3*t*u + 4*pow(u,2)) + 
             pow(b4,3)*((-1 + a4)*t + a4*u)*
              (2*pow(t,2) - t*u + 5*pow(u,2)) + 
             pow(b4,2)*(-6*(-1 + a4)*pow(t,3) + 
                (13 - 23*a4)*pow(t,2)*u + (17 - 28*a4)*t*pow(u,2) - 
                27*a4*pow(u,3)) + 
             2*b4*(2*(-1 + a4)*pow(t,3) + (-5 + 13*a4)*pow(t,2)*u + 
                (-5 + 17*a4)*t*pow(u,2) + 19*a4*pow(u,3)))) + 
       pow(s1,2)*(s2*(4*(1 - 3*b4 - pow(a4,2)*b4 + 6*pow(b4,2) - 
                4*pow(b4,3) + 
                a4*(1 - 3*b4 + pow(b4,2) + 2*pow(b4,3)))*pow(t,2) + 
             (2 - 3*b4 + 11*pow(b4,2) - 4*pow(b4,3) + 
                2*pow(a4,2)*(3 - 4*b4 + 2*pow(b4,2)) + 
                a4*(-1 - 7*b4 - 4*pow(b4,2) + 4*pow(b4,3)))*t*u + 
             4*a4*(-1 + b4)*(b4 + a4*(-3 + 2*b4))*pow(u,2)) - 
          2*t*(-4*a4*u*(t + 2*u) + 
             pow(b4,2)*(-6*(-1 + a4)*pow(t,2) + (5 - 7*a4)*t*u - 
                19*a4*pow(u,2)) + 
             2*b4*(2*(-1 + a4)*pow(t,2) + (-1 + 4*a4)*t*u + 
                11*a4*pow(u,2)) + 
             pow(b4,3)*(2*(-1 + a4)*pow(t,2) + 5*a4*pow(u,2) + 
                t*(u - a4*u)))))*pow(z1,3)*pow(z2,3))/
   ((-1 + b4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) + 
  (8*pow(q,4)*(-(a4*b4*(pow(s1,2)*t*
             ((4 - 6*b4 + 6*pow(b4,2))*t + (2 - 5*b4 + 3*pow(b4,2))*u) \
+ 2*u*(pow(s2,2)*(t + 2*u) + 4*t*(pow(t,2) + t*u + pow(u,2)) - 
               s2*(5*pow(t,2) + 5*t*u + 2*pow(u,2))) + 
            s1*(2*s2*(2*(1 - 3*b4 + 4*pow(b4,2))*pow(t,2) + 
                  (2 - 7*b4 + 5*pow(b4,2))*t*u + 
                  2*pow(-1 + b4,2)*pow(u,2)) - 
               t*(2*(2 - 7*b4 + 7*pow(b4,2))*pow(t,2) + 
                  (10 - 25*b4 + 17*pow(b4,2))*t*u + 
                  (10 - 21*b4 + 11*pow(b4,2))*pow(u,2))) + 
            pow(b4,2)*(pow(s2,2)*
                (10*pow(t,2) + 7*t*u + 4*pow(u,2)) + 
               8*t*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                  pow(u,3)) - 
               s2*(18*pow(t,3) + 25*pow(t,2)*u + 19*t*pow(u,2) + 
                  4*pow(u,3))) - 
            b4*(pow(s2,2)*(6*pow(t,2) + 9*t*u + 8*pow(u,2)) + 
               8*t*(pow(t,3) + 3*pow(t,2)*u + 3*t*pow(u,2) + 
                  2*pow(u,3)) - 
               s2*(14*pow(t,3) + 33*pow(t,2)*u + 29*t*pow(u,2) + 
                  8*pow(u,3))))) + 
       pow(b4,2)*t*(-2*(s2 - t)*
           (s1*(2*t + u) + s2*(2*t + u) - 
             2*(pow(t,2) + t*u + pow(u,2))) + 
          b4*(2*pow(s1,2)*t + 2*pow(s2,2)*(3*t + u) + 
             4*t*(pow(t,2) + t*u + pow(u,2)) - 
             s2*(10*pow(t,2) + 7*t*u + 4*pow(u,2)) + 
             s1*(-3*t*(2*t + u) + 2*s2*(4*t + u)))) + 
       pow(a4,2)*(4*pow(b4,3)*pow(s1 + s2 - t - u,2)*
           (pow(t,2) + t*u + pow(u,2)) - 
          2*(s1 - u)*u*(s1*(t + 2*u) + s2*(t + 2*u) - 
             2*(pow(t,2) + t*u + pow(u,2))) - 
          pow(b4,2)*(3*pow(s1,2)*
              (2*pow(t,2) + 3*t*u + 4*pow(u,2)) + 
             pow(s2,2)*(2*pow(t,2) + 5*t*u + 8*pow(u,2)) - 
             s2*(6*pow(t,3) + 23*pow(t,2)*u + 29*t*pow(u,2) + 
                20*pow(u,3)) + 
             4*(pow(t,4) + 5*pow(t,3)*u + 8*pow(t,2)*pow(u,2) + 
                7*t*pow(u,3) + 3*pow(u,4)) + 
             s1*(-10*pow(t,3) - 31*pow(t,2)*u - 37*t*pow(u,2) - 
                24*pow(u,3) + 2*s2*(4*pow(t,2) + 7*t*u + 10*pow(u,2))\
)) + b4*(pow(s1,2)*(4*pow(t,2) + 7*t*u + 12*pow(u,2)) + 
             s1*(-4*pow(t,3) - 20*pow(t,2)*u - 27*t*pow(u,2) - 
                24*pow(u,3) + 4*s2*(pow(t,2) + 2*t*u + 4*pow(u,2))) + 
             u*(pow(s2,2)*(t + 4*u) - 
                s2*(8*pow(t,2) + 15*t*u + 16*pow(u,2)) + 
                4*(2*pow(t,3) + 5*pow(t,2)*u + 5*t*pow(u,2) + 
                   3*pow(u,3))))))*pow(z1,4)*pow(z2,4))/
   ((-1 + b4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2));
}

double ggggosmeu(double q,double s,double t,double u,double s1,double s2,
		double z1,double z2,double a4,double b4) {
  return    (2.*(-1 + a4)*pow(t,2)*u*(-s1 - s2 + t + u)) + 
  (pow(q,2)*(2*t*pow(u,2)*
        (pow(s2,2)*((-4 + 10*a4 - 7*pow(a4,2) + 3*pow(a4,3))*t + 
             (-2 + a4 + 2*pow(a4,2) - pow(a4,3))*u) + 
          4*(-1 + a4)*t*(pow(t,2) + t*u + pow(u,2)) + 
          s2*((8 - 14*a4 + 5*pow(a4,2))*pow(t,2) - 
             (-6 + 7*a4 + pow(a4,2))*t*u + 
             2*(2 - 3*a4 + pow(a4,2))*pow(u,2))) + 
       pow(s1,3)*(u*(3*(-1 + a4)*(-1 + b4)*pow(t,2) + 
             (3*a4 + 2*b4 - 5*a4*b4)*t*u - 4*a4*pow(u,2)) + 
          s2*(4*(-1 + a4)*pow(-1 + b4,2)*pow(t,2) + 
             (-8 + 11*a4 - 4*pow(a4,2) - 3*b4 + 4*a4*b4)*t*u - 
             4*(3 - 7*a4 + 3*pow(a4,2) + 3*b4 - 4*a4*b4)*pow(u,2))) + 
       pow(s1,2)*(pow(s2,2)*
           (4*(-1 + a4)*(1 - 3*b4 + pow(b4,2))*pow(t,2) + 
             (-13 - 8*pow(a4,2) + 3*b4 + 2*a4*(6 + b4))*t*u + 
             2*(-11 - 14*pow(a4,2) + 2*pow(a4,3) - 5*b4 + 
                8*a4*(3 + b4))*pow(u,2)) + 
          u*(-((-1 + a4)*(-5 + 3*b4 + 2*pow(b4,2))*pow(t,3)) + 
             (-9 + 9*b4 - 2*pow(b4,2) + 
                2*a4*(4 - 6*b4 + 3*pow(b4,2)))*pow(t,2)*u + 
             (a4*(-9 + b4) + 2*(4 + b4))*t*pow(u,2) - 
             4*a4*(-2 + b4)*pow(u,3)) + 
          s2*(-4*(-1 + a4)*pow(-1 + b4,2)*pow(t,3) + 
             (26 + pow(a4,2)*(17 - 2*b4) + 8*a4*(-4 + b4) - 19*b4 + 
                2*pow(b4,2))*pow(t,2)*u + 
             (41 + pow(a4,2)*(35 - 6*b4) + 2*pow(b4,2) - 
                8*a4*(9 + b4))*t*pow(u,2) + 
             2*(8*pow(a4,2) + 8*(1 + b4) - 3*a4*(7 + 3*b4))*pow(u,3))) \
- s1*(u*(pow(s2,3)*((5 - a4 + 4*pow(a4,2))*t + 
                2*(5 - 10*a4 + 8*pow(a4,2) - 2*pow(a4,3))*u) - 
             pow(s2,2)*((20 - 27*a4 + 25*pow(a4,2))*pow(t,2) + 
                (32 - 45*a4 + 27*pow(a4,2))*t*u + 
                2*(7 - 15*a4 + 10*pow(a4,2) - 2*pow(a4,3))*pow(u,2)\
) + s2*((21 - 32*a4 + 17*pow(a4,2))*pow(t,3) + 
                2*(17 - 24*a4 + 11*pow(a4,2))*pow(t,2)*u + 
                (25 - 48*a4 + 27*pow(a4,2))*t*pow(u,2) + 
                2*(2 - 7*a4 + 2*pow(a4,2))*pow(u,3)) + 
             2*(3*(-1 + a4)*pow(t,4) + 4*(-1 + a4)*pow(t,3)*u + 
                3*(-1 + a4)*pow(t,2)*pow(u,2) + 
                (4 - 3*a4)*t*pow(u,3) + 2*a4*pow(u,4))) + 
          b4*(2*pow(s2,3)*((-1 + a4)*t - u)*(2*t + u) + 
             pow(s2,2)*(-4*(-1 + a4)*pow(t,3) + 
                3*(9 - 7*a4 + 2*pow(a4,2))*pow(t,2)*u + 
                (19 - 9*a4 + 2*pow(a4,2))*t*pow(u,2) + 
                2*(-1 + a4)*pow(u,3)) + 
             s2*u*((-23 + 21*a4 - 4*pow(a4,2))*pow(t,3) + 
                2*(-12 + 5*a4 + 4*pow(a4,2))*pow(t,2)*u + 
                (-17 + 17*a4 - 4*pow(a4,2))*t*pow(u,2) + 
                2*(2 + a4)*pow(u,3)) - 
             2*u*(3*(-1 + a4)*pow(t,4) + (-1 + a4)*pow(t,3)*u + 
                2*(-3 + 4*a4)*pow(t,2)*pow(u,2) + 
                2*(-1 + 2*a4)*t*pow(u,3) + 2*a4*pow(u,4)))))*
     pow(z1,2)*pow(z2,2))/
   ((-1 + a4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) - 
  (2*pow(q,3)*(pow(s1,3)*(4*pow(-1 + a4,2)*(-1 + b4)*pow(t,2) + 
          (-1 + a4)*(2 + b4 + a4*(-7 + 4*b4))*t*u + 
          4*a4*(1 + a4*(-2 + b4))*pow(u,2)) + 
       pow(s1,2)*(-4*pow(-1 + a4,2)*(-1 + pow(b4,2))*pow(t,3) + 
          (-1 + a4)*(-8 + 3*b4 + 2*pow(b4,2) + 
             a4*(7 - 6*b4 + 2*pow(b4,2)))*pow(t,2)*u + 
          (-3*b4 + pow(a4,2)*(-3 + 10*b4 - 10*pow(b4,2)) + 
             a4*(1 - b4 + 6*pow(b4,2)))*t*pow(u,2) + 
          4*a4*(-2 + a4*(4 - 3*b4) + b4)*pow(u,3) + 
          2*s2*(2*(-1 + a4)*(1 - b4 - 3*pow(b4,2) + 
                a4*(-1 + 2*b4 + 2*pow(b4,2)))*pow(t,2) + 
             (2*pow(a4,3)*(-1 + b4) + b4*(-1 + 3*b4) + 
                a4*(3 - 5*b4 - 4*pow(b4,2)) + 
                2*pow(a4,2)*(1 + pow(b4,2)))*t*u + 
             2*(1 + 2*pow(a4,3)*(-2 + b4) + b4 + 
                2*pow(a4,2)*(2 + b4) - a4*(2 + 3*b4 + pow(b4,2)))*
              pow(u,2))) + 2*u*
        (a4*u*(pow(s2,2)*(-((-2 + 5*a4 + pow(a4,2))*t) + 
                2*(2 - 3*a4 + pow(a4,2))*u) - 
             8*(-1 + a4)*t*(pow(t,2) + t*u + pow(u,2)) + 
             s2*((-10 + 17*a4 - 5*pow(a4,2))*pow(t,2) + 
                (-10 + 13*a4 + pow(a4,2))*t*u - 
                2*(2 - 3*a4 + pow(a4,2))*pow(u,2))) + 
          b4*(8*(-1 + a4)*t*((-1 + a4)*t + a4*u)*
              (pow(t,2) + t*u + pow(u,2)) + 
             pow(s2,2)*((8 - 22*a4 + 19*pow(a4,2) - 5*pow(a4,3))*
                 pow(t,2) + 
                (4 - 8*a4 + 7*pow(a4,2) + pow(a4,3))*t*u - 
                2*a4*(2 - 3*a4 + pow(a4,2))*pow(u,2)) + 
             s2*((-16 + 38*a4 - 27*pow(a4,2) + 5*pow(a4,3))*
                 pow(t,3) + 
                2*(-6 + 17*a4 - 14*pow(a4,2) + 2*pow(a4,3))*
                 pow(t,2)*u + 
                (-8 + 26*a4 - 23*pow(a4,2) + pow(a4,3))*t*
                 pow(u,2) + 2*a4*(2 - 3*a4 + pow(a4,2))*pow(u,3)))) \
+ s1*(2*u*(-pow(t,3) + t*pow(u,2) - 2*s2*u*(2*t + u) + 
             pow(s2,2)*(t + 2*u)) + 
          2*pow(a4,3)*s2*u*((4 - 7*b4)*pow(t,2) - 5*(-1 + b4)*t*u + 
             2*(5 - 3*b4)*pow(u,2) + 2*s2*((-1 + b4)*t + 2*(-2 + b4)*u)) \
+ 2*pow(b4,2)*t*(3*pow(s2,2)*(2*t + u) + 
             2*t*(pow(t,2) + 2*pow(u,2)) - 
             s2*(8*pow(t,2) + 7*t*u + 6*pow(u,2))) + 
          b4*(pow(s2,2)*u*(-t + 4*u) + 
             s2*(4*pow(t,3) + 9*pow(t,2)*u + 9*t*pow(u,2) - 
                4*pow(u,3)) - 
             4*t*(pow(t,3) + 3*pow(t,2)*u + 3*t*pow(u,2) - 
                pow(u,3))) + 
          pow(a4,2)*(4*b4*(-pow(t,4) - 7*pow(t,3)*u - 
                6*pow(t,2)*pow(u,2) - 6*t*pow(u,3) + 
                2*pow(u,4) + s2*t*u*(11*t + 9*u) + 
                pow(s2,2)*(pow(t,2) - t*u + pow(u,2))) + 
             2*pow(b4,2)*t*(2*pow(t,3) + pow(t,2)*u + 
                4*t*pow(u,2) + 5*pow(u,3) + 
                2*pow(s2,2)*(2*t + u) - 
                s2*(6*pow(t,2) + 5*t*u + 7*pow(u,2))) + 
             u*(pow(s2,2)*(11*t + 24*u) - 
                s2*(31*pow(t,2) + 47*t*u + 28*pow(u,2)) + 
                2*(6*pow(t,3) + 9*pow(t,2)*u + 13*t*pow(u,2) - 
                   4*pow(u,3)))) - 
          a4*(u*(10*pow(t,3) + 16*pow(t,2)*u + 26*t*pow(u,2) - 
                4*pow(u,3) + 3*pow(s2,2)*(t + 4*u) - 
                s2*(13*pow(t,2) + 25*t*u + 16*pow(u,2))) + 
             2*pow(b4,2)*(2*pow(s2,2)*
                 (5*pow(t,2) + 2*t*u + pow(u,2)) - 
                s2*(14*pow(t,3) + 10*pow(t,2)*u + 11*t*pow(u,2) + 
                   2*pow(u,3)) + 
                t*(4*pow(t,3) + pow(t,2)*u + 9*t*pow(u,2) + 
                   4*pow(u,3))) + 
             b4*(pow(s2,2)*(4*pow(t,2) + 7*t*u + 12*pow(u,2)) + 
                s2*(4*pow(t,3) + 25*pow(t,2)*u + 11*t*pow(u,2) - 
                   8*pow(u,3)) - 
                4*(2*pow(t,4) + 10*pow(t,3)*u + 
                   9*pow(t,2)*pow(u,2) + 4*t*pow(u,3) - pow(u,4))))\
))*pow(z1,3)*pow(z2,3))/
   ((-1 + a4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2)) + 
  (8*pow(q,4)*(-2*pow(b4,2)*(s2 - t)*t*
        (s1*(2*t + u) + s2*(2*t + u) - 2*(pow(t,2) + t*u + pow(u,2))) + 
       a4*b4*(pow(s1,2)*t*(4*(-1 + b4)*t + (-2 + b4)*u) - 
          2*u*(pow(s2,2)*(t + 2*u) + 
             4*t*(pow(t,2) + t*u + pow(u,2)) - 
             s2*(5*pow(t,2) + 5*t*u + 2*pow(u,2))) + 
          s1*(t*((4 - 16*b4)*pow(t,2) - 5*(-2 + 3*b4)*t*u + 
                2*(5 - 4*b4)*pow(u,2)) + 
             4*s2*((-1 + 4*b4)*pow(t,2) + (-1 + 2*b4)*t*u + 
                (-1 + b4)*pow(u,2))) + 
          b4*(pow(s2,2)*(12*pow(t,2) + 7*t*u + 4*pow(u,2)) + 
             4*t*(3*pow(t,3) + 5*pow(t,2)*u + 5*t*pow(u,2) + 
                2*pow(u,3)) - 
             s2*(24*pow(t,3) + 27*pow(t,2)*u + 20*t*pow(u,2) + 
                4*pow(u,3)))) + 
       pow(a4,3)*(4*pow(b4,2)*pow(s1 + s2 - t - u,2)*
           (pow(t,2) + t*u + pow(u,2)) + 
          u*(2*pow(s1,2)*(t + 3*u) + 
             s1*(-4*pow(t,2) - 7*t*u - 10*pow(u,2) + 
                2*s2*(t + 4*u)) + 
             u*(2*pow(s2,2) - 3*s2*(t + 2*u) + 
                4*(pow(t,2) + t*u + pow(u,2)))) - 
          b4*(pow(s1,2)*(4*pow(t,2) + 7*t*u + 10*pow(u,2)) + 
             s1*(-4*pow(t,3) - 19*pow(t,2)*u - 25*t*pow(u,2) - 
                18*pow(u,3) + 2*s2*(2*pow(t,2) + 5*t*u + 8*pow(u,2))\
) + u*(3*pow(s2,2)*(t + 2*u) - 
                s2*(11*pow(t,2) + 17*t*u + 14*pow(u,2)) + 
                8*(pow(t,3) + 2*pow(t,2)*u + 2*t*pow(u,2) + 
                   pow(u,3))))) - 
       pow(a4,2)*(2*(s1 - u)*u*
           (s1*(t + 2*u) + s2*(t + 2*u) - 
             2*(pow(t,2) + t*u + pow(u,2))) + 
          pow(b4,2)*(3*pow(s2,2)*
              (4*pow(t,2) + 3*t*u + 2*pow(u,2)) + 
             pow(s1,2)*(8*pow(t,2) + 5*t*u + 2*pow(u,2)) - 
             s2*(24*pow(t,3) + 37*pow(t,2)*u + 31*t*pow(u,2) + 
                10*pow(u,3)) + 
             4*(3*pow(t,4) + 7*pow(t,3)*u + 8*pow(t,2)*pow(u,2) + 
                5*t*pow(u,3) + pow(u,4)) + 
             s1*(-20*pow(t,3) - 29*pow(t,2)*u - 23*t*pow(u,2) - 
                6*pow(u,3) + 2*s2*(10*pow(t,2) + 7*t*u + 4*pow(u,2)))\
) - b4*(pow(s1,2)*(8*pow(t,2) + 9*t*u + 6*pow(u,2)) + 
             s1*(-8*pow(t,3) - 29*pow(t,2)*u - 33*t*pow(u,2) - 
                14*pow(u,3) + 2*s2*(4*pow(t,2) + 7*t*u + 6*pow(u,2))) \
+ u*(pow(s2,2)*(5*t + 6*u) - 
                s2*(21*pow(t,2) + 25*t*u + 14*pow(u,2)) + 
                8*(2*pow(t,3) + 3*pow(t,2)*u + 3*t*pow(u,2) + 
                   pow(u,3))))))*pow(z1,4)*pow(z2,4))/
   ((-1 + a4)*pow(t,2)*pow(s1 + s2 - t - u,2)*pow(u,2));
}

#endif
