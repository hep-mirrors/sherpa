/*=================================================================
   g(1) g(2) -> gamma(3) gamma(4) two-loop amplitudes with a quark 
   in the loop, and with our conventional overall spinor phases, 
   after removing the singular terms a la Catani, and also the 
   (11*N-2*Nf)/6 * ln(s/mu_R^2) terms.
 =================================================================*/

Complex gggamgam2l_pppp(){ return N_c/2. + 3./2./N_c ; }

Complex F_mppp_1(double x, double y){
  double X = log(-x);
  double Y = log(-y);
  return N_c/8. * ( (2. + 4.*x/y/y - 5.*x*x/y/y) 
                     * ( (X + I*PI)*(X + I*PI) + PISQ )
                  - (1. - x*y) * ((X-Y)*(X-Y) + PISQ)
                  + 2. * (9./y - 10.*x) * (X + I*PI) )
   - 1./8./N_c * ( (x*x + 1.)/y/y * ( (X + I*PI)*(X + I*PI) + PISQ )
                 + 1./2. * (x*x + y*y) * ((X-Y)*(X-Y) + PISQ)
	         - 4. * (1./y - x) * (X + I*PI) ) ;
}

Complex gggamgam2l_mppp(){
  double x = sij(2,3)/sij(1,2); 
  double y = sij(1,3)/sij(1,2); 
  return F_mppp_1(x,y) + F_mppp_1(y,x) ;
}

Complex F_ppmp_1(double x, double y){
  double X = log(-x);
  double Y = log(-y);
  return N_c/8. * ( (2. + 6.*x/y/y - 3.*x*x/y/y) 
                     * ( (X + I*PI)*(X + I*PI) + PISQ )
                  - (x-y)*(x-y) * ((X-Y)*(X-Y) + PISQ)
                  + 2. * (9./y - 8.*x) * (X + I*PI) )
   - 1./8./N_c * ( (x*x + 1.)/y/y * ( (X + I*PI)*(X + I*PI) + PISQ )
	         + 1./2. * (x*x + y*y) * ((X-Y)*(X-Y) + PISQ)
	         - 4. * (1./y - x) * (X + I*PI) ) ;
}

Complex gggamgam2l_ppmp(){
  double x = sij(2,3)/sij(1,2); 
  double y = sij(1,3)/sij(1,2); 
  return F_ppmp_1(x,y) + F_ppmp_1(y,x) ;
}

Complex F_mmpp_1(double x, double y){
  double X = log(-x);
  double Y = log(-y);
  double li2mx = li2(-x);
  double li3mx = li3(-x);
  double li4mx = li4(-x);
  double li3my = li3(-y);
  double li4my = li4(-y);
  double li4mxy = li4(-x/y);
  double x2 = x*x;
  double y2 = y*y;
  double X2 = X*X;
  double Y2 = Y*Y;
  double Y3 = Y2*Y;
  double Y4 = Y2*Y2;
  double XmY2 = (X-Y)*(X-Y);
  double XmY3 = XmY2*(X-Y);
  double XmY4 = XmY2*XmY2;
  return     
    N_c * ( - (x2+y2) * ( 4. * li4mx + (Y - 3.*X) * li3mx 
                        + X2 * ( li2mx - PISQ/12. ) 
                        + 1./48. * (X+Y)*(X+Y)*(X+Y)*(X+Y)
                        - 109./720. * PI4 )
          - 1./2. * x2 * ( li3(-x/y) - (X-Y) * li2(-x/y) - ZETA3
                        - 1./2. * XmY3 - Y * (XmY2 + PISQ) )
          - 2. * x*y * ( li3(-x/y) - (X-Y) * li2(-x/y) - ZETA3
                      + 1./2. * Y * (XmY2 + PISQ) )
          + 1./8. * ( y2 + 12.*x*y - 27.*x2 - 8./y + 9./y2 ) * X2
          - 1./8. * ( 38.*x*y - 13. ) * X * Y + PISQ/48. * ( 114.*x*y - 43. ) 
          - 9./4. * (1./y + 2.*x) * X + 1./4. )
+ 1./N_c * ( 2. * x2 * ( li4mx + li4my - li4mxy
                    - X * ( li3mx + li3my )
                    + PISQ/6. * ( li2mx + X*Y - 1./4. * XmY2 )
                    + 1./16. * XmY4 - 1./2. * XmY2 * Y2 
                    - 2./3. * (X-Y) * Y3 - 1./4. * Y4 - 49./720. * PI4 )
         + x * ( 2. * li3mx - 2. * X * li2mx
               - li3(-x/y) + (X-Y) * li2(-x/y) - 3. * ZETA3
               - 2./3. * X * (X2 + PISQ) + 5./12. * (X-Y) * (XmY2 + PISQ)
               + X2 * (X-Y) )
         + 1./4. * ( 2. - 3.*x2 + x2*x2/y2 ) * X2 
         - 1./4. * ( 2.*x*y + 3. ) * X * Y
         + 1./24. * ( 6.*x*y + 7. ) * PISQ
         - ( 1./2./y + x ) * X + 1./4. )
+ I * PI * ( 
    N_c * ( - (x2+y2) * ( li3(-x/y) - (X-Y) * li2(-x/y) - ZETA3
                        + 1./2. * (Y - 3./4.) * (XmY2 + PISQ) )
          + 1./4. * ( 14.*(x-y) - 8./y + 9./y2 ) * X
          - 9./4. * ( 1./y - 1. ) )
   + 1./N_c * ( 2. * x2 * ( li3(-y/x) - ZETA3 )
           - 2. * x * ( li2(-y/x) - 1./4. * XmY2 - PISQ/12 ) 
           + 1./2. * (2.*x2-1.)/y2 * X     
	     - 1./2./y + 1./2. ) ) ;
}

Complex gggamgam2l_mmpp(){
  double x = sij(2,3)/sij(1,2); 
  double y = sij(1,3)/sij(1,2); 
  return F_mmpp_1(x,y) + F_mmpp_1(y,x) ;
}

Complex F_mpmp(double x, double y){
  double X = log(-x);
  double Y = log(-y);
  double li2mx = li2(-x);
  double li3mx = li3(-x);
  double li4mx = li4(-x);
  double li3my = li3(-y);
  double li4my = li4(-y);
  double li4mxy = li4(-x/y);
  double x2 = x*x;
  double y2 = y*y;
  double X2 = X*X;
  double Y2 = Y*Y;
  double X3 = X2*X;
  double Y3 = Y2*Y;
  double X4 = X2*X2;
  double Y4 = Y2*Y2;
  double XmY2 = (X-Y)*(X-Y);
  return     
    N_c * (  2. * (x2+1.)/y2 * ( li4mxy - li4my
                  + 1./2. * (X - 2*Y) * ( li3mx - ZETA3 )
                  - 1./48. * X4 - 1./6. * X * Y3 + 1./24. * Y4
                  + PISQ/24. * (7.*X2 + 2.*Y2) + 7./360. * PI4 )
         + 4. * x*(x-3.)/y2 * ( li4mx + li4mxy - li4my
                  - Y * ( li3mx - ZETA3 ) 
                  + PISQ/6. * ( li2mx + 1./2. * Y2 )
                  - 1./6. * X * Y3 + 1./24. * Y4 - 7./360. * PI4 )
         + 2./3. * ( 8. - x + 30.*x/y ) 
               * ( li3mx - li3my 
                 - (X + Y) * li2mx - 1./2. * X * Y2 )
         - 1./6. * ( 47. + 154.*x/y - 4.*x2/y2 )
               * ( li3mx - X * li2mx - ZETA3 )
         + 1./12. * ( 3. - 2./y - 12.*x/y2 ) * X * ( X2 + 3.*PISQ )
         - 1./3. * y * X * Y2 + PISQ/3. * ( y - x + 2.*x/y ) * Y
         + 5./9. * PISQ * x * ( 1. - 2./y - 10.*x/y2 ) * X
         + 2. * ( 1. + 2./y ) * ZETA3
         + 1./24. * ( y2 - 24.*y + 44. - 8.*x2*x/y ) * (XmY2+PISQ)
         - 1./24. * ( 15. - 14.*x/y - 48.*x/y2 ) * X2
         + 1./24. * ( 8.*x/y + 60. - 24.*y/x + 27.*y2/x2 ) * Y2
         + 4./9. * PISQ * x/y
         + 1./12. * ( 2.*x2 - 54.*x - 27.*y2 ) * ( 1./y * X + 1./x * Y ) ) 
 + 1./N_c * ( 2. * (x2+1.)/y2 * ( 
                li4mxy - li4my 
              - Y * ( li3mx - ZETA3 )
              + 1./2. * X * ( li3mx - ZETA3 )
              + 1./24. * X4 - 1./6. * X * Y3 + 1./24. * Y4 
              + PISQ/12. * Y2 + 7./360. * PI4 )
         + 2. * (x-1.)/y * ( li4mx 
              - 1./2. * X * ( li3mx - ZETA3 )
              + PISQ/6. * ( li2mx - 1./2. * X2 )
              - 1./48. * X4 - 7./180. * PI4 )
         - (2.*x/y - 1.) * ( li3mx - X * li2mx 
                         + ZETA3 - 1./6. * X3 - PISQ/3. * (X + Y) )
         - 2. * (2.*x/y + 1.) * ( li3my + Y * li2mx - ZETA3 
                             + 1./4. * X * ( 2. * Y2 + PISQ ) - 1./8. * X3 )
         + 1./4. * ( 2.*x2 - y2 ) * (XmY2 + PISQ)
         + 1./4. * ( 3. + 2.*x/y2 ) * X2 + 1./4. * (2.-y2)/x2 * Y2 - PISQ/6.
         - 1./2. * ( 2.*x + y2 ) * ( 1./y * X + 1./x * Y ) + 1./2. )
 + I * PI * ( 
      N_c * ( ( 18.*x/y2 + 4.*x/y - 1. ) * ( li3mx - ZETA3 )
            - 1./6. * (x2+1.)/y2 * X * ( 2. * X2 - PISQ )
            - 1./6. * ( 13. - 8.*x + 78.*x/y + 4./y2 ) * li2mx
            - 1./3. * x * ( 3. + 2.*x/y - 10.*x/y2 ) * X2
            - 1./3. * y * Y * ( Y + 2.*X )
            - PISQ/36. * ( 7. + 4.*x + 18.*x/y - 4./y2 )
            - 1./12. * ( 15. - 14.*x/y - 48.*x/y2 ) * X
            + 1./12. * ( 8.*x/y + 60. - 24.*y/x + 27.*y2/x2 ) * Y
            - 1./12. * ( 2.*x/y - 54./y - 27.*y/x ) )
     + 1./N_c * ( - 2./y2 * ( li3mx - ZETA3 ) + 1./6. * (x2+1.)/y2 * X3 
               + (2./y - 1.) * li2mx + 3./4. * (x-1.)/y * X2
               + ( x/y2 + 3./2. ) * X + 1./2. * (2.-y2)/x2 * Y
               + 1./y + 1./2. * y/x ) ) ;
}

Complex gggamgam2l_mpmp(){
  double x = sij(2,3)/sij(1,2); 
  double y = sij(1,3)/sij(1,2); 
  return F_mpmp(x,y) ;
}

Complex gggamgam2l_mppm(){
  double x = sij(2,3)/sij(1,2); 
  double y = sij(1,3)/sij(1,2); 
  return F_mpmp(y,x) ;
}

/*=====================================================================
   Get g(1,h1) g(2,h2) -> gamma(3,h3) gamma(4,h4) 
   2-loop amplitudes for all helicities, using symmetries
======================================================================*/
 
Complex gggamgam2l(int h1, int h2, int h3, int h4){
  int hsum = h1+h2+h3+h4;
//  more (+)'s than (-)'s
  if (hsum == 4){ return gggamgam2l_pppp(); }
  else if (hsum == 2){
    if ((h1 == -1) || (h2 == -1)){ return gggamgam2l_mppp(); }
    else { return gggamgam2l_ppmp(); } }
  else if (hsum == 0){
    if (((h1 == -1) && (h2 == -1)) || ((h3 == -1) && (h4 == -1))){
      return gggamgam2l_mmpp(); } 
    else if (((h2 == -1) && (h3 == -1)) || ((h4 == -1) && (h1 == -1))){
      return gggamgam2l_mppm(); }
    else { return gggamgam2l_mpmp(); } }
//  more (-)'s than (+)'s
  else if (hsum == -2){
    if ((h1 == 1) || (h2 == 1)){ return gggamgam2l_mppp(); }
    else { return gggamgam2l_ppmp(); } }
  else { return gggamgam2l_pppp(); }
}
