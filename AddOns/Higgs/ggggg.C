//=== FIVE-POINT AMPLITUDES ===

Complex ggggamgam_ppppp(int i1, int i2, int i3, int i4, int i5) {
  return spb(i4,i5)*spb(i4,i5)/spa(i1,i2)/spa(i1,i3)/spa(i2,i3) ; }

//

Complex ggggamgam_mpppp(int i1, int i2, int i3, int i4, int i5) {
  return - ( spb(i2,i3)*spb(i2,i3)*spb(i2,i3)/spb(i1,i2)/spa(i4,i5)/spa(i4,i5)/spb(i1,i3)
  + spa(i1,i5)*spa(i1,i5)*spa(i1,i5)*spb(i3,i5)/spa(i1,i2)/spa(i3,i5)/spa(i2,i5)/spa(i4,i5)/spa(i4,i5)
  + spa(i1,i4)*spa(i1,i4)*spa(i1,i4)*spb(i3,i4)/spa(i1,i2)/spa(i3,i4)/spa(i2,i4)/spa(i4,i5)/spa(i4,i5)
  + spa(i1,i2)*spa(i1,i4)*spa(i1,i5)*spb(i4,i5)
     /spa(i2,i3)/spa(i1,i3)/spa(i2,i4)/spa(i2,i5)/spa(i4,i5) ) ; }

Complex ggggamgam_pppmp(int i1, int i2, int i3, int i4, int i5) {
  return - spa(i4,i5)/spa(i1,i3)/spa(i1,i5)/spa(i3,i5) * (
  spa(i1,i4)*spa(i3,i4)*spb(i4,i5)/spa(i1,i2)/spa(i2,i3)
  - spb(i2,i5)*spa(i4,i5)/spa(i2,i5) ) ; }

//

Complex ggggamgam_pppmm_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_pppmm_Ls1 = spb(i1,i2)*spb(i1,i2)/spa(i1,i2)/spa(i1,i3)/spa(i2,i3)/pow(sij(i3,i5),2)
          * ( 2. * spa(i1,i4)*spa(i1,i5)*spa(i2,i4)*spa(i2,i5) + spa(i1,i2)*spa(i1,i2) * spa(i4,i5)*spa(i4,i5) ) ;
  Complex c_pppmm_L0 = 1./3. * spa(i4,i5)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1) /sij(i2,i5)
           * ( spa(i4,i5) * ( spa(i1,i3)*spb(i2,i3)*spa(i3,i4)*spa(i2,i5)/spa(i3,i5)/spa(i1,i4)
                          + spa(i2,i3)*spb(i1,i3)*spa(i3,i5)*spa(i1,i4)/spa(i3,i4)/spa(i2,i5)
                          -  sij(i1,i3) - sij(i2,i3) )
             + 3. * (spb(i1,i3)*spa(i3,i4)*spa(i1,i5)-spb(i2,i3)*spa(i3,i5)*spa(i2,i4)) ) ;
  Complex c_pppmm_ln = - 1./3. * spa(i4,i5)*spa(i4,i5)*spa(i4,i5) * (
            1./spa(i1,i4)/spa(i1,i2)/spa(i2,i3)/spa(i3,i5)
          + 1./spa(i2,i5)/spa(i2,i3)/spa(i1,i3)/spa(i1,i4)
          + spa(i4,i5)*spa(i2,i3)/spa(i1,i2)/spa(i1,i3)/spa(i2,i4)/spa(i2,i5)/spa(i3,i4)/spa(i3,i5) ) ;
  return c_pppmm_Ls1/2. * Ls1(sij(i1,i4),sij(i3,i5),sij(i2,i4),sij(i3,i5))
       + c_pppmm_L0/2. * L0(sij(i1,i4),sij(i2,i5))
       + c_pppmm_ln/2. * Clog1(sij(i1,i4)) ; }

Complex ggggamgam_pppmm(int i1, int i2, int i3, int i4, int i5) {
  return spa(i4,i5)*spa(i4,i5)/spa(i1,i2)/spa(i2,i3)/spa(i3,i1)
  + ggggamgam_pppmm_A(i1,i2,i3,i4,i5) + ggggamgam_pppmm_A(i2,i3,i1,i4,i5)
  + ggggamgam_pppmm_A(i3,i1,i2,i4,i5) + ggggamgam_pppmm_A(i1,i2,i3,i5,i4)
  + ggggamgam_pppmm_A(i2,i3,i1,i5,i4) + ggggamgam_pppmm_A(i3,i1,i2,i5,i4) ; }

//

Complex ggggamgam_mmppp_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_mmppp_Ls1A = spb(i3,i4)*spb(i3,i4)/spa(i3,i4)/spa(i3,i5)/spa(i4,i5)/pow(sij(i2,i5),2)
           * ( 2. * spa(i1,i3)*spa(i2,i3)*spa(i1,i4)*spa(i2,i4) + spa(i1,i2)*spa(i1,i2) * spa(i3,i4)*spa(i3,i4) ) ;
  Complex c_mmppp_L1A = 2. * spa(i1,i5)*spa(i1,i5)*spa(i1,i5)*spa(i2,i4)*spa(i2,i4)*spb(i4,i5)*spb(i4,i5)
           /spa(i1,i3)/spa(i3,i5)/spa(i4,i5)/spa(i4,i5)/pow(sij(i2,i4),2) ;
  Complex c_mmppp_L0A = - 1./sij(i2,i4) * (
            + 1./3. * spa(i1,i2)*spa(i1,i2)*spa(i1,i2)/spa(i3,i4) * ( spb(i4,i5)/spa(i1,i3)/spa(i2,i5)
            - 5. * spb(i1,i5)/spa(i2,i4)/spa(i3,i5) - spb(i3,i5)/spa(i1,i5)/spa(i2,i4) )
            + 2. * spa(i1,i2)*spa(i1,i4)*spa(i1,i5)*spa(i2,i3)*spb(i4,i5)
                    /spa(i1,i3)/spa(i3,i4)/spa(i3,i5)/spa(i4,i5)
            - spa(i1,i2)*spa(i1,i2) * spb(i2,i5) * ( 2. * spa(i1,i2)*spa(i3,i5) - spa(i1,i5)*spa(i2,i3) )
                    /spa(i1,i3)/spa(i3,i4)/spa(i3,i5)/spa(i4,i5)
            - spa(i1,i2)*spa(i1,i2) * spb(i3,i5) * ( 2. * spa(i2,i5)*spa(i3,i4)+spa(i2,i4)*spa(i3,i5))
                    /spa(i2,i4)/spa(i3,i4)/spa(i3,i5)/spa(i4,i5) ) ;
  Complex c_mmppp_lnA = - 1./3. * spa(i1,i2)*spa(i1,i2)*spa(i1,i2)/spa(i2,i3)/spa(i2,i5)/spa(i4,i5)
             * ( spa(i1,i2)/spa(i1,i3)/spa(i1,i4) + spa(i2,i5)/spa(i1,i5)/spa(i3,i4) ) ;
  return c_mmppp_Ls1A/2. * Ls1(sij(i1,i3),sij(i2,i5),sij(i1,i4),sij(i2,i5))
       + c_mmppp_L1A/2. * L1(sij(i1,i3),sij(i2,i4))
       + c_mmppp_L0A/2. * L0(sij(i1,i3),sij(i2,i4))
       + c_mmppp_lnA/2. * Clog1(sij(i1,i4)) ; }

Complex ggggamgam_mmppp_B(int i1, int i2, int i3, int i4, int i5) {
  Complex c_mmppp_Ls1B = spb(i4,i5)*spb(i4,i5) * (spa(i2,i4)*spa(i3,i5)+spa(i2,i5)*spa(i3,i4))
            /spa(i2,i3)/spa(i3,i4)/spa(i3,i5)/spa(i4,i5)/spa(i4,i5)/sij(i2,i3)/sij(i2,i3)
           * ( 2. * spa(i1,i4)*spa(i2,i4)*spa(i1,i5)*spa(i2,i5) + spa(i1,i2)*spa(i1,i2) * spa(i4,i5)*spa(i4,i5) ) ;
  Complex c_mmppp_L0B = 1./3. * spa(i1,i2)/spa(i4,i5)/spa(i5,i3)/spa(i3,i4) /sij(i5,i2)
           * ( spa(i1,i2) * ( spa(i4,i3)*spb(i5,i3)*spa(i3,i1)*spa(i5,i2)/spa(i3,i2)/spa(i4,i1)
                          + spa(i5,i3)*spb(i4,i3)*spa(i3,i2)*spa(i4,i1)/spa(i3,i1)/spa(i5,i2)
                          -  sij(i4,i3) - sij(i5,i3) )
             + 3. * (spb(i4,i3)*spa(i3,i1)*spa(i4,i2)-spb(i5,i3)*spa(i3,i2)*spa(i5,i1)) ) ;
  Complex c_mmppp_lnB = 1./3. * spa(i1,i2)*spa(i1,i2)*spa(i1,i2)/spa(i2,i4)/spa(i2,i5)/spa(i3,i5)
             * ( spa(i1,i2)/spa(i1,i3)/spa(i1,i4) + spa(i2,i5)/spa(i1,i5)/spa(i3,i4) ) ;
  return c_mmppp_Ls1B/2. * Ls1(sij(i1,i4),sij(i2,i3),sij(i1,i5),sij(i2,i3))
       + c_mmppp_L0B/2. * L0(sij(i1,i4),sij(i2,i5))
       + c_mmppp_lnB/2. * Clog1(sij(i1,i3)) ; }

Complex ggggamgam_mmppp(int i1, int i2, int i3, int i4, int i5) {
  return ( spa(i1,i2)*spa(i1,i2)*spa(i1,i2)/spa(i1,i3)/spa(i2,i3)/spa(i4,i5)/spa(i4,i5)
      + spb(i1,i3)*spb(i3,i4)*spb(i3,i5)*spb(i4,i5)
           /spb(i1,i2)/spb(i2,i3)/spb(i1,i4)/spb(i1,i5)/spa(i4,i5)
      - spa(i2,i4)*spb(i3,i4)*spb(i3,i4)*spb(i3,i4)/spb(i1,i3)/spb(i1,i4)/spb(i2,i4)/spa(i4,i5)/spa(i4,i5)
      - spa(i2,i5)*spb(i3,i5)*spb(i3,i5)*spb(i3,i5)/spb(i1,i3)/spb(i1,i5)/spb(i2,i5)/spa(i4,i5)/spa(i4,i5) )
      + ggggamgam_mmppp_A(i1,i2,i3,i4,i5) + ggggamgam_mmppp_A(i1,i2,i3,i5,i4)
      - ggggamgam_mmppp_A(i2,i1,i3,i4,i5) - ggggamgam_mmppp_A(i2,i1,i3,i5,i4)
      + ggggamgam_mmppp_B(i1,i2,i3,i4,i5) - ggggamgam_mmppp_B(i2,i1,i3,i4,i5) ; }

//

inline Complex c_mppmp_L0A2(int i1, int i2, int i3, int i4, int i5) {
  return 1./3. * (
         spa(i1,i4)*spa(i1,i4)*spb(i3,i5)*(spa(i1,i3)*spa(i4,i5)+spa(i1,i4)*spa(i3,i5))
            /spa(i1,i2)/spa(i2,i3)/spa(i3,i5)/spa(i4,i5)
         + spb(i2,i5)*spa(i1,i4)*spa(i1,i4)
           * (spa(i1,i4)*spa(i2,i5)*spa(i3,i5)+spa(i4,i5)*spa(i1,i5)*spa(i2,i3))
           /spa(i2,i3)/spa(i3,i4)/spa(i1,i5)/spa(i2,i5)/spa(i3,i5)
         - 2. * spa(i1,i4)*spa(i1,i4) * (-spa(i2,i3)*spb(i3,i5)*spa(i1,i5)+spb(i2,i5)*spa(i1,i2)*spa(i2,i5))
           /spa(i1,i2)/spa(i2,i3)/spa(i2,i5)/spa(i3,i5)
         + 3. * spa(i1,i5)*spb(i2,i5)*spa(i2,i4)*spa(i1,i4)
           /spa(i2,i3)/spa(i2,i5)/spa(i3,i5) ) ; }
inline Complex c_mppmp_L1A(int i1, int i2, int i3, int i4, int i5) {
  return spa(i1,i2)*spa(i1,i5)*spb(i2,i5)*spa(i3,i4)*spb(i3,i5)*spa(i4,i5)
           /spa(i2,i3)/spa(i2,i5)/spa(i3,i5)
         - (spa(i1,i5)*spb(i5,i3)*spa(i3,i4))*(spa(i1,i5)*spb(i5,i3)*spa(i3,i4))
           *(spa(i1,i3)*spa(i2,i5)+spa(i1,i5)*spa(i2,i3))
           /spa(i1,i2)/spa(i2,i3)/spa(i2,i5)/spa(i3,i5)/spa(i3,i5) ; }
Complex ggggamgam_mppmp_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_mppmp_Ls1A1 = spb(i3,i5)*spb(i3,i5)*(spa(i1,i3)*spa(i2,i5)+spa(i1,i5)*spa(i2,i3))
               /spa(i3,i5)/spa(i3,i5)/spa(i1,i2)/spa(i2,i3)/spa(i2,i5)
             * ( spa(i1,i4)*spa(i1,i4)*spa(i3,i5)*spa(i3,i5) + 2. * spa(i1,i5)*spa(i5,i4)*spa(i1,i3)*spa(i3,i4) ) ;
  Complex c_mppmp_Ls1A2 = spb(i3,i5)*spb(i3,i5)/spa(i2,i3)/spa(i2,i5)/spa(i3,i5)
             * ( spa(i1,i4)*spa(i1,i4)*spa(i3,i5)*spa(i3,i5) + 2. * spa(i1,i5)*spa(i5,i4)*spa(i1,i3)*spa(i3,i4) ) ;
  Complex c_mppmp_L0A1 = 1./3. * spa(i1,i4) * (
               spa(i1,i4)*spa(i1,i3)*spb(i2,i3)/spa(i2,i3)/spa(i3,i5)/spa(i1,i5)
             + spa(i1,i4)*spa(i3,i4)*spb(i3,i5)/spa(i2,i4)/spa(i2,i3)/spa(i3,i5)
             + spa(i1,i4)*spa(i1,i4)*spb(i2,i3)/spa(i3,i4)/spa(i1,i5)/spa(i2,i5)
             + spa(i1,i4)*spa(i1,i4)*spb(i3,i5)/spa(i1,i3)/spa(i2,i5)/spa(i2,i4)
             + spa(i1,i4)*(spa(i3,i5)*spb(i3,i5)+spa(i2,i3)*spb(i2,i3))
                  /spa(i2,i5)/spa(i2,i3)/spa(i3,i5)
             - 3. * (spa(i1,i3)*spb(i2,i3)*spa(i2,i4)+spa(i3,i4)*spb(i3,i5)*spa(i1,i5))
                  /spa(i2,i5)/spa(i2,i3)/spa(i3,i5) ) ;
  Complex c_mppmp_ln12 = 1./3. * spa(i1,i4)*spa(i1,i4)*spa(i1,i4)/spa(i2,i3)/spa(i4,i5)
            * ( spa(i1,i4)/spa(i1,i2)/spa(i3,i4)/spa(i1,i5) - 1./spa(i1,i3)/spa(i2,i5) ) ;
  Complex c_mppmp_ln24 = - 1./3. * spa(i1,i4)*spa(i1,i4)*spa(i1,i4)*spa(i3,i5)/spa(i2,i3)/spa(i1,i5)/spa(i2,i5)
            * ( spa(i1,i2)/spa(i1,i3)/spa(i2,i4)/spa(i3,i5)
              - spa(i1,i4)/spa(i3,i4)/spa(i1,i3)/spa(i4,i5) ) ;
  return c_mppmp_Ls1A1/2. * Ls1(sij(i3,i4),sij(i1,i2),sij(i4,i5),sij(i1,i2))/pow(sij(i1,i2),2)
       + c_mppmp_Ls1A2/2. * Ls1(sij(i1,i3),sij(i2,i4),sij(i1,i5),sij(i2,i4))/pow(sij(i2,i4),2)
       + c_mppmp_L0A1/2. * L0(sij(i1,i5),sij(i2,i4))/sij(i2,i4)
       + c_mppmp_L1A(i1,i2,i3,i4,i5)/2. * L1(sij(i1,i2),sij(i3,i4))/pow(sij(i3,i4),2)
       + c_mppmp_L1A(i1,i2,i5,i4,i3)/2. * L1(sij(i1,i2),sij(i4,i5))/pow(sij(i4,i5),2)
       + c_mppmp_L0A2(i1,i2,i3,i4,i5)/2. * L0(sij(i1,i2),sij(i3,i4))/sij(i3,i4)
       + c_mppmp_L0A2(i1,i2,i5,i4,i3)/2. * L0(sij(i1,i2),sij(i4,i5))/sij(i4,i5)
       + c_mppmp_ln12/2. * Clog1(sij(i1,i2)) + c_mppmp_ln24/2. * Clog1(sij(i2,i4)) ; }
                  
Complex ggggamgam_mppmp_B(int i1, int i2, int i3, int i4, int i5) {
  Complex c_mppmp_Ls1B = spb(i2,i3)*spb(i2,i3)/spa(i2,i3)/spa(i2,i5)/spa(i3,i5)
             * ( spa(i1,i4)*spa(i1,i4)*spa(i2,i3)*spa(i2,i3) + 2. * spa(i1,i2)*spa(i2,i4)*spa(i1,i3)*spa(i3,i4) ) ;
  Complex c_mppmp_ln15 = - 1./3. * spa(i1,i4)*spa(i1,i4)*spa(i1,i4)/spa(i1,i5)/spa(i2,i3)
            * ( spa(i1,i3)/spa(i3,i5)/spa(i1,i2)/spa(i3,i4)
              + spa(i1,i2)/spa(i2,i5)/spa(i1,i3)/spa(i2,i4) ) ;
  Complex c_mppmp_ln45 = 1./3. * spa(i1,i4)*spa(i1,i4)*spa(i1,i4)/spa(i1,i2)/spa(i3,i5)/spa(i1,i3)/spa(i2,i5)
             * ( spa(i1,i4)*spa(i2,i3)/spa(i2,i4)/spa(i3,i4)
              - (spa(i2,i5)*spa(i1,i3)+spa(i1,i2)*spa(i3,i5))/spa(i2,i3)/spa(i4,i5) ) ;
  return c_mppmp_Ls1B/2. * ( - Ls1(sij(i1,i2),sij(i4,i5),sij(i1,i3),sij(i4,i5))/pow(sij(i4,i5),2)
                        + Ls1(sij(i2,i4),sij(i1,i5),sij(i3,i4),sij(i1,i5))/pow(sij(i1,i5),2) )
       + c_mppmp_ln15/2. * Clog1(sij(i1,i5)) + c_mppmp_ln45/2. * Clog1(sij(i4,i5)) ; }

Complex ggggamgam_mppmp(int i1, int i2, int i3, int i4, int i5) {
  return 1./2. * spa(i1,i4)*spa(i1,i4) /spa(i2,i3)/spa(i2,i5)/spa(i3,i5)
         - 1./2. * spa(i1,i4) * spb(i2,i3) * (spa(i1,i2)*spb(i2,i4)-spa(i1,i3)*spb(i3,i4))
                      /spa(i2,i3)/spb(i2,i4)/spb(i3,i4)/spa(i2,i5)/spa(i3,i5)
         - spa(i1,i2) * spa(i1,i5)*spa(i1,i5) * spb(i2,i5)*spb(i2,i5)
                      /spa(i1,i3)/spa(i3,i5)/spb(i2,i4)/spb(i4,i5)/spa(i2,i5)/spa(i2,i5)
         + spa(i1,i3) * spa(i1,i5)*spa(i1,i5) * spb(i3,i5)*spb(i3,i5)
                      /spa(i1,i2)/spa(i2,i5)/spb(i3,i4)/spb(i4,i5)/spa(i3,i5)/spa(i3,i5)
         - spb(i2,i3)*spb(i2,i3) * spb(i2,i5) * spb(i3,i5)
                      /spb(i1,i2)/spb(i1,i3)/spb(i3,i4)/spa(i3,i5)/spb(i4,i5)
         + spb(i2,i3)*spb(i2,i3) * spb(i2,i5) * spb(i3,i5)
                      /spb(i1,i2)/spb(i1,i3)/spb(i2,i4)/spa(i2,i5)/spb(i4,i5)
         + spb(i2,i3) * spb(i2,i5) * spb(i3,i5) * spa(i4,i5)
                      /spb(i1,i2)/spb(i1,i3)/spa(i2,i5)/spa(i3,i5)/spb(i4,i5)
         + ggggamgam_mppmp_A(i1,i2,i3,i4,i5) - ggggamgam_mppmp_A(i1,i3,i2,i4,i5)
         + ggggamgam_mppmp_B(i1,i2,i3,i4,i5) ; }

//---

Complex ggggamgam_mmmmm(int i1, int i2, int i3, int i4, int i5) {
  return - spa(i4,i5)*spa(i4,i5)/spb(i1,i2)/spb(i1,i3)/spb(i2,i3) ; }

//

Complex ggggamgam_pmmmm(int i1, int i2, int i3, int i4, int i5) {
  return ( spa(i2,i3)*spa(i2,i3)*spa(i2,i3)/spa(i1,i2)/spb(i4,i5)/spb(i4,i5)/spa(i1,i3)
  + spb(i1,i5)*spb(i1,i5)*spb(i1,i5)*spa(i3,i5)/spb(i1,i2)/spb(i3,i5)/spb(i2,i5)/spb(i4,i5)/spb(i4,i5)
  + spb(i1,i4)*spb(i1,i4)*spb(i1,i4)*spa(i3,i4)/spb(i1,i2)/spb(i3,i4)/spb(i2,i4)/spb(i4,i5)/spb(i4,i5)
  + spb(i1,i2)*spb(i1,i4)*spb(i1,i5)*spa(i4,i5)
     /spb(i2,i3)/spb(i1,i3)/spb(i2,i4)/spb(i2,i5)/spb(i4,i5) ) ; }

//

Complex ggggamgam_mmmpm(int i1, int i2, int i3, int i4, int i5) {
  return spb(i4,i5)/spb(i1,i3)/spb(i1,i5)/spb(i3,i5) * (
  spb(i1,i4)*spb(i3,i4)*spa(i4,i5)/spb(i1,i2)/spb(i2,i3)
  - spa(i2,i5)*spb(i4,i5)/spb(i2,i5) ) ; }

//

Complex ggggamgam_mmmpp_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_mmmpp_Ls1 = - spa(i1,i2)*spa(i1,i2)/spb(i1,i2)/spb(i1,i3)/spb(i2,i3)/pow(sij(i3,i5),2)
          * ( 2. * spb(i1,i4)*spb(i1,i5)*spb(i2,i4)*spb(i2,i5) + spb(i1,i2)*spb(i1,i2) * spb(i4,i5)*spb(i4,i5) ) ;
  Complex c_mmmpp_L0 = - 1./3. * spb(i4,i5)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1) /sij(i2,i5)
           * ( spb(i4,i5) * ( spb(i1,i3)*spa(i2,i3)*spb(i3,i4)*spb(i2,i5)/spb(i3,i5)/spb(i1,i4)
                          + spb(i2,i3)*spa(i1,i3)*spb(i3,i5)*spb(i1,i4)/spb(i3,i4)/spb(i2,i5)
                          -  sij(i1,i3) - sij(i2,i3) )
             + 3. * (spa(i1,i3)*spb(i3,i4)*spb(i1,i5)-spa(i2,i3)*spb(i3,i5)*spb(i2,i4)) ) ;
  Complex c_mmmpp_ln = 1./3. * spb(i4,i5)*spb(i4,i5)*spb(i4,i5) * (
            1./spb(i1,i4)/spb(i1,i2)/spb(i2,i3)/spb(i3,i5)
          + 1./spb(i2,i5)/spb(i2,i3)/spb(i1,i3)/spb(i1,i4)
          + spb(i4,i5)*spb(i2,i3)/spb(i1,i2)/spb(i1,i3)/spb(i2,i4)/spb(i2,i5)/spb(i3,i4)/spb(i3,i5) ) ;
  return c_mmmpp_Ls1/2. * Ls1(sij(i1,i4),sij(i3,i5),sij(i2,i4),sij(i3,i5))
       + c_mmmpp_L0/2. * L0(sij(i1,i4),sij(i2,i5))
       + c_mmmpp_ln/2. * Clog1(sij(i1,i4)) ; }

Complex ggggamgam_mmmpp(int i1, int i2, int i3, int i4, int i5) {
  return - spb(i4,i5)*spb(i4,i5)/spb(i1,i2)/spb(i2,i3)/spb(i3,i1)
  + ggggamgam_mmmpp_A(i1,i2,i3,i4,i5) + ggggamgam_mmmpp_A(i2,i3,i1,i4,i5)
  + ggggamgam_mmmpp_A(i3,i1,i2,i4,i5) + ggggamgam_mmmpp_A(i1,i2,i3,i5,i4)
  + ggggamgam_mmmpp_A(i2,i3,i1,i5,i4) + ggggamgam_mmmpp_A(i3,i1,i2,i5,i4) ; }

//

Complex ggggamgam_ppmmm_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_ppmmm_Ls1A = - spa(i3,i4)*spa(i3,i4)/spb(i3,i4)/spb(i3,i5)/spb(i4,i5)/pow(sij(i2,i5),2)
           * ( 2. * spb(i1,i3)*spb(i2,i3)*spb(i1,i4)*spb(i2,i4) + spb(i1,i2)*spb(i1,i2) * spb(i3,i4)*spb(i3,i4) ) ;
  Complex c_ppmmm_L1A = - 2. * spb(i1,i5)*spb(i1,i5)*spb(i1,i5)*spb(i2,i4)*spb(i2,i4)*spa(i4,i5)*spa(i4,i5)
           /spb(i1,i3)/spb(i3,i5)/spb(i4,i5)/spb(i4,i5)/pow(sij(i2,i4),2) ;
  Complex c_ppmmm_L0A = 1./sij(i2,i4) * (
            + 1./3. * spb(i1,i2)*spb(i1,i2)*spb(i1,i2)/spb(i3,i4) * ( spa(i4,i5)/spb(i1,i3)/spb(i2,i5)
            - 5. * spa(i1,i5)/spb(i2,i4)/spb(i3,i5) - spa(i3,i5)/spb(i1,i5)/spb(i2,i4) )
            + 2. * spb(i1,i2)*spb(i1,i4)*spb(i1,i5)*spb(i2,i3)*spa(i4,i5)
                    /spb(i1,i3)/spb(i3,i4)/spb(i3,i5)/spb(i4,i5)
            - spb(i1,i2)*spb(i1,i2) * spa(i2,i5) * ( 2. * spb(i1,i2)*spb(i3,i5) - spb(i1,i5)*spb(i2,i3) )
                    /spb(i1,i3)/spb(i3,i4)/spb(i3,i5)/spb(i4,i5)
            - spb(i1,i2)*spb(i1,i2) * spa(i3,i5) * ( 2. * spb(i2,i5)*spb(i3,i4)+spb(i2,i4)*spb(i3,i5))
                    /spb(i2,i4)/spb(i3,i4)/spb(i3,i5)/spb(i4,i5) ) ;
  Complex c_ppmmm_lnA = 1./3. * spb(i1,i2)*spb(i1,i2)*spb(i1,i2)/spb(i2,i3)/spb(i2,i5)/spb(i4,i5)
             * ( spb(i1,i2)/spb(i1,i3)/spb(i1,i4) + spb(i2,i5)/spb(i1,i5)/spb(i3,i4) ) ;
  return c_ppmmm_Ls1A/2. * Ls1(sij(i1,i3),sij(i2,i5),sij(i1,i4),sij(i2,i5))
       + c_ppmmm_L1A/2. * L1(sij(i1,i3),sij(i2,i4))
       + c_ppmmm_L0A/2. * L0(sij(i1,i3),sij(i2,i4))
       + c_ppmmm_lnA/2. * Clog1(sij(i1,i4)) ; }

Complex ggggamgam_ppmmm_B(int i1, int i2, int i3, int i4, int i5) {
  Complex c_ppmmm_Ls1B = - spa(i4,i5)*spa(i4,i5) * (spb(i2,i4)*spb(i3,i5)+spb(i2,i5)*spb(i3,i4))
            /spb(i2,i3)/spb(i3,i4)/spb(i3,i5)/spb(i4,i5)/spb(i4,i5)/sij(i2,i3)/sij(i2,i3)
           * ( 2. * spb(i1,i4)*spb(i2,i4)*spb(i1,i5)*spb(i2,i5) + spb(i1,i2)*spb(i1,i2) * spb(i4,i5)*spb(i4,i5) ) ;
  Complex c_ppmmm_L0B = - 1./3. * spb(i1,i2)/spb(i4,i5)/spb(i5,i3)/spb(i3,i4) /sij(i5,i2)
           * ( spb(i1,i2) * ( spb(i4,i3)*spa(i5,i3)*spb(i3,i1)*spb(i5,i2)/spb(i3,i2)/spb(i4,i1)
                          + spb(i5,i3)*spa(i4,i3)*spb(i3,i2)*spb(i4,i1)/spb(i3,i1)/spb(i5,i2)
                          -  sij(i4,i3) - sij(i5,i3) )
             + 3. * (spa(i4,i3)*spb(i3,i1)*spb(i4,i2)-spa(i5,i3)*spb(i3,i2)*spb(i5,i1)) ) ;
  Complex c_ppmmm_lnB = - 1./3. * spb(i1,i2)*spb(i1,i2)*spb(i1,i2)/spb(i2,i4)/spb(i2,i5)/spb(i3,i5)
             * ( spb(i1,i2)/spb(i1,i3)/spb(i1,i4) + spb(i2,i5)/spb(i1,i5)/spb(i3,i4) ) ;
  return c_ppmmm_Ls1B/2. * Ls1(sij(i1,i4),sij(i2,i3),sij(i1,i5),sij(i2,i3))
       + c_ppmmm_L0B/2. * L0(sij(i1,i4),sij(i2,i5))
       + c_ppmmm_lnB/2. * Clog1(sij(i1,i3)) ; }

Complex ggggamgam_ppmmm(int i1, int i2, int i3, int i4, int i5) {
  return - ( spb(i1,i2)*spb(i1,i2)*spb(i1,i2)/spb(i1,i3)/spb(i2,i3)/spb(i4,i5)/spb(i4,i5)
      + spa(i1,i3)*spa(i3,i4)*spa(i3,i5)*spa(i4,i5)
           /spa(i1,i2)/spa(i2,i3)/spa(i1,i4)/spa(i1,i5)/spb(i4,i5)
      - spb(i2,i4)*spa(i3,i4)*spa(i3,i4)*spa(i3,i4)/spa(i1,i3)/spa(i1,i4)/spa(i2,i4)/spb(i4,i5)/spb(i4,i5)
      - spb(i2,i5)*spa(i3,i5)*spa(i3,i5)*spa(i3,i5)/spa(i1,i3)/spa(i1,i5)/spa(i2,i5)/spb(i4,i5)/spb(i4,i5) )
      + ggggamgam_ppmmm_A(i1,i2,i3,i4,i5) + ggggamgam_ppmmm_A(i1,i2,i3,i5,i4)
      - ggggamgam_ppmmm_A(i2,i1,i3,i4,i5) - ggggamgam_ppmmm_A(i2,i1,i3,i5,i4)
      + ggggamgam_ppmmm_B(i1,i2,i3,i4,i5) - ggggamgam_ppmmm_B(i2,i1,i3,i4,i5) ; }

//

inline Complex c_pmmpm_L0A2(int i1, int i2, int i3, int i4, int i5) {
  return - 1./3. * (
         spb(i1,i4)*spb(i1,i4)*spa(i3,i5)*(spb(i1,i3)*spb(i4,i5)+spb(i1,i4)*spb(i3,i5))
            /spb(i1,i2)/spb(i2,i3)/spb(i3,i5)/spb(i4,i5)
         + spa(i2,i5)*spb(i1,i4)*spb(i1,i4)
           * (spb(i1,i4)*spb(i2,i5)*spb(i3,i5)+spb(i4,i5)*spb(i1,i5)*spb(i2,i3))
           /spb(i2,i3)/spb(i3,i4)/spb(i1,i5)/spb(i2,i5)/spb(i3,i5)
         - 2. * spb(i1,i4)*spb(i1,i4) * (-spb(i2,i3)*spa(i3,i5)*spb(i1,i5)+spa(i2,i5)*spb(i1,i2)*spb(i2,i5))
           /spb(i1,i2)/spb(i2,i3)/spb(i2,i5)/spb(i3,i5)
         + 3. * spb(i1,i5)*spa(i2,i5)*spb(i2,i4)*spb(i1,i4)
           /spb(i2,i3)/spb(i2,i5)/spb(i3,i5) ) ; }
inline Complex c_pmmpm_L1A(int i1, int i2, int i3, int i4, int i5) {
  return - spb(i1,i2)*spb(i1,i5)*spa(i2,i5)*spb(i3,i4)*spa(i3,i5)*spb(i4,i5)
           /spb(i2,i3)/spb(i2,i5)/spb(i3,i5)
         + (spb(i1,i5)*spa(i5,i3)*spb(i3,i4))*(spb(i1,i5)*spa(i5,i3)*spb(i3,i4))
           *(spb(i1,i3)*spb(i2,i5)+spb(i1,i5)*spb(i2,i3))
           /spb(i1,i2)/spb(i2,i3)/spb(i2,i5)/spb(i3,i5)/spb(i3,i5) ; }
Complex ggggamgam_pmmpm_A(int i1, int i2, int i3, int i4, int i5) {
  Complex c_pmmpm_Ls1A1 = - spa(i3,i5)*spa(i3,i5)*(spb(i1,i3)*spb(i2,i5)+spb(i1,i5)*spb(i2,i3))
               /spb(i3,i5)/spb(i3,i5)/spb(i1,i2)/spb(i2,i3)/spb(i2,i5)
             * ( spb(i1,i4)*spb(i1,i4)*spb(i3,i5)*spb(i3,i5) + 2. * spb(i1,i5)*spb(i5,i4)*spb(i1,i3)*spb(i3,i4) ) ;
  Complex c_pmmpm_Ls1A2 = - spa(i3,i5)*spa(i3,i5)/spb(i2,i3)/spb(i2,i5)/spb(i3,i5)
             * ( spb(i1,i4)*spb(i1,i4)*spb(i3,i5)*spb(i3,i5) + 2. * spb(i1,i5)*spb(i5,i4)*spb(i1,i3)*spb(i3,i4) ) ;
  Complex c_pmmpm_L0A1 = - 1./3. * spb(i1,i4) * (
               spb(i1,i4)*spb(i1,i3)*spa(i2,i3)/spb(i2,i3)/spb(i3,i5)/spb(i1,i5)
             + spb(i1,i4)*spb(i3,i4)*spa(i3,i5)/spb(i2,i4)/spb(i2,i3)/spb(i3,i5)
             + spb(i1,i4)*spb(i1,i4)*spa(i2,i3)/spb(i3,i4)/spb(i1,i5)/spb(i2,i5)
             + spb(i1,i4)*spb(i1,i4)*spa(i3,i5)/spb(i1,i3)/spb(i2,i5)/spb(i2,i4)
             + spb(i1,i4)*(spb(i3,i5)*spa(i3,i5)+spb(i2,i3)*spa(i2,i3))
                  /spb(i2,i5)/spb(i2,i3)/spb(i3,i5)
             - 3. * (spb(i1,i3)*spa(i2,i3)*spb(i2,i4)+spb(i3,i4)*spa(i3,i5)*spb(i1,i5))
                  /spb(i2,i5)/spb(i2,i3)/spb(i3,i5) ) ;
  Complex c_pmmpm_ln12 = - 1./3. * spb(i1,i4)*spb(i1,i4)*spb(i1,i4)/spb(i2,i3)/spb(i4,i5)
            * ( spb(i1,i4)/spb(i1,i2)/spb(i3,i4)/spb(i1,i5) - 1./spb(i1,i3)/spb(i2,i5) ) ;
  Complex c_pmmpm_ln24 = 1./3. * spb(i1,i4)*spb(i1,i4)*spb(i1,i4)*spb(i3,i5)/spb(i2,i3)/spb(i1,i5)/spb(i2,i5)
            * ( spb(i1,i2)/spb(i1,i3)/spb(i2,i4)/spb(i3,i5)
              - spb(i1,i4)/spb(i3,i4)/spb(i1,i3)/spb(i4,i5) ) ;
  return c_pmmpm_Ls1A1/2. * Ls1(sij(i3,i4),sij(i1,i2),sij(i4,i5),sij(i1,i2))/pow(sij(i1,i2),2)
       + c_pmmpm_Ls1A2/2. * Ls1(sij(i1,i3),sij(i2,i4),sij(i1,i5),sij(i2,i4))/pow(sij(i2,i4),2)
       + c_pmmpm_L0A1/2. * L0(sij(i1,i5),sij(i2,i4))/sij(i2,i4)
       + c_pmmpm_L1A(i1,i2,i3,i4,i5)/2. * L1(sij(i1,i2),sij(i3,i4))/pow(sij(i3,i4),2)
       + c_pmmpm_L1A(i1,i2,i5,i4,i3)/2. * L1(sij(i1,i2),sij(i4,i5))/pow(sij(i4,i5),2)
       + c_pmmpm_L0A2(i1,i2,i3,i4,i5)/2. * L0(sij(i1,i2),sij(i3,i4))/sij(i3,i4)
       + c_pmmpm_L0A2(i1,i2,i5,i4,i3)/2. * L0(sij(i1,i2),sij(i4,i5))/sij(i4,i5)
       + c_pmmpm_ln12/2. * Clog1(sij(i1,i2)) + c_pmmpm_ln24/2. * Clog1(sij(i2,i4)) ; }

Complex ggggamgam_pmmpm_B(int i1, int i2, int i3, int i4, int i5) {
  Complex c_pmmpm_Ls1B = - spa(i2,i3)*spa(i2,i3)/spb(i2,i3)/spb(i2,i5)/spb(i3,i5)
             * ( spb(i1,i4)*spb(i1,i4)*spb(i2,i3)*spb(i2,i3) + 2. * spb(i1,i2)*spb(i2,i4)*spb(i1,i3)*spb(i3,i4) ) ;
  Complex c_pmmpm_ln15 = 1./3. * spb(i1,i4)*spb(i1,i4)*spb(i1,i4)/spb(i1,i5)/spb(i2,i3)
            * ( spb(i1,i3)/spb(i3,i5)/spb(i1,i2)/spb(i3,i4)
              + spb(i1,i2)/spb(i2,i5)/spb(i1,i3)/spb(i2,i4) ) ;
  Complex c_pmmpm_ln45 = - 1./3. * spb(i1,i4)*spb(i1,i4)*spb(i1,i4)/spb(i1,i2)/spb(i3,i5)/spb(i1,i3)/spb(i2,i5)
             * ( spb(i1,i4)*spb(i2,i3)/spb(i2,i4)/spb(i3,i4)
              - (spb(i2,i5)*spb(i1,i3)+spb(i1,i2)*spb(i3,i5))/spb(i2,i3)/spb(i4,i5) ) ;
  return c_pmmpm_Ls1B/2. * ( - Ls1(sij(i1,i2),sij(i4,i5),sij(i1,i3),sij(i4,i5))/pow(sij(i4,i5),2)
                        + Ls1(sij(i2,i4),sij(i1,i5),sij(i3,i4),sij(i1,i5))/pow(sij(i1,i5),2) )
       + c_pmmpm_ln15/2. * Clog1(sij(i1,i5)) + c_pmmpm_ln45/2. * Clog1(sij(i4,i5)) ; }

Complex ggggamgam_pmmpm(int i1, int i2, int i3, int i4, int i5) {
  return - 1./2. * spb(i1,i4)*spb(i1,i4) /spb(i2,i3)/spb(i2,i5)/spb(i3,i5)
         + 1./2. * spb(i1,i4) * spa(i2,i3) * (spb(i1,i2)*spa(i2,i4)-spb(i1,i3)*spa(i3,i4))
                      /spb(i2,i3)/spa(i2,i4)/spa(i3,i4)/spb(i2,i5)/spb(i3,i5)
         + spb(i1,i2) * spb(i1,i5)*spb(i1,i5) * spa(i2,i5)*spa(i2,i5)
                      /spb(i1,i3)/spb(i3,i5)/spa(i2,i4)/spa(i4,i5)/spb(i2,i5)/spb(i2,i5)
         - spb(i1,i3) * spb(i1,i5)*spb(i1,i5) * spa(i3,i5)*spa(i3,i5)
                      /spb(i1,i2)/spb(i2,i5)/spa(i3,i4)/spa(i4,i5)/spb(i3,i5)/spb(i3,i5)
         + spa(i2,i3)*spa(i2,i3) * spa(i2,i5) * spa(i3,i5)
                      /spa(i1,i2)/spa(i1,i3)/spa(i3,i4)/spb(i3,i5)/spa(i4,i5)
         - spa(i2,i3)*spa(i2,i3) * spa(i2,i5) * spa(i3,i5)
                      /spa(i1,i2)/spa(i1,i3)/spa(i2,i4)/spb(i2,i5)/spa(i4,i5)
         - spa(i2,i3) * spa(i2,i5) * spa(i3,i5) * spb(i4,i5)
                      /spa(i1,i2)/spa(i1,i3)/spb(i2,i5)/spb(i3,i5)/spa(i4,i5)
         + ggggamgam_pmmpm_A(i1,i2,i3,i4,i5) - ggggamgam_pmmpm_A(i1,i3,i2,i4,i5)
         + ggggamgam_pmmpm_B(i1,i2,i3,i4,i5) ; }

//===

Complex ggggamgam_gen(int i1, int h1, int i2, int h2, int i3, int h3,
              int i4, int h4, int i5, int h5){
  int hsum = h1+h2+h3+h4+h5;
//  more (+)'s than (-)'s
  if (hsum == 5){ return ggggamgam_ppppp(i1,i2,i3,i4,i5); }
  else if (hsum == 3){
    if (h1 == -1){ return ggggamgam_mpppp(i1,i2,i3,i4,i5); }
    else if (h2 == -1){ return ggggamgam_mpppp(i2,i3,i1,i4,i5); }
    else if (h3 == -1){ return ggggamgam_mpppp(i3,i1,i2,i4,i5); }
    else if (h4 == -1){ return ggggamgam_pppmp(i1,i2,i3,i4,i5); }
    else { return ggggamgam_pppmp(i1,i2,i3,i5,i4); } }
  else if (hsum == 1){
    if ((h1 == -1) && (h2 == -1)){ return ggggamgam_mmppp(i1,i2,i3,i4,i5); }
    else if ((h2 == -1) && (h3 == -1)){ return ggggamgam_mmppp(i2,i3,i1,i4,i5); }
    else if ((h1 == -1) && (h3 == -1)){ return ggggamgam_mmppp(i3,i1,i2,i4,i5); }
    else if ((h4 == -1) && (h1 == -1)){ return ggggamgam_mppmp(i1,i2,i3,i4,i5); }
    else if ((h5 == -1) && (h1 == -1)){ return ggggamgam_mppmp(i1,i2,i3,i5,i4); }
    else if ((h4 == -1) && (h2 == -1)){ return ggggamgam_mppmp(i2,i3,i1,i4,i5); }
    else if ((h5 == -1) && (h2 == -1)){ return ggggamgam_mppmp(i2,i3,i1,i5,i4); }
    else if ((h4 == -1) && (h3 == -1)){ return ggggamgam_mppmp(i3,i1,i2,i4,i5); }
    else if ((h5 == -1) && (h3 == -1)){ return ggggamgam_mppmp(i3,i1,i2,i5,i4); }
    else { return ggggamgam_pppmm(i1,i2,i3,i4,i5); } }
//  more (-)'s than (+)'s
  else if (hsum == -1){
    if ((h1 == 1) && (h2 == 1)){ return ggggamgam_ppmmm(i1,i2,i3,i4,i5); }
    else if ((h2 == 1) && (h3 == 1)){ return ggggamgam_ppmmm(i2,i3,i1,i4,i5); }
    else if ((h1 == 1) && (h3 == 1)){ return ggggamgam_ppmmm(i3,i1,i2,i4,i5); }
    else if ((h4 == 1) && (h1 == 1)){ return ggggamgam_pmmpm(i1,i2,i3,i4,i5); }
    else if ((h5 == 1) && (h1 == 1)){ return ggggamgam_pmmpm(i1,i2,i3,i5,i4); }
    else if ((h4 == 1) && (h2 == 1)){ return ggggamgam_pmmpm(i2,i3,i1,i4,i5); }
    else if ((h5 == 1) && (h2 == 1)){ return ggggamgam_pmmpm(i2,i3,i1,i5,i4); }
    else if ((h4 == 1) && (h3 == 1)){ return ggggamgam_pmmpm(i3,i1,i2,i4,i5); }
    else if ((h5 == 1) && (h3 == 1)){ return ggggamgam_pmmpm(i3,i1,i2,i5,i4); }
    else { return ggggamgam_mmmpp(i1,i2,i3,i4,i5); } }
  else if (hsum == -3){
    if (h1 == 1){ return ggggamgam_pmmmm(i1,i2,i3,i4,i5); }
    else if (h2 == 1){ return ggggamgam_pmmmm(i2,i3,i1,i4,i5); }
    else if (h3 == 1){ return ggggamgam_pmmmm(i3,i1,i2,i4,i5); }
    else if (h4 == 1){ return ggggamgam_mmmpm(i1,i2,i3,i4,i5); }
    else { return ggggamgam_mmmpm(i1,i2,i3,i5,i4); } }
  else { return ggggamgam_mmmmm(i1,i2,i3,i4,i5); }
}

//================================================================

/* Version without i1, ..., i5 labels,
 for g(1) g(2) -> gamma(3) gamma(4) g(5)
 since there is only one independent color structure, 
 Tr(T^{a_1} T^{a_2} T^{a_3}) 
 (thanks to Furry's theorem, the reversed one is NEGATIVE). */

Complex gggamgamg(int h1, int h2, int h3, int h4, int h5){
 return ggggamgam_gen(1,h1,2,h2,5,h5,3,h3,4,h4);
}
