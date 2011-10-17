*********************************
* FORM version 3.2(Apr 16 2007) *
*********************************

#procedure colorbar(cf,ct,nc)
.sort
c `ct', `cf'; s `nc';
i ctmpa, ctmpf1, ctmpf2;
id `ct'(ctmpa?,ctmpf1?,ctmpf2?) = `ct'(ctmpa,ctmpf2,ctmpf1);
.sort
#endprocedure

#procedure colortrace(cf,ct,nc)
.sort
c `ct', `cf', ctmpf(c); s `nc', ctcf;
i coli1, ..., coli3, colf1, ..., colf4, cola1, ..., cola3;
#do j = 1,1
  id once   `cf'(cola1?,cola2?,cola3?)
          = 2*i_*`ct'(cola1,coli1,coli2)*
            ( `ct'(cola2,coli2,coli3)*`ct'(cola3,coli3,coli1)
             -`ct'(cola3,coli2,coli3)*`ct'(cola2,coli3,coli1));
  sum coli1, coli2, coli3;
  id once   `cf'(cola1?,cola2?,cola3?)*`ct'(cola3?,colf1?,colf2?)
          = i_*( `ct'(cola1,colf1,coli1)*`ct'(cola2,coli1,colf2)
                -`ct'(cola2,colf1,coli2)*`ct'(cola1,coli2,colf2));
  sum coli1, coli2;
  id once   `cf'(cola1?,cola2?,cola3?)*`ct'(cola2?,colf1?,colf2?)
          = i_*( `ct'(cola3,colf1,coli1)*`ct'(cola1,coli1,colf2)
                -`ct'(cola1,colf1,coli2)*`ct'(cola3,coli2,colf2));
  sum coli1, coli2;
  id once   `cf'(cola1?,cola2?,cola3?)*`ct'(cola1?,colf1?,colf2?)
          = i_*( `ct'(cola2,colf1,coli1)*`ct'(cola3,coli1,colf2)
                -`ct'(cola3,colf1,coli2)*`ct'(cola2,coli2,colf2));
  sum coli1; sum coli2;
  id   `ct'(cola1?,colf1?,colf2?)*`ct'(cola1?,colf2?,colf3?)
     = 1/2*ctcf*ctmpf(colf1,colf3);
  id   `ct'(cola1?,colf1?,colf2?)*`ct'(cola1?,colf3?,colf1?)
     = 1/2*ctcf*ctmpf(colf2,colf3);
  id   `ct'(cola1?,colf1?,colf2?)*`ct'(cola1?,colf3?,colf4?)
     = 1/2*( ctmpf(colf1,colf4)*ctmpf(colf2,colf3)
            -1/`nc'*ctmpf(colf1,colf2)*ctmpf(colf3,colf4));
  id ctmpf(colf1?,colf3?)*ctmpf(colf3?,colf2?) = ctmpf(colf1,colf2);
  id ctmpf(colf1?,colf3?)*ctmpf(colf2?,colf3?) = ctmpf(colf1,colf2);
  id ctmpf(colf3?,colf1?)*ctmpf(colf3?,colf2?) = ctmpf(colf1,colf2);
  id ctmpf(colf3?,colf1?)*ctmpf(colf2?,colf3?) = ctmpf(colf1,colf2);
  id ctmpf(colf1?,colf1?) = `nc';
  if ( count(`cf',1)>0 ) redefine j "0";
  .sort
#enddo
id ctcf = `nc'-1/`nc';
.sort
#endprocedure

* color structures

c T, F(c);
i i1, ..., i100, j1, ..., j100, a1, ..., a200;

l TCC  =  T(a1,i1,j1)*T(a1,i2,j2);

l FCC  =  F(a1,a2,a3)
         *T(a1,i1,j1)*T(a2,i2,j2)*T(a3,i3,j3);

l FFCC =  F(a1,a2,a101)*F(a101,a3,a4)
         *T(a1,i1,j1)*T(a2,i2,j2)*T(a3,i3,j3)*T(a4,i4,j4);

#call colortrace(F,T,NC)

print;
.end
