#include "DISinclu.H"
#include "Tools.H"
#include <assert.h>
#include <iostream>
namespace SHNNLO {

// h=1 corresponds to ++ or -- scattering, h=-1 corresponds to +- or -+ scattering
double DIScoeff(int h, double x, double y, double F2, double FL, double F3) {
	return F3*x*(-2 + y)*y*h - FL*x*pow(y,2) + F2*x*(2 - 2*y + pow(y,2));
}

// Jz to adjust for Jacobian transformation used for z generation
double DD(int n, double z, double x, double f, double f0, double Jz=1.) {
	if (n<0) { return f0/Jz; }
	double Lx=log(1.-x);
	double Lxnp=(n==-1?1.:(n==0?Lx:pow(Lx,n+1)));
	if (z<x) return (f0*Lxnp)/(n+1.)/Jz;
	else {
	   double Lz=log(1.-z);
	   double Lzn=(n==0?1:(n==1?Lz:pow(Lz,n)));
	   return (f0*Lxnp)/(n+1.)/Jz + (f-f0)*Lzn/(1.-z);
	}
}


  double DISinclu_0(int i, int h, double x, double y, double r, double Q, double mu) {
    double fq0 = PDF(i,x,mu); // PDF of parton i with argument x
    return DIScoeff(h,x,y,fq0,0.,fq0);
  }

  double DISinclu_1(int i, int h, double x, double y, double r, double Q, double mu) {
    double L = 2.*log(Q/mu);
    double Jz = 1.-x, z=Jz*r+x;
    double fq0 = PDF(i,x,mu); // PDF of parton i with argument x
    double fq = PDF(i,x/z,mu)/z; // PDF of parton i with argument x/z
    double fg = PDF(0,x/z,mu)/z; // PDF of gluon with argument x/z

    double F2qiqi   = -2*CF*(-3 + L - 2*z + L*z)*fq +
                      CF*(-3 + 4*L)*DD(0,z,x,fq,fq0,Jz) +
                      4*CF*DD(1,z,x,fq,fq0,Jz) -
                      2*CF*(1 + z)*log(1 - z)*fq + 
                      (CF*DD(-1,z,x,fq,fq0,Jz)*(-27 + 9*L - 2*pow(Pi,2)))/3. +
                      (2*CF*log(z)*(1 + pow(z,2)))/(-1 + z)*fq;
    
    double F2qg = 2*TF*log(1 - z)*(1 - 2*z + 2*pow(z,2)) - 2*TF*log(z)*(1 - 2*z + 2*pow(z,2)) + 2*TF*(-1 + L + 8*z - 2*L*z - 8*pow(z,2) + 2*L*pow(z,2));

    double FLqiqi = 4*CF*z;
    double FLqg = -8*TF*(-1 + z)*z;

    double F3qiqi  = -2*CF*(-2 + L - z + L*z)*fq +
                      CF*(-3 + 4*L)*DD(0,z,x,fq,fq0,Jz) +
                      4*CF*DD(1,z,x,fq,fq0,Jz) -
                      2*CF*(1 + z)*log(1 - z)*fq +
                      (CF*DD(-1,z,x,fq,fq0,Jz)*(-27 + 9*L - 2*pow(Pi,2)))/3. +
                      (2*CF*log(z)*(1 + pow(z,2)))/(-1 + z)*fq;
    double F3qg = 0.;
    
    return DIScoeff(h,x,y,F2qiqi+F2qg*fg,FLqiqi*fq+FLqg*fg,F3qiqi+F3qg*fg)*Jz;
  }

/*
  double DISinclu_my(int i, int h, double x, double y, double r, double Q, double mu) {
    // translate r to range [x:1]
    double z   = r*(1.0-x)+x;
    double fq0 = PDF(i,x,mu);
    double fq  = PDF(i,x/z,mu)/z;
    double fg  = PDF(0,x/z,mu)/z;

    double omz = 1.0-z; double omx = 1.0-x;
    double c2qxf = (fq-fq0)*( 2.0*log(omz)/omz  - 1.5/omz                              )*omx; // plus distribution part
    c2qxf       -=  fq0    *( -pow(log(omx),2)  + 1.5*log(omx)                         )    ; // remainder -fq0 \int_0^x 'DISTRIBUTION(x)' dx
    c2qxf       +=  fq     *( -(1.0+z)*log(omz) - (1.0+pow(z,2))/omz*log(z)+3.0+2.0*z  )*omx; // regular part
    c2qxf       +=  fq0    *( - pow(M_PI,2)/3.0 - 4.5                                  )    ; // delta distribution part

    // Pink book convention
    // c2qxf       +=  (fq0*(1.5+2.0*log(omx))+omx*(((1.0+pow(z,2))*fq-2.0*fq0)/omz))*2.0*log(Q/mu);     // running of PDF
    // Catani-Seymour convention
    c2qxf       +=  (-fq0*(-0.5*x*x-x-2.0*log(omx))+omx*((1.0+pow(z,2))/omz*(fq-fq0)))*2.0*log(Q/mu); // running of PDF

    c2qxf       *=  2.0*CF                                                                  ; // overall factor of 2CF

    double c2gxf = fg*2.0*TF*((pow(omz,2) + pow(z,2))*log(omz/z)-8.0*pow(z,2)+8.0*z-1.0)*omx; // gluon contribution

    double cLqxf = fq*4.0*CF*z                                                          *omx;
    double cLgxf = fg*8.0*TF*z*omz                                                      *omx;
    return DIScoeff(h,x,y,c2qxf+c2gxf,cLqxf+cLgxf,0.);
  }

*/
  
double DISinclu_2(int i, int h, double x, double y, double r, double Q, double mu) {
	double L = 2.*log(Q/mu);
	double Jz = 1.-x, z=Jz*r+x;
	double fq0 = PDF(i,x,mu); // PDF of parton i with argument x
	double fq = PDF(i,x/z,mu)/z; // PDF of parton i with argument x/z
	double fqb = PDF(-i,x/z,mu)/z; // PDF of parton -i with argument x/z
	double fqj = 0.; // sum of PDF of quark j with j!=i and argument x/z
	double fg = PDF(0,x/z,mu)/z; // PDF of gluon with argument x/z
	for (int id=1; id<=Nf; id++) { if (id!=i&&id!=-i) { fqj += PDF(id,x/z,mu)+PDF(-id,x/z,mu); } }
	fqj = fqj/z;

	// -------- F2 --------
	double C2PSQ = (CF*TF*(474 + 344/z - 1266*z - 216*(1 + z)*Li3(1 - z) + 448*pow(z,2) + (6*pow(Pi,2)*(-8 - 3*z - 9*pow(z,2) + 4*pow(z,3)))/z + 
        (36*Li2(1 - z)*(4 + 3*z - 3*pow(z,2) + 8*pow(z,3)))/z - (144*Li2(-z)*pow(1 + z,3))/z))/27. - 
   (4*CF*TF*log(z)*(18*z*(1 + z)*Li2(1 - z) + z*(-126 + 66*z + 3*(1 + z)*pow(Pi,2) + 32*pow(z,2)) + 12*log(1 + z)*pow(1 + z,3)))/(9.*z) - 
   (2*CF*TF*(-4 - 3*z - 6*z*(1 + z)*log(z) + 3*pow(z,2) + 4*pow(z,3))*pow(log(1 - z),2))/(3.*z) - 
   (CF*TF*(3 - 45*z + 32*pow(z,2))*pow(log(z),2))/3. + (8*CF*TF*log(1 - z)*
      (13 - 39*z + 9*z*(1 + z)*Li2(1 - z) + 30*pow(z,2) - 4*pow(z,3) + 18*log(z)*pow(z,3) - 9*z*(1 + z)*pow(log(z),2)))/(9.*z) + 
   (10*CF*TF*(1 + z)*pow(log(z),3))/3.; // already divided by 2 nf
	//
	double C2NSpQ = (-2*CF*(11*CA + 27*CF - 4*nf*TF)*DD(2,z,x,1.,fq0/fq,Jz))/3. + 8*DD(3,z,x,1.,fq0/fq,Jz)*pow(CF,2) - 
   (CF*DD(1,z,x,1.,fq0/fq,Jz)*(116*nf*TF + CA*(-367 + 12*pow(Pi,2)) + 3*CF*(81 + 16*pow(Pi,2))))/9. + 
   (CF*DD(0,z,x,1.,fq0/fq,Jz)*(4*nf*TF*(247 - 12*pow(Pi,2)) + 27*CF*(51 - 16*Zeta3 + 12*pow(Pi,2)) + CA*(-3155 + 2160*Zeta3 + 132*pow(Pi,2))))/
    54. + (CF*DD(-1,z,x,1.,fq0/fq,Jz)*(CA*(-27325 + 16800*Zeta3 - 5020*pow(Pi,2) + 142*pow(Pi,4)) + 
        5*(4*nf*TF*(457 + 48*Zeta3 + 76*pow(Pi,2)) + 3*CF*(993 - 1872*Zeta3 + 276*pow(Pi,2) + 4*pow(Pi,4)))))/360. - 
   (CF*log(z)*(12*CA*z - 24*CF*z - 12*CA*log(1 + z) + 24*CF*log(1 + z) + 12*CA*z*log(1 + z) - 24*CF*z*log(1 + z) - 603*CA*pow(z,2) + 
        501*CF*pow(z,2) + 120*nf*TF*pow(z,2) + 240*CA*Li2(z)*pow(z,2) - 300*CF*Li2(z)*pow(z,2) - 300*CA*log(1 + z)*pow(z,2) + 
        600*CF*log(1 + z)*pow(z,2) + 50*CF*pow(Pi,2)*pow(z,2) + 120*(CA - 2*CF)*Li2(-z)*pow(z,2)*(2 - 3*z + 4*pow(z,2)) + 
        60*Li2(1 - z)*pow(z,2)*(CF*(-1 + pow(z,2)) + CA*(1 + pow(z,2))) + 870*CA*pow(z,3) - 630*CF*pow(z,3) - 120*nf*TF*pow(z,3) - 
        360*CA*Li2(z)*pow(z,3) + 720*CF*Li2(z)*pow(z,3) + 240*CA*log(1 + z)*pow(z,3) - 480*CF*log(1 + z)*pow(z,3) - 1582*CA*pow(z,4) + 
        1284*CF*pow(z,4) + 380*nf*TF*pow(z,4) + 480*CA*Li2(z)*pow(z,4) - 780*CF*Li2(z)*pow(z,4) + 420*CA*log(1 + z)*pow(z,4) - 
        840*CF*log(1 + z)*pow(z,4) + 70*CF*pow(Pi,2)*pow(z,4) + 108*CA*pow(z,5) - 216*CF*pow(z,5) - 252*CA*log(1 + z)*pow(z,5) + 
        504*CF*log(1 + z)*pow(z,5) - 108*CA*log(1 + z)*pow(z,6) + 216*CF*log(1 + z)*pow(z,6)))/(15.*(-1 + z)*pow(z,2)) - 
   (CF*(90*Li2(1 - z)*pow(z,2)*(-4*nf*TF*(1 + pow(z,2)) + CA*(17 + 5*pow(z,2)) + 6*CF*(-5 - 4*z + 6*pow(z,2))) + 
        z*(-108*CA + 216*CF + 1630*CA*z + 810*CF*z - 1580*nf*TF*z - 2970*CA*z*Zeta3 + 2700*CF*z*Zeta3 - 2160*CA*z*Li3(z) + 2700*CF*z*Li3(z) - 
           210*CA*z*pow(Pi,2) - 630*CF*z*pow(Pi,2) + 60*nf*TF*z*pow(Pi,2) + 17211*CA*pow(z,2) - 15687*CF*pow(z,2) - 3300*nf*TF*pow(z,2) - 
           3240*CA*Zeta3*pow(z,2) + 6480*CF*Zeta3*pow(z,2) + 3240*CA*Li3(z)*pow(z,2) - 6480*CF*Li3(z)*pow(z,2) - 180*CA*pow(Pi,2)*pow(z,2) + 
           360*CF*pow(Pi,2)*pow(z,2) - 2160*(CA - 2*CF)*z*Li3(-z)*(2 - 3*z + 4*pow(z,2)) - 
           540*z*Li3(1 - z)*(3*CA*(1 + pow(z,2)) - 2*(CF + 2*CF*pow(z,2))) - 19705*CA*pow(z,3) + 16605*CF*pow(z,3) + 4880*nf*TF*pow(z,3) + 
           4590*CA*Zeta3*pow(z,3) - 2700*CF*Zeta3*pow(z,3) - 4320*CA*Li3(z)*pow(z,3) + 7020*CF*Li3(z)*pow(z,3) + 930*CA*pow(Pi,2)*pow(z,3) - 
           810*CF*pow(Pi,2)*pow(z,3) - 60*nf*TF*pow(Pi,2)*pow(z,3) + 972*CA*pow(z,4) - 1944*CF*pow(z,4) - 378*CA*pow(Pi,2)*pow(z,4) + 
           756*CF*pow(Pi,2)*pow(z,4) - 162*CA*pow(Pi,2)*pow(z,5) + 324*CF*pow(Pi,2)*pow(z,5)) - 
        108*(CA - 2*CF)*Li2(-z)*(1 - z + 25*pow(z,2) - 20*pow(z,3) - 35*pow(z,4) + 21*pow(z,5) + 9*pow(z,6))))/(135.*(-1 + z)*pow(z,2)) + 
   (CF*((-1 + z)*(11*CA*(1 + z) - 4*nf*TF*(1 + z) + 6*CF*(7 + 9*z)) + 12*CF*log(z)*(3 + 4*pow(z,2)))*pow(log(1 - z),2))/(3.*(-1 + z)) - 
   4*(1 + z)*pow(CF,2)*pow(log(1 - z),3) + (CF*(CA*(335 - 120*z + 695*pow(z,2) - 252*pow(z,3) - 108*pow(z,4)) + 
        2*(-50*nf*TF*(1 + pow(z,2)) + 3*CF*(-65 + 10*z - 50*pow(z,2) + 84*pow(z,3) + 36*pow(z,4))))*pow(log(z),2))/(30.*(-1 + z)) - 
   (2*CF*log(1 - z)*(-2*(-1 + z)*(2*nf*TF*(8 + 17*z) - 9*CF*(8 - 21*z + (-1 + z)*pow(Pi,2)) + CA*(4 + z*(-185 + 9*pow(Pi,2)))) + 
        18*Li2(1 - z)*(CA + CA*pow(z,2) - 2*CF*pow(z,2)) + 
        12*log(z)*(-2*nf*TF*(1 + pow(z,2)) + CA*(7 + 4*pow(z,2)) + 3*CF*(-2 - 2*z + 7*pow(z,2))) + 
        9*(CF*(-1 + 12*z - 5*pow(z,2)) + CA*(5 - 6*z + 9*pow(z,2)))*pow(log(z),2)))/(9.*(-1 + z)) + 
   (CF*(3*CA*(1 + pow(z,2)) + CF*(-1 + 9*pow(z,2)))*pow(log(z),3))/(3.*(-1 + z));
	//
	double C2NSmQ = (4*CF*(-CA + 2*CF)*log(1 + z)*pow(Pi,2)*(1 + 3*z + 3*pow(z,2)))/(3.*(1 + z)) + 
   (CF*(-CA/2. + CF)*log(z)*(-26 - 8/z - 106*z + 72*pow(z,2) - (40*(1 + pow(z,2)))/(1 + z) + (40*Li2(1 - z)*(1 + pow(z,2)))/(1 + z) + 
        (80*Li2(-z)*(1 + pow(z,2)))/(1 + z) + (40*Li2(z)*(1 + pow(z,2)))/(1 + z) + 8*log(1 + z)*(20 + 20*z + 1/pow(z,2) + 30*pow(z,2) - 9*pow(z,3))
        ))/5. + (CF*(-CA/2. + CF)*(-120*(1 + 5*z)*Zeta3 + 120*(1 + z)*Li2(1 - z) + 120*(1 + 5*z)*Li3(-z) + 240*(1 + 5*z)*Li3(1/(1 + z)) - 
        (240*Li3(1 - z)*(1 + pow(z,2)))/(1 + z) + (120*Li3(-z)*(1 + pow(z,2)))/(1 + z) - (120*Li3(z)*(1 + pow(z,2)))/(1 + z) + 
        (240*Li3(1/(1 + z))*(1 + pow(z,2)))/(1 + z) + (240*Li3((1 - z)/(1 + z))*(1 + pow(z,2)))/(1 + z) - 
        (240*Li3((-1 + z)/(1 + z))*(1 + pow(z,2)))/(1 + z) + 3*(-162 + 8/z + 82*z + 72*pow(z,2)) + 
        24*Li2(-z)*(20 + 20*z + 1/pow(z,2) + 30*pow(z,2) - 9*pow(z,3)) - 2*pow(Pi,2)*(5 - 25*z - 60*pow(z,2) + 18*pow(z,3))))/15. + 
   (2*CF*(-CA + 2*CF)*(-5 - 25*z - 50*pow(z,2) + 10*log(1 + z)*(3 + 3*z + 5*pow(z,2)) - 21*pow(z,3) + 9*pow(z,4))*pow(log(z),2))/(5.*(1 + z)) + 
   (2*(CA - 2*CF)*CF*log(1 - z)*(-12 + pow(Pi,2) + 12*pow(z,2) + pow(Pi,2)*pow(z,2) + 12*Li2(-z)*(1 + pow(z,2)) + 
        6*log(z)*(2*log(1 + z)*(1 + pow(z,2)) - pow(1 + z,2)) - 6*(1 + pow(z,2))*pow(log(z),2)))/(3.*(1 + z)) - 
   (CF*(-CA + 2*CF)*(1 + pow(z,2))*pow(log(z),3))/(1 + z) - (8*CF*(-CA + 2*CF)*(1 + 3*z + 3*pow(z,2))*pow(log(1 + z),3))/(3.*(1 + z));
	//
	double C2G = (-4*(CA - 2*CF)*TF*log(1 + z)*pow(Pi,2)*pow(1 + z,2))/3. + 
   TF*log(z)*((8*CA*(-3 + z)*z*pow(Pi,2))/3. - 4*CF*Li2(z)*pow(-1 + z,2) + 2*CF*(2*Li2(1 - z) - 8*Li2(-z) + pow(Pi,2))*pow(-1 + z,2) + 
      (CF*(-236 - 8/z + 339*z - 648*pow(z,2)))/15. + CA*(58 + (584*z)/3. - (2090*pow(z,2))/9.) + 8*CA*Li2(1 - z)*pow(z,2) + 
      12*CF*Li2(1 - z)*pow(z,2) - 4*CA*Li2(z)*pow(z,2) + 12*CF*Li2(z)*pow(z,2) + (10*CF*pow(Pi,2)*pow(z,2))/3. + 
      (8*log(1 + z)*((5*CA*(-2 - 9*z + 10*pow(z,3)))/z + CF*(90 + 40*z + 1/pow(z,2) + 36*pow(z,3))))/15. - 8*CA*Li2(1 - z)*pow(1 + z,2) + 
      8*CA*Li2(-z)*pow(1 + z,2) + 4*CA*Li2(z)*pow(1 + z,2)) + 
   TF*((239*CA)/9. - (647*CF)/15. + (344*CA)/(27.*z) + (8*CF)/(15.*z) + (1072*CA*z)/9. + (239*CF*z)/5. + 4*CF*z*(16 + 9*z)*Zeta3 + 
      2*CF*(-5 + 12*z)*Li2(1 - z) + 64*CF*z*Li3(-z) + 48*CF*Zeta3*pow(-1 + z,2) - 16*CF*Li3(1 - z)*pow(-1 + z,2) + 48*CF*Li3(-z)*pow(-1 + z,2) + 
      4*CF*Li3(z)*pow(-1 + z,2) + (13*CF*pow(Pi,2)*pow(-1 + z,2))/3. - (4493*CA*pow(z,2))/27. - (36*CF*pow(z,2))/5. - 4*CA*Zeta3*pow(z,2) + 
      12*CF*Zeta3*pow(z,2) - 8*CF*Li3(1 - z)*pow(z,2) + 16*CA*Li3(-z)*pow(z,2) + 4*CA*Li3(z)*pow(z,2) - 12*CF*Li3(z)*pow(z,2) + 
      8*CA*Li3((1 - z)/(1 + z))*(1 + 2*z + 2*pow(z,2)) - 8*CA*Li3((-1 + z)/(1 + z))*(1 + 2*z + 2*pow(z,2)) - 2*CA*Zeta3*(5 + 6*z + 6*pow(z,2)) + 
      CA*Li3(1 - z)*(-4 - 72*z + 8*pow(z,2)) + (CF*z*pow(Pi,2)*(-50 + 345*z + 144*pow(z,2)))/45. + 
      (8*CA*Li2(-z)*(-2 - 9*z + 10*pow(z,3)))/(3.*z) + (8*CF*Li2(-z)*(90 + 40*z + 1/pow(z,2) + 36*pow(z,3)))/15. + 
      (4*CA*Li2(1 - z)*(4 + 3*z - 48*pow(z,2) + 44*pow(z,3)))/(3.*z) + (2*CA*pow(Pi,2)*(-8 + 3*z - 60*pow(z,2) + 67*pow(z,3)))/(9.*z) + 
      20*CA*Zeta3*pow(1 + z,2) - 32*CF*Zeta3*pow(1 + z,2) - 8*CA*Li3(-z)*pow(1 + z,2) - 4*CA*Li3(z)*pow(1 + z,2) - 
      16*CA*Li3(1/(1 + z))*pow(1 + z,2) + 32*CF*Li3(1/(1 + z))*pow(1 + z,2)) - 
   (TF*(3*CF*z*(13 - 40*z + 36*pow(z,2)) + 12*z*log(z)*(2*CA*(-3 + z)*z + CF*(2 - 4*z + 5*pow(z,2))) + 
        2*CA*(-4 + 3*z - 54*pow(z,2) + 61*pow(z,3)))*pow(log(1 - z),2))/(3.*z) + (2*(CA + 5*CF)*TF*(1 - 2*z + 2*pow(z,2))*pow(log(1 - z),3))/3. + 
   TF*(CA*(-1 + 88*z - (194*pow(z,2))/3.) + CF*(-1.5 + (22*z)/3. - 36*pow(z,2) - (48*pow(z,3))/5.) + 
      log(1 + z)*(4*CA*(1 + 2*z + 3*pow(z,2)) + 8*CF*pow(1 + z,2)))*pow(log(z),2) + 
   TF*log(1 - z)*(2*CF*z*(-6 + 5*z) + 14*CF*pow(-1 + z,2) + 8*CF*Li2(1 - z)*pow(-1 + z,2) + CA*Li2(1 - z)*(4 + 40*z - 8*pow(z,2)) - 
      (8*CF*pow(Pi,2)*pow(z,2))/3. - 8*CA*Li2(-z)*(1 + 2*z + 2*pow(z,2)) - (2*CA*pow(Pi,2)*(5 - 6*z + 8*pow(z,2)))/3. + 
      4*log(z)*(-2*CA*log(1 + z)*(1 + 2*z + 2*pow(z,2)) + 2*CF*(2 - 7*z + 9*pow(z,2)) + CA*(2 - 36*z + 37*pow(z,2))) + 
      (2*CA*(52 - 93*z - 681*pow(z,2) + 785*pow(z,3)))/(9.*z) + 2*(CA*(-1 - 14*z + 4*pow(z,2)) + CF*(3 - 6*z + 14*pow(z,2)))*pow(log(z),2)) + 
   (TF*(2*CA*(5 + 14*z) - 5*CF*(1 - 2*z + 4*pow(z,2)))*pow(log(z),3))/3. + (8*(CA - 2*CF)*TF*pow(1 + z,2)*pow(log(1 + z),3))/3.; // already divided by 2 nf
	// assemble for each channel
	double F2qiqi = C2NSpQ+C2PSQ;
	double F2qiqbi = C2NSmQ+C2PSQ;
	double F2qiqj = C2PSQ;
	double F2qg = C2G;
	// add mu dependent part
	if (L!=0) {
	   F2qiqi += (4*CF*L*(-11*CA + 3*CF*(-3 + 4*L) + 4*nf*TF)*DD(1,z,x,1.,fq0/fq,Jz))/3. + 24*L*DD(2,z,x,1.,fq0/fq,Jz)*pow(CF,2) - 
   (CF*L*DD(0,z,x,1.,fq0/fq,Jz)*(4*(29 - 6*L)*nf*TF + CA*(-367 + 66*L + 12*pow(Pi,2)) + CF*(405 - 108*L + 48*pow(Pi,2))))/9. - 
   (CF*L*DD(-1,z,x,1.,fq0/fq,Jz)*(CA*(-645 + 99*L + 216*Zeta3 - 88*pow(Pi,2)) + 4*nf*TF*(57 - 9*L + 8*pow(Pi,2)) + 
        3*CF*(153 - 27*L - 240*Zeta3 + 12*pow(Pi,2) + 8*L*pow(Pi,2))))/18. + 
   (2*CF*L*log(z)*(-4*CA*(7 + 4*pow(z,2)) + 3*CF*(10 + L + 8*z - 12*pow(z,2) + 3*L*pow(z,2)) + 
        2*TF*(3*L*(-1 + pow(z,2)) + 4*(3*(-1 + z)*pow(z,2) + nf*(1 + pow(z,2))))))/(3.*(-1 + z)) - 
   (CF*L*(-104*TF - 24*L*TF + 92*CA*z - 198*CF*z - 33*CA*L*z + 90*CF*L*z + 312*TF*z - 18*L*TF*z - 64*nf*TF*z + 12*L*nf*TF*z - 
        72*TF*z*(1 + z)*Li2(1 - z) + 36*CF*z*(1 + z)*Li2(z) - 6*CA*z*pow(Pi,2) - 30*CF*z*pow(Pi,2) + 506*CA*pow(z,2) - 396*CF*pow(z,2) - 
        33*CA*L*pow(z,2) + 18*CF*L*pow(z,2) - 240*TF*pow(z,2) + 18*L*TF*pow(z,2) - 136*nf*TF*pow(z,2) + 12*L*nf*TF*pow(z,2) - 
        6*CA*pow(Pi,2)*pow(z,2) - 30*CF*pow(Pi,2)*pow(z,2) + 32*TF*pow(z,3) + 24*L*TF*pow(z,3)))/(9.*z) - 
   (2*CF*L*log(1 - z)*(-12*z*log(z)*(TF*(-1 + pow(z,2)) + 3*CF*(1 + pow(z,2))) + 
        (-1 + z)*(-(z*(11*CA*(1 + z) - 12*CF*(-2 + L - 3*z + L*z))) + TF*(-8 + (-6 + 4*nf)*z + (6 + 4*nf)*pow(z,2) + 8*pow(z,3)))))/(3.*(-1 + z)*z)
     - 12*L*(1 + z)*pow(CF,2)*pow(log(1 - z),2) - (2*CF*L*(CA - 4*TF + CA*pow(z,2) + 4*CF*pow(z,2) + 4*TF*pow(z,2))*pow(log(z),2))/(-1 + z);
	   F2qiqbi += (2*CF*L*(-36*(CA - 2*CF) - 156*TF + 9*L*TF + (52*TF)/z + (12*L*TF)/z + 36*(CA - 2*CF)*z + 120*TF*z - 9*L*TF*z + 36*TF*(1 + z)*Li2(1 - z) - 
        16*TF*pow(z,2) - 12*L*TF*pow(z,2) + (36*(CA - 2*CF)*Li2(-z)*(1 + pow(z,2)))/(1 + z) + (3*(CA - 2*CF)*pow(Pi,2)*(1 + pow(z,2)))/(1 + z)))/9.
     + (4*CF*L*log(z)*(2*(CA - 2*CF)*log(1 + z)*(1 + pow(z,2)) + (1 + z)*(-(CA*(1 + z)) + 2*CF*(1 + z) + TF*(L + L*z + 4*pow(z,2)))))/(1 + z) - 
   (4*CF*L*TF*log(1 - z)*(-4 - 3*z - 6*z*(1 + z)*log(z) + 3*pow(z,2) + 4*pow(z,3)))/(3.*z) - 
   (2*CF*L*(CA*(1 + pow(z,2)) - 2*CF*(1 + pow(z,2)) + 4*TF*pow(1 + z,2))*pow(log(z),2))/(1 + z);
	   F2qiqj += 4*CF*L*TF*log(z)*(L + L*z + 4*pow(z,2)) - (2*CF*L*TF*(-36*z*(1 + z)*Li2(1 - z) + 
        (-1 + z)*(4*(13 - 26*z + 4*pow(z,2)) + 3*L*(4 + 7*z + 4*pow(z,2)))))/(9.*z) - 
   (4*CF*L*TF*log(1 - z)*(-4 - 3*z - 6*z*(1 + z)*log(z) + 3*pow(z,2) + 4*pow(z,3)))/(3.*z) - 8*CF*L*TF*(1 + z)*pow(log(z),2);
	   F2qg += -2*L*TF*log(z)*(-2*CA*(L + 4*L*z + z*(-24 + 25*z)) + 4*CA*log(1 + z)*(1 + 2*z + 2*pow(z,2)) + 
      CF*(-2 + L + 12*z - 2*L*z - 20*pow(z,2) + 4*L*pow(z,2))) + 
   (L*TF*(104*CA + 24*CA*L - 330*CA*z + 162*CF*z + 18*CA*L*z - 9*CF*L*z + 72*CA*z*(1 + 4*z)*Li2(1 - z) - 36*CF*z*Li2(z) - 12*CA*z*pow(Pi,2) - 
        18*CF*z*pow(Pi,2) - 552*CA*pow(z,2) - 306*CF*pow(z,2) + 144*CA*L*pow(z,2) + 36*CF*L*pow(z,2) + 72*CF*Li2(z)*pow(z,2) + 
        36*CF*pow(Pi,2)*pow(z,2) - 72*CA*z*Li2(-z)*(1 + 2*z + 2*pow(z,2)) + 814*CA*pow(z,3) + 72*CF*pow(z,3) - 186*CA*L*pow(z,3) - 
        24*CA*pow(Pi,2)*pow(z,3) - 48*CF*pow(Pi,2)*pow(z,3)))/(9.*z) + 
   (2*L*TF*log(1 - z)*(-24*z*log(z)*(CA*(-3 + z)*z + CF*(1 - 2*z + 2*pow(z,2))) + 3*CF*z*(-7 + 24*z - 20*pow(z,2) + L*(2 - 4*z + 4*pow(z,2))) + 
        2*CA*(4 + 3*(-1 + L)*z - 6*(-10 + L)*pow(z,2) + (-67 + 6*L)*pow(z,3))))/(3.*z) + 
   4*(CA + 2*CF)*L*TF*(1 - 2*z + 2*pow(z,2))*pow(log(1 - z),2) + 4*L*TF*(-2*CA + CF - 6*CA*z - 2*CF*z + 4*CF*pow(z,2))*pow(log(z),2);
	}

	// -------- FL --------
	double CLPSQ = -16*CF*TF*log(z)*(-1 + z + 2*pow(z,2)) + (16*CF*TF*(-1 + z)*log(1 - z)*(-1 + 2*z + 2*pow(z,2)))/(3.*z) + 
   (8*CF*TF*(-2 + 6*z + 18*Li2(z)*pow(z,2) - 3*(8 + pow(Pi,2))*pow(z,2) + 20*pow(z,3)))/(9.*z) + 16*CF*TF*z*pow(log(z),2); // already divided by 2 nf
	//
	double CLNSpQ = (4*CF*(-23*CA*z + 3*CF*(2 + 7*z))*log(1 - z))/3. - (16*CF*nf*TF*z*log(pow(z,2)/(1 - z)))/3. - (16*(CA - 2*CF)*CF*z*log(1 + z)*pow(Pi,2))/3. + 
   (8*(CA - 2*CF)*CF*z*log(1 - pow(z,2))*pow(Pi,2))/3. + log(z)*
    ((4*CF*z*(44*CF - (6*CF)/z - 12*(CA - 2*CF)*Li2(-z) - 12*(CA - 2*CF)*Li2(z) + 
           (4*(CA - 2*CF)*(6 - 3*z + 47*pow(z,2) - 9*pow(z,3)))/(5.*pow(z,2))))/3. - 
      (16*(CA - 2*CF)*CF*log(1 + z)*(2 + 10*pow(z,2) + 5*pow(z,3) - 3*pow(z,5)))/(5.*pow(z,2))) + 
   (4*CF*z*((355*CF)/6. - (50*nf*TF)/3. - (13*CF)/z + (4*nf*TF)/z + 6*CF*Li2(z) + 12*(CA - 2*CF)*Li3(-z) + 12*(CA - 2*CF)*Li3(z) - 
        24*(CA - 2*CF)*Li3(1/(1 + z)) - 2*(CA - 2*CF)*pow(Pi,2) - 3*CF*pow(Pi,2) + (6*(CA - 2*CF)*pow(Pi,2)*pow(z,2))/5. - 
        ((CA - 2*CF)*(144 + 294*z - 1729*pow(z,2) + 216*pow(z,3)))/(30.*pow(z,2)) - 
        (12*(CA - 2*CF)*Li2(-z)*(2 + 10*pow(z,2) + 5*pow(z,3) - 3*pow(z,5)))/(5.*pow(z,3))))/3. - 
   (8*(CA - 2*CF)*CF*z*(-5 + 5*log(1 - pow(z,2)) + 3*pow(z,2))*pow(log(z),2))/5. + 8*z*pow(CF,2)*pow(log(z/(1 - z)),2) + 
   (16*(CA - 2*CF)*CF*z*pow(log(1 + z),3))/3.;
	//
	double CLNSmQ = 0.;
	//
	double CLG = (8*TF*log(1 - z)*(6*(CF + 2*CA*(-3 + z))*log(z)*pow(z,2) + (-1 + z)*(-3*CF*z*(1 + 4*z) + CA*(-2 + 4*z + 58*pow(z,2)))))/(3.*z) + 
   (16*TF*(-45*(4*CA - CF)*Li2(1 - z)*pow(z,3) + z*(5*CA*(-1 + 3*z + 51*pow(z,2) + (-53 + 3*pow(Pi,2))*pow(z,3)) + 
           CF*(6 - 24*z - (171 + 5*pow(Pi,2))*pow(z,2) + 189*pow(z,3) + 6*pow(Pi,2)*pow(z,4))) + 
        6*Li2(-z)*(15*CA*(1 + z)*pow(z,3) + CF*(1 - 5*pow(z,3) + 6*pow(z,5)))))/(45.*pow(z,2)) + 
   (8*TF*log(z)*(z*(30*CA*z*(1 + 8*z - 13*pow(z,2)) + CF*(-4 - 13*z - 78*pow(z,2) + 36*pow(z,3))) + 
        4*log(1 + z)*(15*CA*(1 + z)*pow(z,3) + CF*(1 - 5*pow(z,3) + 6*pow(z,5)))))/(15.*pow(z,2)) - 16*CA*TF*(-1 + z)*z*pow(log(1 - z),2) + 
   TF*(48*CA*z - (32*CF*z*(5 + 3*pow(z,2)))/15.)*pow(log(z),2); // already divided by 2 nf
	// assemble for each channel
	double FLqiqi = CLNSpQ+CLPSQ;
	double FLqiqbi = CLNSmQ+CLPSQ;
	double FLqiqj = CLPSQ;
	double FLqg = CLG;
	// add mu dependent part
	if (L!=0) {
	   FLqiqi += -8*CF*L*(CF + 2*TF)*z*log(z) + 16*L*z*log(1 - z)*pow(CF,2) + 
   (4*CF*L*(z*(-11*CA*z + 3*CF*(2 + z)) + 4*TF*(1 - 3*z + nf*pow(z,2) + 2*pow(z,3))))/(3.*z);
	   FLqiqbi += -16*CF*L*TF*z*log(z) + (16*CF*L*TF*(1 - 3*z + 2*pow(z,3)))/(3.*z);
	   FLqiqj += -16*CF*L*TF*z*log(z) + (16*CF*L*TF*(1 - 3*z + 2*pow(z,3)))/(3.*z);
	   FLqg += -32*CA*L*TF*(-1 + z)*z*log(1 - z) - 16*(4*CA - CF)*L*TF*z*log(z) + (8*L*TF*(-1 + z)*(-3*CF*z*(1 + 2*z) + CA*(-2 + 4*z + 34*pow(z,2))))/(3.*z);
	}

	// -------- F3 --------
	double C3PSQ = 0.;
	//
	double C3NSpQ = C2NSpQ;
	C3NSpQ += (CF*(180*CF*(1 + z)*Li2(1 - z)*pow(z,2) + z*(36*CA - 72*CF + 701*CA*z - 2187*CF*z + 140*nf*TF*z - 360*CA*z*Zeta3 + 720*CF*z*Zeta3 - 
           720*(CA - 2*CF)*z*(-1 + 3*z)*Li3(-z) - 360*(CA - 2*CF)*z*(-1 + 3*z)*Li3(z) + 30*CA*z*pow(Pi,2) - 3211*CA*pow(z,2) + 4437*CF*pow(z,2) + 
           620*nf*TF*pow(z,2) + 1080*CA*Zeta3*pow(z,2) - 2160*CF*Zeta3*pow(z,2) + 90*CA*pow(Pi,2)*pow(z,2) - 120*CF*pow(Pi,2)*pow(z,2) + 
           324*CA*pow(z,3) - 648*CF*pow(z,3) - 150*CA*pow(Pi,2)*pow(z,3) + 300*CF*pow(Pi,2)*pow(z,3) - 54*CA*pow(Pi,2)*pow(z,4) + 
           108*CF*pow(Pi,2)*pow(z,4)) - 36*(CA - 2*CF)*Li2(-z)*(-1 - 5*z - 30*pow(z,2) - 10*pow(z,3) + 25*pow(z,4) + 9*pow(z,5))))/(45.*pow(z,2))\
    + (4*CF*log(z)*(20*nf*TF*(1 + z) + 60*CF*(1 - 3*z)*Li2(-z) + 30*CA*(-1 + 3*z)*Li2(-z) + 60*CF*(1 - 3*z)*Li2(z) + 30*CA*(-1 + 3*z)*Li2(z) + 
        (3*CF*(2 + 4*z + 49*pow(z,2) - 18*pow(z,3)))/z + (CA*(-3 - 31*z - 121*pow(z,2) + 27*pow(z,3)))/z - 
        (3*(CA - 2*CF)*log(1 + z)*(-1 - 5*z - 30*pow(z,2) - 10*pow(z,3) + 25*pow(z,4) + 9*pow(z,5)))/pow(z,2)))/15. - 
   4*(1 + z)*pow(CF,2)*pow(log(1 - z),2) - (2*CF*(2*CF*z*(-10 + 25*z + 9*pow(z,2)) + CA*(5 + 15*z - 25*pow(z,2) - 9*pow(z,3)))*pow(log(z),2))/5. + 
   (2*CF*log(1 - z)*(63*CF - 4*nf*TF - 117*CF*z - 4*nf*TF*z + 18*CF*(1 + z)*log(z) - 4*CF*pow(Pi,2) + 12*CF*z*pow(Pi,2) + 
        CA*(-25 + 71*z + (2 - 6*z)*pow(Pi,2)) + 6*(CA - 2*CF)*(-1 + 3*z)*pow(log(z),2)))/3.;
	//
	double C3NSmQ = C2NSmQ;
        C3NSmQ += (-4*CF*(-CA + 2*CF)*(1 + 3*z)*log(1 + z)*pow(Pi,2))/3. + (4*CF*(-CA/2. + CF)*
      (78 - 2/z - 58*z + 20*(1 + 3*z)*Zeta3 - 20*(1 + 3*z)*Li3(-z) - 40*(1 + 3*z)*Li3(1/(1 + z)) - 18*pow(z,2) + 
        10*Li2(-z)*(-4 + 1/z - 4*z - 1/(5.*pow(z,2)) - 5*pow(z,2) + (9*pow(z,3))/5.) + (pow(Pi,2)*(5 - 15*z - 25*pow(z,2) + 9*pow(z,3)))/3.))/5. + 
   (4*CF*(-CA + 2*CF)*log(z)*(z + 7*pow(z,2) + 37*pow(z,3) - 9*pow(z,4) + 
        log(1 + z)*(-1 + 5*z - 20*pow(z,2) - 20*pow(z,3) - 25*pow(z,4) + 9*pow(z,5))))/(5.*pow(z,2)) + 
   (2*(CA - 2*CF)*CF*(5 - 15*z + 10*(1 + 3*z)*log(1 + z) - 25*pow(z,2) + 9*pow(z,3))*pow(log(z),2))/5. + 
   (8*CF*(-CA + 2*CF)*(1 + 3*z)*pow(log(1 + z),3))/3.;
	//
	double C3G = 0.;
	// assemble each channel
	double F3qiqi = C3NSpQ+C3PSQ;
	double F3qiqbi = C3NSmQ+C3PSQ;
	double F3qiqj = C3PSQ;
	double F3qg = C3G;
	// add mu dependent part
	if (L!=0) {
	   F3qiqi += (4*CF*L*(-11*CA + 3*CF*(-3 + 4*L) + 4*nf*TF)*DD(1,z,x,1.,fq0/fq,Jz))/3. + 24*L*DD(2,z,x,1.,fq0/fq,Jz)*pow(CF,2) - 
   (CF*L*DD(0,z,x,1.,fq0/fq,Jz)*(4*(29 - 6*L)*nf*TF + CA*(-367 + 66*L + 12*pow(Pi,2)) + CF*(405 - 108*L + 48*pow(Pi,2))))/9. - 
   (CF*L*DD(-1,z,x,1.,fq0/fq,Jz)*(CA*(-645 + 99*L + 216*Zeta3 - 88*pow(Pi,2)) + 4*nf*TF*(57 - 9*L + 8*pow(Pi,2)) + 
        3*CF*(153 - 27*L - 240*Zeta3 + 12*pow(Pi,2) + 8*L*pow(Pi,2))))/18. + 
   (2*CF*L*log(1 - z)*((-1 + z)*(11*CA*(1 + z) - 4*(nf*TF*(1 + z) + 3*CF*(-1 + L - 2*z + L*z))) + 36*CF*log(z)*(1 + pow(z,2))))/(3.*(-1 + z)) + 
   (2*CF*L*log(z)*(-4*CA*(7 + 4*pow(z,2)) + 3*CF*(8 + L + 8*z - 10*pow(z,2) + 3*L*pow(z,2)) + 
        2*TF*(-3 - 12*z + 3*L*(-1 + pow(z,2)) + 7*pow(z,2) + 4*nf*(1 + pow(z,2)) + 8*pow(z,3))))/(3.*(-1 + z)) - 
   (CF*L*(36*CF*z*(1 + z)*Li2(z) - z*(6*CF*(24 + 57*z - 3*L*(5 + z) + 5*(1 + z)*pow(Pi,2)) + 
           CA*(-26 - 440*z + 33*L*(1 + z) + 6*(1 + z)*pow(Pi,2))) + 
        2*TF*(-4*(10 + (-9 + 5*nf)*z + (27 + 14*nf)*pow(z,2) - 28*pow(z,3)) + 3*L*(-4 + (-3 + 2*nf)*z + (3 + 2*nf)*pow(z,2) + 4*pow(z,3)))))/(9.*z)
     - 12*L*(1 + z)*pow(CF,2)*pow(log(1 - z),2) - (2*CF*L*(2*TF*(-1 + pow(z,2)) + 4*CF*pow(z,2) + CA*(1 + pow(z,2)))*pow(log(z),2))/(-1 + z);
	   F3qiqbi += 4*CF*L*log(z)*(-((CA - 2*CF)*(1 + z)) + L*TF*(1 + z) + (2*(CA - 2*CF)*log(1 + z)*(1 + pow(z,2)))/(1 + z) + TF*(1 + 5*z + (8*pow(z,2))/3.)) - 
   (2*CF*L*(-36*(CA - 2*CF)*z*Li2(-z)*(1 + pow(z,2)) - 3*(CA - 2*CF)*z*(12*(-1 + pow(z,2)) + pow(Pi,2)*(1 + pow(z,2))) + 
        TF*(-1 + pow(z,2))*(3*L*(4 + 7*z + 4*pow(z,2)) + 4*(10 + z + 28*pow(z,2)))))/(9.*z*(1 + z)) + 
   2*CF*L*(-2*TF*(1 + z) - ((CA - 2*CF)*(1 + pow(z,2)))/(1 + z))*pow(log(z),2);
	}

	return DIScoeff(h,x,y,F2qiqi*fq+F2qiqbi*fqb+F2qiqj*fqj+F2qg*fg,
	                      FLqiqi*fq+FLqiqbi*fqb+FLqiqj*fqj+FLqg*fg,
	                      F3qiqi*fq+F3qiqbi*fqb+F3qiqj*fqj+F3qg*fg)*Jz;
}

}
