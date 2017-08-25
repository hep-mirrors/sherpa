#include "hard.H"
#include "Tools.H"

namespace SHNNLO {

// For DY or DY-like part of VH production

double hf0qq(double mu, double Q) {
	return 1.;
}

// hard function at nlo without the overall factor as/4/Pi
double hf1qq(double mu, double Q) {
	double L = log(mu/Q);
	return 2*cH1qq + 2*gH0qq*L - 2*CF*G0*pow(L,2);
}

// hard function at nnlo without the overall factor (as/4/Pi)^2
double hf2qq(double mu, double Q) {
	double L = log(mu/Q);
	return 2*cH2qq + 4*cH1qq*gH0qq*L + 2*gH1qq*L - 2*CF*G1*pow(L,2) + (-4*CF*cH1qq*G0 + 2*pow(gH0qq,2))*pow(L,2) - 
               4*CF*G0*gH0qq*pow(L,3) + 2*pow(CF,2)*pow(G0,2)*pow(L,4) + 
               2*beta0*(2*cH1qq*L + gH0qq*pow(L,2) - (2*CF*G0*pow(L,3))/3.);
}

//-----------------------------------------------

// For process of quark quark to diphoton process

void updatecH1qq(int i, int vm) {
#ifdef USING__VV
	// i==0 for dd and i==1 for uu
	if (vm==1) { // WW offshell
	  static WWhard ww;
	  ww.calc();
	  cH1qq = ww.cH1qq(i);
	  return;
	}
	if (vm==2) { // WW onshell
	  cH1qq = ww_os_cH1qq(i,sij(0,1),sij(0,3));
	  return;
	  //static WWhard_os ww_os; 
          //ww_os.calc();
	  //cH1qq = ww_os.cH1qq(i);
	  //return;
	}
        if (vm==3) { // WZ offshell
          static WZhard wz;
          wz.calc();
          cH1qq = wz.cH1qq(i);
          return;
        }
#endif
        // YY onsehll
	double t=sij(0,2),u=sij(0,3),s=-t-u;
	double X=log(-t/s),Y=log(-u/s);
	cH1qq = (CF*(6*u*(2*t + 3*u)*X + 6*t*(3*t + 2*u)*Y + 7*(-6 + pow(Pi,2))*(pow(t,2) + pow(u,2)) + 6*(pow(s,2) + pow(t,2))*pow(X,2) + 6*(pow(s,2) + pow(u,2))*pow(Y,2)))/
   (6.*(pow(t,2) + pow(u,2)));
}

void updatecH2qq(int i, int vm) {
#ifdef USING__VV
	// i==0 for dd and i==1 for uu
        if (vm==1) { // WW offshell
	  static WWhard ww;
          ww.calc();
          cH2qq = ww.cH2qq(i);
	  return;
        }
        if (vm==2) { // WW onshell
	  static WWhard_os ww_os;
          ww_os.calc();
	  cH2qq = ww_os.cH2qq(i);
	  return;
        }
        if (vm==3) { // WZ offshell
          static WZhard wz;
          wz.calc();
          cH2qq = wz.cH2qq(i);
          return;
        }
#endif
        // YY onshell
	double t=sij(0,2),u=sij(0,3),s=-t-u;
	double X=log(-t/s),Y=log(-u/s);
	double X2=Li2(-t/s),Y2=Li2(-u/s);
	double X3=Li3(-t/s),Y3=Li3(-u/s);
	double X4=Li4(-t/s),Y4=Li4(-u/s);
	double V=Li4(-t/u),W=Li4(-u/t);
	double EE=(288*t*u*pow(Pi,2) + 3401*pow(t,2) + 72*Zeta3*pow(t,2) - 330*pow(Pi,2)*pow(t,2) + 36*X*(-(u*(38*t + 45*u)) + 4*pow(Pi,2)*(pow(s,2) + pow(t,2))) + 
     3401*pow(u,2) + 72*Zeta3*pow(u,2) - 330*pow(Pi,2)*pow(u,2) + 36*Y*(2*t*u*(-19 + 4*pow(Pi,2)) + (-45 + 4*pow(Pi,2))*pow(t,2) + 8*pow(Pi,2)*pow(u,2)) - 
     36*(26*t*u + 38*pow(t,2) + pow(u,2))*pow(X,2) + 144*(pow(s,2) + pow(t,2))*pow(X,3) - 36*(26*t*u + pow(t,2) + 38*pow(u,2))*pow(Y,2) + 
     144*(pow(s,2) + pow(u,2))*pow(Y,3))/(162.*(pow(t,2) + pow(u,2)));
	double AA=(4*(720*u*W*pow(t,5) + 270*u*X*pow(t,5) - 120*u*X3*pow(t,5) - 720*u*X4*pow(t,5) + 720*u*X3*Y*pow(t,5) + 720*u*Y4*pow(t,5) + 120*u*Zeta3*pow(t,5) - 
       720*u*Y*Zeta3*pow(t,5) + 100*u*X*pow(Pi,2)*pow(t,5) - 120*u*X*Y*pow(Pi,2)*pow(t,5) + 28*u*pow(Pi,4)*pow(t,5) + 1440*W*pow(t,4)*pow(u,2) + 
       570*X*pow(t,4)*pow(u,2) + 600*X3*pow(t,4)*pow(u,2) - 1440*X4*pow(t,4)*pow(u,2) + 270*Y*pow(t,4)*pow(u,2) - 30*X*Y*pow(t,4)*pow(u,2) + 
       1800*X3*Y*pow(t,4)*pow(u,2) + 720*Y3*pow(t,4)*pow(u,2) + 360*X*Y3*pow(t,4)*pow(u,2) + 1440*Y4*pow(t,4)*pow(u,2) - 1320*Zeta3*pow(t,4)*pow(u,2) - 
       360*X*Zeta3*pow(t,4)*pow(u,2) - 1800*Y*Zeta3*pow(t,4)*pow(u,2) - 185*pow(Pi,2)*pow(t,4)*pow(u,2) + 220*X*pow(Pi,2)*pow(t,4)*pow(u,2) + 
       130*Y*pow(Pi,2)*pow(t,4)*pow(u,2) - 240*X*Y*pow(Pi,2)*pow(t,4)*pow(u,2) + 67*pow(Pi,4)*pow(t,4)*pow(u,2) - 
       60*t*Y2*pow(u,2)*(2*u*pow(Pi,2)*pow(s,2) + t*X*pow(u,2) + Y*(24*u*pow(t,2) + 13*pow(t,3) + 10*t*pow(u,2) - 2*pow(u,3))) + 720*V*pow(t,3)*pow(u,3) + 
       720*W*pow(t,3)*pow(u,3) + 570*X*pow(t,3)*pow(u,3) + 1440*X3*pow(t,3)*pow(u,3) + 570*Y*pow(t,3)*pow(u,3) - 60*X*Y*pow(t,3)*pow(u,3) + 
       1440*X3*Y*pow(t,3)*pow(u,3) + 1440*Y3*pow(t,3)*pow(u,3) + 1440*X*Y3*pow(t,3)*pow(u,3) - 2880*Zeta3*pow(t,3)*pow(u,3) - 
       1440*X*Zeta3*pow(t,3)*pow(u,3) - 1440*Y*Zeta3*pow(t,3)*pow(u,3) - 370*pow(Pi,2)*pow(t,3)*pow(u,3) + 240*X*pow(Pi,2)*pow(t,3)*pow(u,3) + 
       240*Y*pow(Pi,2)*pow(t,3)*pow(u,3) - 240*X*Y*pow(Pi,2)*pow(t,3)*pow(u,3) + 78*pow(Pi,4)*pow(t,3)*pow(u,3) - 
       60*u*X2*pow(t,2)*(t*(t*u*Y + 2*pow(Pi,2)*pow(s,2)) + X*(10*u*pow(t,2) - 2*pow(t,3) + 24*t*pow(u,2) + 13*pow(u,3))) + 1440*V*pow(t,2)*pow(u,4) + 
       270*X*pow(t,2)*pow(u,4) + 720*X3*pow(t,2)*pow(u,4) + 1440*X4*pow(t,2)*pow(u,4) + 570*Y*pow(t,2)*pow(u,4) - 30*X*Y*pow(t,2)*pow(u,4) + 
       360*X3*Y*pow(t,2)*pow(u,4) + 600*Y3*pow(t,2)*pow(u,4) + 1800*X*Y3*pow(t,2)*pow(u,4) - 1440*Y4*pow(t,2)*pow(u,4) - 1320*Zeta3*pow(t,2)*pow(u,4) - 
       1800*X*Zeta3*pow(t,2)*pow(u,4) - 360*Y*Zeta3*pow(t,2)*pow(u,4) - 185*pow(Pi,2)*pow(t,2)*pow(u,4) + 130*X*pow(Pi,2)*pow(t,2)*pow(u,4) + 
       220*Y*pow(Pi,2)*pow(t,2)*pow(u,4) - 240*X*Y*pow(Pi,2)*pow(t,2)*pow(u,4) + 67*pow(Pi,4)*pow(t,2)*pow(u,4) + 720*t*V*pow(u,5) + 720*t*X4*pow(u,5) + 
       270*t*Y*pow(u,5) - 120*t*Y3*pow(u,5) + 720*t*X*Y3*pow(u,5) - 720*t*Y4*pow(u,5) + 120*t*Zeta3*pow(u,5) - 720*t*X*Zeta3*pow(u,5) + 
       100*t*Y*pow(Pi,2)*pow(u,5) - 120*t*X*Y*pow(Pi,2)*pow(u,5) + 28*t*pow(Pi,4)*pow(u,5) + 180*u*pow(t,5)*pow(X,2) + 60*u*pow(Pi,2)*pow(t,5)*pow(X,2) + 
       135*pow(t,6)*pow(X,2) + 45*pow(t,4)*pow(u,2)*pow(X,2) - 360*Y*pow(t,4)*pow(u,2)*pow(X,2) + 120*pow(Pi,2)*pow(t,4)*pow(u,2)*pow(X,2) + 
       90*pow(t,3)*pow(u,3)*pow(X,2) - 720*Y*pow(t,3)*pow(u,3)*pow(X,2) + 60*pow(Pi,2)*pow(t,3)*pow(u,3)*pow(X,2) + 90*pow(t,2)*pow(u,4)*pow(X,2) - 
       420*Y*pow(t,2)*pow(u,4)*pow(X,2) - 120*u*Y*pow(t,5)*pow(X,3) - 240*Y*pow(t,4)*pow(u,2)*pow(X,3) - 120*Y*pow(t,3)*pow(u,3)*pow(X,3) + 
       30*u*pow(t,5)*pow(X,4) + 60*pow(t,4)*pow(u,2)*pow(X,4) + 30*pow(t,3)*pow(u,3)*pow(X,4) + 90*pow(t,4)*pow(u,2)*pow(Y,2) - 
       420*X*pow(t,4)*pow(u,2)*pow(Y,2) + 90*pow(t,3)*pow(u,3)*pow(Y,2) - 720*X*pow(t,3)*pow(u,3)*pow(Y,2) + 60*pow(Pi,2)*pow(t,3)*pow(u,3)*pow(Y,2) + 
       45*pow(t,2)*pow(u,4)*pow(Y,2) - 360*X*pow(t,2)*pow(u,4)*pow(Y,2) + 120*pow(Pi,2)*pow(t,2)*pow(u,4)*pow(Y,2) + 180*t*pow(u,5)*pow(Y,2) + 
       60*t*pow(Pi,2)*pow(u,5)*pow(Y,2) + 135*pow(u,6)*pow(Y,2) + 180*u*pow(t,5)*pow(X,2)*pow(Y,2) + 450*pow(t,4)*pow(u,2)*pow(X,2)*pow(Y,2) + 
       540*pow(t,3)*pow(u,3)*pow(X,2)*pow(Y,2) + 450*pow(t,2)*pow(u,4)*pow(X,2)*pow(Y,2) + 180*t*pow(u,5)*pow(X,2)*pow(Y,2) - 
       120*X*pow(t,3)*pow(u,3)*pow(Y,3) - 240*X*pow(t,2)*pow(u,4)*pow(Y,3) - 120*t*X*pow(u,5)*pow(Y,3) + 30*pow(t,3)*pow(u,3)*pow(Y,4) + 
       60*pow(t,2)*pow(u,4)*pow(Y,4) + 30*t*pow(u,5)*pow(Y,4)))/(45.*t*u*(pow(t,2) + pow(u,2))*pow(s,2));
	double DD=-(217085*u*pow(t,5) - 77760*u*W*pow(t,5) - 32400*u*X3*pow(t,5) + 12960*u*X*X3*pow(t,5) + 25920*u*X4*pow(t,5) - 92340*u*Y*pow(t,5) - 90720*u*X3*Y*pow(t,5) + 
      19440*u*Y3*pow(t,5) - 12960*u*X*Y3*pow(t,5) - 25920*u*Y*Y3*pow(t,5) - 38880*u*Y4*pow(t,5) - 80280*u*Zeta3*pow(t,5) + 116640*u*Y*Zeta3*pow(t,5) - 
      22110*u*pow(Pi,2)*pow(t,5) + 3960*u*X*pow(Pi,2)*pow(t,5) + 4680*u*Y*pow(Pi,2)*pow(t,5) + 12960*u*X*Y*pow(Pi,2)*pow(t,5) + 612*u*pow(Pi,4)*pow(t,5) + 
      434170*pow(t,4)*pow(u,2) - 155520*W*pow(t,4)*pow(u,2) - 77400*X*pow(t,4)*pow(u,2) - 168480*X3*pow(t,4)*pow(u,2) + 25920*X*X3*pow(t,4)*pow(u,2) + 
      25920*X4*pow(t,4)*pow(u,2) - 300960*Y*pow(t,4)*pow(u,2) + 64800*X*Y*pow(t,4)*pow(u,2) - 259200*X3*Y*pow(t,4)*pow(u,2) - 64800*Y3*pow(t,4)*pow(u,2) - 
      103680*X*Y3*pow(t,4)*pow(u,2) - 51840*Y*Y3*pow(t,4)*pow(u,2) - 103680*Y4*pow(t,4)*pow(u,2) - 134640*Zeta3*pow(t,4)*pow(u,2) + 
      77760*X*Zeta3*pow(t,4)*pow(u,2) + 311040*Y*Zeta3*pow(t,4)*pow(u,2) - 75900*pow(Pi,2)*pow(t,4)*pow(u,2) + 34560*X*pow(Pi,2)*pow(t,4)*pow(u,2) + 
      36000*Y*pow(Pi,2)*pow(t,4)*pow(u,2) + 43200*X*Y*pow(Pi,2)*pow(t,4)*pow(u,2) + 2304*pow(Pi,4)*pow(t,4)*pow(u,2) + 434170*pow(t,3)*pow(u,3) - 
      77760*V*pow(t,3)*pow(u,3) - 77760*W*pow(t,3)*pow(u,3) - 286020*X*pow(t,3)*pow(u,3) - 220320*X3*pow(t,3)*pow(u,3) - 12960*X*X3*pow(t,3)*pow(u,3) - 
      64800*X4*pow(t,3)*pow(u,3) - 286020*Y*pow(t,3)*pow(u,3) + 90720*X*Y*pow(t,3)*pow(u,3) - 259200*X3*Y*pow(t,3)*pow(u,3) - 220320*Y3*pow(t,3)*pow(u,3) - 
      259200*X*Y3*pow(t,3)*pow(u,3) - 12960*Y*Y3*pow(t,3)*pow(u,3) - 64800*Y4*pow(t,3)*pow(u,3) - 108720*Zeta3*pow(t,3)*pow(u,3) + 
      272160*X*Zeta3*pow(t,3)*pow(u,3) + 272160*Y*Zeta3*pow(t,3)*pow(u,3) - 88140*pow(Pi,2)*pow(t,3)*pow(u,3) + 61920*X*pow(Pi,2)*pow(t,3)*pow(u,3) + 
      61920*Y*pow(Pi,2)*pow(t,3)*pow(u,3) + 60480*X*Y*pow(Pi,2)*pow(t,3)*pow(u,3) + 3384*pow(Pi,4)*pow(t,3)*pow(u,3) + 434170*pow(t,2)*pow(u,4) - 
      155520*V*pow(t,2)*pow(u,4) - 300960*X*pow(t,2)*pow(u,4) - 64800*X3*pow(t,2)*pow(u,4) - 51840*X*X3*pow(t,2)*pow(u,4) - 103680*X4*pow(t,2)*pow(u,4) - 
      77400*Y*pow(t,2)*pow(u,4) + 64800*X*Y*pow(t,2)*pow(u,4) - 103680*X3*Y*pow(t,2)*pow(u,4) - 168480*Y3*pow(t,2)*pow(u,4) - 259200*X*Y3*pow(t,2)*pow(u,4) + 
      25920*Y*Y3*pow(t,2)*pow(u,4) + 25920*Y4*pow(t,2)*pow(u,4) - 134640*Zeta3*pow(t,2)*pow(u,4) + 311040*X*Zeta3*pow(t,2)*pow(u,4) + 
      77760*Y*Zeta3*pow(t,2)*pow(u,4) - 75900*pow(Pi,2)*pow(t,2)*pow(u,4) + 36000*X*pow(Pi,2)*pow(t,2)*pow(u,4) + 34560*Y*pow(Pi,2)*pow(t,2)*pow(u,4) + 
      43200*X*Y*pow(Pi,2)*pow(t,2)*pow(u,4) + 2304*pow(Pi,4)*pow(t,2)*pow(u,4) + 217085*t*pow(u,5) - 77760*t*V*pow(u,5) - 92340*t*X*pow(u,5) + 
      19440*t*X3*pow(u,5) - 25920*t*X*X3*pow(u,5) - 38880*t*X4*pow(u,5) - 12960*t*X3*Y*pow(u,5) - 32400*t*Y3*pow(u,5) - 90720*t*X*Y3*pow(u,5) + 
      12960*t*Y*Y3*pow(u,5) + 25920*t*Y4*pow(u,5) - 80280*t*Zeta3*pow(u,5) + 116640*t*X*Zeta3*pow(u,5) - 22110*t*pow(Pi,2)*pow(u,5) + 
      4680*t*X*pow(Pi,2)*pow(u,5) + 3960*t*Y*pow(Pi,2)*pow(u,5) + 12960*t*X*Y*pow(Pi,2)*pow(u,5) + 612*t*pow(Pi,4)*pow(u,5) - 32040*u*pow(t,5)*pow(X,2) + 
      25920*u*Y*pow(t,5)*pow(X,2) + 6480*u*pow(Pi,2)*pow(t,5)*pow(X,2) + 9720*pow(t,6)*pow(X,2) - 156600*pow(t,4)*pow(u,2)*pow(X,2) + 
      90720*Y*pow(t,4)*pow(u,2)*pow(X,2) + 23760*pow(Pi,2)*pow(t,4)*pow(u,2)*pow(X,2) - 165060*pow(t,3)*pow(u,3)*pow(X,2) + 
      113400*Y*pow(t,3)*pow(u,3)*pow(X,2) + 30240*pow(Pi,2)*pow(t,3)*pow(u,3)*pow(X,2) - 75960*pow(t,2)*pow(u,4)*pow(X,2) + 
      58320*Y*pow(t,2)*pow(u,4)*pow(X,2) + 15120*pow(Pi,2)*pow(t,2)*pow(u,4)*pow(X,2) - 6300*t*pow(u,5)*pow(X,2) + 9720*t*Y*pow(u,5)*pow(X,2) + 
      2160*t*pow(Pi,2)*pow(u,5)*pow(X,2) + 12600*u*pow(t,5)*pow(X,3) + 12960*u*Y*pow(t,5)*pow(X,3) + 45360*pow(t,4)*pow(u,2)*pow(X,3) + 
      25920*Y*pow(t,4)*pow(u,2)*pow(X,3) + 60840*pow(t,3)*pow(u,3)*pow(X,3) + 6480*Y*pow(t,3)*pow(u,3)*pow(X,3) + 36000*pow(t,2)*pow(u,4)*pow(X,3) - 
      12960*Y*pow(t,2)*pow(u,4)*pow(X,3) + 7920*t*pow(u,5)*pow(X,3) - 6480*t*Y*pow(u,5)*pow(X,3) - 540*u*pow(t,5)*pow(X,4) + 
      1080*pow(t,4)*pow(u,2)*pow(X,4) + 3780*pow(t,3)*pow(u,3)*pow(X,4) + 2160*pow(t,2)*pow(u,4)*pow(X,4) - 6300*u*pow(t,5)*pow(Y,2) + 
      9720*u*X*pow(t,5)*pow(Y,2) + 2160*u*pow(Pi,2)*pow(t,5)*pow(Y,2) - 75960*pow(t,4)*pow(u,2)*pow(Y,2) + 58320*X*pow(t,4)*pow(u,2)*pow(Y,2) + 
      15120*pow(Pi,2)*pow(t,4)*pow(u,2)*pow(Y,2) - 165060*pow(t,3)*pow(u,3)*pow(Y,2) + 113400*X*pow(t,3)*pow(u,3)*pow(Y,2) + 
      30240*pow(Pi,2)*pow(t,3)*pow(u,3)*pow(Y,2) - 156600*pow(t,2)*pow(u,4)*pow(Y,2) + 90720*X*pow(t,2)*pow(u,4)*pow(Y,2) + 
      23760*pow(Pi,2)*pow(t,2)*pow(u,4)*pow(Y,2) - 32040*t*pow(u,5)*pow(Y,2) + 25920*t*X*pow(u,5)*pow(Y,2) + 6480*t*pow(Pi,2)*pow(u,5)*pow(Y,2) + 
      9720*pow(u,6)*pow(Y,2) - 22680*u*pow(t,5)*pow(X,2)*pow(Y,2) - 64800*pow(t,4)*pow(u,2)*pow(X,2)*pow(Y,2) - 84240*pow(t,3)*pow(u,3)*pow(X,2)*pow(Y,2) - 
      64800*pow(t,2)*pow(u,4)*pow(X,2)*pow(Y,2) - 22680*t*pow(u,5)*pow(X,2)*pow(Y,2) + 
      2160*t*Y2*pow(s,2)*pow(u,2)*(9*u*X + 3*(16*t + 5*u)*Y - 4*u*pow(Pi,2) - 3*u*pow(X,2) - 6*s*pow(Y,2)) - 
      2160*u*X2*pow(s,2)*pow(t,2)*(-3*(5*t + 16*u)*X + 6*s*pow(X,2) + t*(-9*Y + 4*pow(Pi,2) + 3*pow(Y,2))) + 7920*u*pow(t,5)*pow(Y,3) - 
      6480*u*X*pow(t,5)*pow(Y,3) + 36000*pow(t,4)*pow(u,2)*pow(Y,3) - 12960*X*pow(t,4)*pow(u,2)*pow(Y,3) + 60840*pow(t,3)*pow(u,3)*pow(Y,3) + 
      6480*X*pow(t,3)*pow(u,3)*pow(Y,3) + 45360*pow(t,2)*pow(u,4)*pow(Y,3) + 25920*X*pow(t,2)*pow(u,4)*pow(Y,3) + 12600*t*pow(u,5)*pow(Y,3) + 
      12960*t*X*pow(u,5)*pow(Y,3) + 2160*pow(t,4)*pow(u,2)*pow(Y,4) + 3780*pow(t,3)*pow(u,3)*pow(Y,4) + 1080*pow(t,2)*pow(u,4)*pow(Y,4) - 
      540*t*pow(u,5)*pow(Y,4))/(3240.*t*u*(pow(t,2) + pow(u,2))*pow(s,2));
	double BB=(1440*u*X*pow(Pi,2)*pow(t,7) + 17595*pow(t,6)*pow(u,2) - 20160*W*pow(t,6)*pow(u,2) - 1440*X*pow(t,6)*pow(u,2) - 11520*X3*pow(t,6)*pow(u,2) + 
     8640*X*X3*pow(t,6)*pow(u,2) - 14580*Y*pow(t,6)*pow(u,2) - 23040*X3*Y*pow(t,6)*pow(u,2) + 2160*Y3*pow(t,6)*pow(u,2) - 2880*X*Y3*pow(t,6)*pow(u,2) - 
     2880*Y*Y3*pow(t,6)*pow(u,2) - 15840*Y4*pow(t,6)*pow(u,2) - 3600*Zeta3*pow(t,6)*pow(u,2) - 5760*X*Zeta3*pow(t,6)*pow(u,2) + 
     25920*Y*Zeta3*pow(t,6)*pow(u,2) - 4440*pow(Pi,2)*pow(t,6)*pow(u,2) + 5280*X*pow(Pi,2)*pow(t,6)*pow(u,2) + 1440*Y*pow(Pi,2)*pow(t,6)*pow(u,2) + 
     6240*X*Y*pow(Pi,2)*pow(t,6)*pow(u,2) + 506*pow(Pi,4)*pow(t,6)*pow(u,2) + 33750*pow(t,5)*pow(u,3) - 40320*W*pow(t,5)*pow(u,3) - 
     14040*X*pow(t,5)*pow(u,3) - 50400*X3*pow(t,5)*pow(u,3) + 23040*X*X3*pow(t,5)*pow(u,3) - 14400*X4*pow(t,5)*pow(u,3) - 48960*Y*pow(t,5)*pow(u,3) + 
     14400*X*Y*pow(t,5)*pow(u,3) - 69120*X3*Y*pow(t,5)*pow(u,3) - 23040*Y3*pow(t,5)*pow(u,3) - 28800*X*Y3*pow(t,5)*pow(u,3) - 46080*Y4*pow(t,5)*pow(u,3) - 
     1440*Zeta3*pow(t,5)*pow(u,3) + 5760*X*Zeta3*pow(t,5)*pow(u,3) + 69120*Y*Zeta3*pow(t,5)*pow(u,3) - 17040*pow(Pi,2)*pow(t,5)*pow(u,3) + 
     13200*X*pow(Pi,2)*pow(t,5)*pow(u,3) + 9840*Y*pow(Pi,2)*pow(t,5)*pow(u,3) + 22080*X*Y*pow(Pi,2)*pow(t,5)*pow(u,3) + 1268*pow(Pi,4)*pow(t,5)*pow(u,3) + 
     32310*pow(t,4)*pow(u,4) - 20160*V*pow(t,4)*pow(u,4) - 20160*W*pow(t,4)*pow(u,4) - 46980*X*pow(t,4)*pow(u,4) - 64080*X3*pow(t,4)*pow(u,4) + 
     17280*X*X3*pow(t,4)*pow(u,4) - 44640*X4*pow(t,4)*pow(u,4) - 46980*Y*pow(t,4)*pow(u,4) + 20160*X*Y*pow(t,4)*pow(u,4) - 72000*X3*Y*pow(t,4)*pow(u,4) - 
     64080*Y3*pow(t,4)*pow(u,4) - 72000*X*Y3*pow(t,4)*pow(u,4) + 17280*Y*Y3*pow(t,4)*pow(u,4) - 44640*Y4*pow(t,4)*pow(u,4) + 4320*Zeta3*pow(t,4)*pow(u,4) + 
     54720*X*Zeta3*pow(t,4)*pow(u,4) + 54720*Y*Zeta3*pow(t,4)*pow(u,4) - 20880*pow(Pi,2)*pow(t,4)*pow(u,4) + 17760*X*pow(Pi,2)*pow(t,4)*pow(u,4) + 
     17760*Y*pow(Pi,2)*pow(t,4)*pow(u,4) + 31680*X*Y*pow(Pi,2)*pow(t,4)*pow(u,4) + 1524*pow(Pi,4)*pow(t,4)*pow(u,4) + 33750*pow(t,3)*pow(u,5) - 
     40320*V*pow(t,3)*pow(u,5) - 48960*X*pow(t,3)*pow(u,5) - 23040*X3*pow(t,3)*pow(u,5) - 46080*X4*pow(t,3)*pow(u,5) - 14040*Y*pow(t,3)*pow(u,5) + 
     14400*X*Y*pow(t,3)*pow(u,5) - 28800*X3*Y*pow(t,3)*pow(u,5) - 50400*Y3*pow(t,3)*pow(u,5) - 69120*X*Y3*pow(t,3)*pow(u,5) + 23040*Y*Y3*pow(t,3)*pow(u,5) - 
     14400*Y4*pow(t,3)*pow(u,5) - 1440*Zeta3*pow(t,3)*pow(u,5) + 69120*X*Zeta3*pow(t,3)*pow(u,5) + 5760*Y*Zeta3*pow(t,3)*pow(u,5) - 
     17040*pow(Pi,2)*pow(t,3)*pow(u,5) + 9840*X*pow(Pi,2)*pow(t,3)*pow(u,5) + 13200*Y*pow(Pi,2)*pow(t,3)*pow(u,5) + 22080*X*Y*pow(Pi,2)*pow(t,3)*pow(u,5) + 
     1268*pow(Pi,4)*pow(t,3)*pow(u,5) + 17595*pow(t,2)*pow(u,6) - 20160*V*pow(t,2)*pow(u,6) - 14580*X*pow(t,2)*pow(u,6) + 2160*X3*pow(t,2)*pow(u,6) - 
     2880*X*X3*pow(t,2)*pow(u,6) - 15840*X4*pow(t,2)*pow(u,6) - 1440*Y*pow(t,2)*pow(u,6) - 2880*X3*Y*pow(t,2)*pow(u,6) - 11520*Y3*pow(t,2)*pow(u,6) - 
     23040*X*Y3*pow(t,2)*pow(u,6) + 8640*Y*Y3*pow(t,2)*pow(u,6) - 3600*Zeta3*pow(t,2)*pow(u,6) + 25920*X*Zeta3*pow(t,2)*pow(u,6) - 
     5760*Y*Zeta3*pow(t,2)*pow(u,6) - 4440*pow(Pi,2)*pow(t,2)*pow(u,6) + 1440*X*pow(Pi,2)*pow(t,2)*pow(u,6) + 5280*Y*pow(Pi,2)*pow(t,2)*pow(u,6) + 
     6240*X*Y*pow(Pi,2)*pow(t,2)*pow(u,6) + 506*pow(Pi,4)*pow(t,2)*pow(u,6) + 1440*t*Y*pow(Pi,2)*pow(u,7) + 1440*u*pow(t,7)*pow(X,2) + 
     4320*u*pow(Pi,2)*pow(t,7)*pow(X,2) + 720*pow(Pi,2)*pow(t,8)*pow(X,2) - 1080*pow(t,6)*pow(u,2)*pow(X,2) + 9000*Y*pow(t,6)*pow(u,2)*pow(X,2) + 
     12480*pow(Pi,2)*pow(t,6)*pow(u,2)*pow(X,2) - 15840*pow(t,5)*pow(u,3)*pow(X,2) + 29520*Y*pow(t,5)*pow(u,3)*pow(X,2) + 
     19920*pow(Pi,2)*pow(t,5)*pow(u,3)*pow(X,2) - 20160*pow(t,4)*pow(u,4)*pow(X,2) + 33120*Y*pow(t,4)*pow(u,4)*pow(X,2) + 
     17040*pow(Pi,2)*pow(t,4)*pow(u,4)*pow(X,2) - 12960*pow(t,3)*pow(u,5)*pow(X,2) + 13680*Y*pow(t,3)*pow(u,5)*pow(X,2) + 
     6960*pow(Pi,2)*pow(t,3)*pow(u,5)*pow(X,2) - 1800*pow(t,2)*pow(u,6)*pow(X,2) + 1080*Y*pow(t,2)*pow(u,6)*pow(X,2) + 
     960*pow(Pi,2)*pow(t,2)*pow(u,6)*pow(X,2) + 720*u*pow(t,7)*pow(X,3) + 3960*pow(t,6)*pow(u,2)*pow(X,3) + 10080*Y*pow(t,6)*pow(u,2)*pow(X,3) + 
     11280*pow(t,5)*pow(u,3)*pow(X,3) + 20160*Y*pow(t,5)*pow(u,3)*pow(X,3) + 15720*pow(t,4)*pow(u,4)*pow(X,3) + 9360*Y*pow(t,4)*pow(u,4)*pow(X,3) + 
     9840*pow(t,3)*pow(u,5)*pow(X,3) - 1440*Y*pow(t,3)*pow(u,5)*pow(X,3) + 2160*pow(t,2)*pow(u,6)*pow(X,3) - 720*Y*pow(t,2)*pow(u,6)*pow(X,3) + 
     1080*u*pow(t,7)*pow(X,4) + 180*pow(t,8)*pow(X,4) + 3000*pow(t,6)*pow(u,2)*pow(X,4) + 5040*pow(t,5)*pow(u,3)*pow(X,4) + 4860*pow(t,4)*pow(u,4)*pow(X,4) + 
     2280*pow(t,3)*pow(u,5)*pow(X,4) + 360*pow(t,2)*pow(u,6)*pow(X,4) - 1800*pow(t,6)*pow(u,2)*pow(Y,2) + 1080*X*pow(t,6)*pow(u,2)*pow(Y,2) + 
     960*pow(Pi,2)*pow(t,6)*pow(u,2)*pow(Y,2) - 12960*pow(t,5)*pow(u,3)*pow(Y,2) + 13680*X*pow(t,5)*pow(u,3)*pow(Y,2) + 
     6960*pow(Pi,2)*pow(t,5)*pow(u,3)*pow(Y,2) - 20160*pow(t,4)*pow(u,4)*pow(Y,2) + 33120*X*pow(t,4)*pow(u,4)*pow(Y,2) + 
     17040*pow(Pi,2)*pow(t,4)*pow(u,4)*pow(Y,2) - 15840*pow(t,3)*pow(u,5)*pow(Y,2) + 29520*X*pow(t,3)*pow(u,5)*pow(Y,2) + 
     19920*pow(Pi,2)*pow(t,3)*pow(u,5)*pow(Y,2) - 1080*pow(t,2)*pow(u,6)*pow(Y,2) + 9000*X*pow(t,2)*pow(u,6)*pow(Y,2) + 
     12480*pow(Pi,2)*pow(t,2)*pow(u,6)*pow(Y,2) + 1440*t*pow(u,7)*pow(Y,2) + 4320*t*pow(Pi,2)*pow(u,7)*pow(Y,2) + 720*pow(Pi,2)*pow(u,8)*pow(Y,2) - 
     5040*pow(t,6)*pow(u,2)*pow(X,2)*pow(Y,2) - 14400*pow(t,5)*pow(u,3)*pow(X,2)*pow(Y,2) - 18720*pow(t,4)*pow(u,4)*pow(X,2)*pow(Y,2) - 
     14400*pow(t,3)*pow(u,5)*pow(X,2)*pow(Y,2) - 5040*pow(t,2)*pow(u,6)*pow(X,2)*pow(Y,2) - 
     720*X2*pow(s,2)*pow(t,3)*pow(u,2)*(-2*(8*t + 19*u)*X + 2*s*pow(X,2) + t*(-3*Y + 2*pow(Pi,2) + pow(Y,2))) + 
     720*Y2*pow(s,2)*pow(t,2)*pow(u,3)*(3*u*X - u*pow(X,2) + 2*((19*t + 8*u)*Y - u*pow(Pi,2) - s*pow(Y,2))) + 2160*pow(t,6)*pow(u,2)*pow(Y,3) - 
     720*X*pow(t,6)*pow(u,2)*pow(Y,3) + 9840*pow(t,5)*pow(u,3)*pow(Y,3) - 1440*X*pow(t,5)*pow(u,3)*pow(Y,3) + 15720*pow(t,4)*pow(u,4)*pow(Y,3) + 
     9360*X*pow(t,4)*pow(u,4)*pow(Y,3) + 11280*pow(t,3)*pow(u,5)*pow(Y,3) + 20160*X*pow(t,3)*pow(u,5)*pow(Y,3) + 3960*pow(t,2)*pow(u,6)*pow(Y,3) + 
     10080*X*pow(t,2)*pow(u,6)*pow(Y,3) + 720*t*pow(u,7)*pow(Y,3) + 360*pow(t,6)*pow(u,2)*pow(Y,4) + 2280*pow(t,5)*pow(u,3)*pow(Y,4) + 
     4860*pow(t,4)*pow(u,4)*pow(Y,4) + 5040*pow(t,3)*pow(u,5)*pow(Y,4) + 3000*pow(t,2)*pow(u,6)*pow(Y,4) + 1080*t*pow(u,7)*pow(Y,4) + 180*pow(u,8)*pow(Y,4))/
   (360.*(pow(t,2) + pow(u,2))*pow(s,2)*pow(t,2)*pow(u,2));
	cH2qq = EE*CF*nf*TF + AA*CF*sumQsq/(i==0?1./9.:4./9.)*TF + DD*CA*CF + BB*pow(CF,2);
}


//===============================================

// Higgs production
// the corresponding 5-flavor QCD (i.e. Higgs effective theory) calculation should use hard function:
//	hf0tt(mu,Q,mt) * | 1 + as/4/Pi * hf1tt(mu,Q,mt) + (as/4/Pi)^2 * hf2tt(mu,Q,mt) |^2
// which is from integrating out top quark, i.e. matching from full QCD to 5-flavor QCD
// in heavy top limit, call these functions with Q=0, i.e. hf0tt(mu,0,mt),hf1tt(mu,0,mt),hf2tt(mu,0,mt)

double hf0tt(double mu, double Q, double mt) {
	double z = pow(Q/2./mt,2);
	if (z==0) return 1;
	double re = 1 - ( z<1 ? (1.-z)/z*pow(asin(sqrt(z)),2) : (z-1.)/z*(pow(log(sqrt(z)+sqrt(z-1)),2)-pow(Pi,2)/4.) ) ;
	double im = ( z<1 ? 0 : (z-1.)/z*(log(sqrt(z)+sqrt(z-1))*Pi) ) ;
	return pow(3./2./z,2)*(pow(re,2) + pow(im,2));
}

// hard function at nlo without the overall factor as/4/Pi
double hf1tt(double mu, double Q, double mt) {
        double z = pow(Q/2./mt,2);
	return (5.-z*38./45.-pow(z,2)*1289./4725.-pow(z,3)*155./1134.-pow(z,4)*5385047./65488500.)*CA +
               (-3.+z*307./9.+pow(z,2)*25813./18900.+pow(z,3)*3055907./3969000.+pow(z,4)*659504801./1309770000.)*CF;
}

// hard function at nnlo without the overall factor (as/4/Pi)^2
double hf2tt(double mu, double Q, double mt) {
	double Lmt = 2.*log(mt/mu);
	//double z = pow(Q/2./mt,2);
	return - (-beta1 + beta0*hf1tt(mu,Q,mt))*Lmt
               - (100.*CA*CF)/3.- (5.*CA*TF)/6. - (4.*CF*TF)/3. - (47.*CA*nf*TF)/9. - 5.*CF*nf*TF + (1063.*pow(CA,2))/36. + (27.*pow(CF,2))/2.;
}

// hard function at nnnlo without overall factor (as/4/pi)^3
double hf3tt(double mu, double Q, double mt) {
	double Lmt = 2.*log(mt/mu);
	//double z = pow(Q/2./mt,2);
	return 64*(-2892659/41472. - (1733*Lmt)/288. + (897943*Zeta3)/9216. + nf*(40291/20736. - (55*Lmt)/54. - (110779*Zeta3)/13824. + (23*pow(Lmt,2))/32.) + 
     (209*pow(Lmt,2))/64. + (-6865/31104. - (77*Lmt)/1728. - pow(Lmt,2)/18.)*pow(nf,2));
}

//-----------------------------------------------

// For process of gluon gluon to color neutral particle such as Higgs

double hf0gg(double mu, double Q) {
	return 1.;
}

// hard function from matching 5-flavor QCD to SCET
double hf1gg(double mu, double Q) {
	double L = log(mu/Q);
	return 2*(cH1gg + (2*beta0 + gH0gg)*L - CA*G0*pow(L,2));
}

// hard function from matching 5-flavor QCD to SCET
double hf2gg(double mu, double Q) {
	double L = log(mu/Q);
	return 12*beta0*cH1gg*L + 4*cH1gg*gH0gg*L - 4*CA*cH1gg*G0*pow(L,2) + 10*beta0*gH0gg*pow(L,2) + 12*pow(beta0,2)*pow(L,2) + 
               2*pow(gH0gg,2)*pow(L,2) + 2*(cH2gg + (4*beta1 + gH1gg)*L - CA*G1*pow(L,2)) - (28*beta0*CA*G0*pow(L,3))/3. - 
               4*CA*G0*gH0gg*pow(L,3) + 2*pow(CA,2)*pow(G0,2)*pow(L,4);
}

}
