#include "coeffqt.H"
#include "Tools.H"

namespace SHNNLO {

//-----------------------------------------------

double S2(double z) {
	return -pow(Pi,2)/6. - 2*Li2(-z) - 2*log(z)*log(1 + z);
}

//-----------------------------------------------

double Pgq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (f*(2 - 2*z + pow(z,2)))/pow(z,2);
}
double Pqg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return f*(-2 + 1/z + 2*z);
}
double Pqq(double z, double x, double f, double f0) {
	if (z<x) return 2*f0*log(1 - x);
	return 2*f0*log(1 - x) - (-2*f0*z + f*(1 + pow(z,2)))/((-1 + z)*z);
}
double Pgg(double z, double x, double f, double f0) {
	if (z<x) return 2*f0*log(1 - x);
	return 2*f0*log(1 - x) - (2*(-(f0*pow(z,2)) + f*pow(1 - z + pow(z,2),2)))/((-1 + z)*pow(z,2));
}

double PPqqV(double z, double x, double f, double f0) {
	if (z<x) return (f0*(3*beta0 + 9*CF + CA*(6 - 72*Zeta3) + 144*CF*Zeta3 + 6*G1*log(1 - x) + 4*beta0*pow(Pi,2) - 12*CF*pow(Pi,2)))/24.;
	return (f0*(3*beta0 + 9*CF + CA*(6 - 72*Zeta3) + 144*CF*Zeta3 + 6*G1*log(1 - x) + 4*beta0*pow(Pi,2) - 12*CF*pow(Pi,2)))/24. - 
   (-2*f0*G1*z + f*(8*beta0 - 40*CF + G1 - 16*beta0*z + 80*CF*z + 24*CA*pow(-1 + z,2) + 8*beta0*pow(z,2) - 40*CF*pow(z,2) + 
         G1*pow(z,2) - 4*log(z)*(-beta0 + 6*CF + 4*CF*z + 2*CA*(-1 + pow(z,2)) - beta0*pow(z,2) - 4*CF*pow(z,2) + 
            4*CF*log(1 - z)*(1 + pow(z,2))) + 4*(CF*(-1 + pow(z,2)) + CA*(1 + pow(z,2)))*pow(log(z),2)))/(8.*(-1 + z)*z);
}

double PPqqbV(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return -((CA - 2*CF)*f*(4 - 4*pow(z,2) + 2*log(z)*pow(1 + z,2) + (1 + pow(z,2))*pow(log(z),2) + 2*(1 + pow(z,2))*S2(z)))/(2.*z*(1 + z));
}

double PPqqS(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (f*TF*(20 - 18*z + 54*pow(z,2) + 3*z*log(z)*(3 + 15*z + 8*pow(z,2)) - 56*pow(z,3) - 9*z*(1 + z)*pow(log(z),2)))/(9.*pow(z,2));
}

double PPqg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (f*(40*CA - 36*CA*z + 126*CF*z + 18*CA*z*log(z) + 27*CF*z*log(z) + 3*CA*z*pow(Pi,2) - 6*CF*z*pow(Pi,2) + 450*CA*pow(z,2) - 
       261*CF*pow(z,2) + 144*CA*log(z)*pow(z,2) - 36*CF*log(z)*pow(z,2) - 6*CA*pow(Pi,2)*pow(z,2) + 12*CF*pow(Pi,2)*pow(z,2) + 
       36*z*log(1 - z)*(2*(CA - CF)*(-1 + z)*z - CF*log(z)*(1 - 2*z + 2*pow(z,2))) - 436*CA*pow(z,3) + 180*CF*pow(z,3) + 
       264*CA*log(z)*pow(z,3) + 72*CF*log(z)*pow(z,3) + 6*CA*pow(Pi,2)*pow(z,3) - 12*CF*pow(Pi,2)*pow(z,3) - 
       18*(CA - CF)*z*(1 - 2*z + 2*pow(z,2))*pow(log(1 - z),2) - 18*CA*z*pow(log(z),2) + 9*CF*z*pow(log(z),2) - 
       36*CA*pow(z,2)*pow(log(z),2) - 18*CF*pow(z,2)*pow(log(z),2) + 36*CF*pow(z,3)*pow(log(z),2) + 
       18*CA*z*(1 + 2*z + 2*pow(z,2))*S2(z)))/(18.*pow(z,2));
}

double PPgg(double z, double x, double f, double f0) {
	if (z<x) return f0*(beta0 - (CF*nf*TF)/CA + CA*(-1 + 3*Zeta3) + (G1*log(1 - x))/4.);
	return f0*(beta0 - (CF*nf*TF)/CA + CA*(-1 + 3*Zeta3) + (G1*log(1 - x))/4.) - 
   (-9*CA*f0*G1*(1 + z)*pow(z,2) + f*(72*z*pow(log(z),2)*(CF*nf*TF*(-1 + z)*pow(1 + z,2) + pow(CA,2)*pow(1 + z - pow(z,2),2)) - 
         12*(1 + z)*log(z)*(-((-1 + z)*z*(-3*beta0*CA*(1 + z) + 6*CF*nf*TF*(3 + 5*z) + 4*pow(CA,2)*(9 + 11*pow(z,2)))) + 
            12*log(1 - z)*pow(CA,2)*pow(1 - z + pow(z,2),2)) - 
         (1 + z)*(48*CF*nf*TF*pow(-1 + z,2)*(-1 + 11*z + 5*pow(z,2)) + 2*pow(CA,2)*pow(-1 + z,2)*(277 - 65*z + 277*pow(z,2)) - 
            3*CA*(2*beta0*pow(-1 + z,2)*(13 + 4*z + 13*pow(z,2)) + 3*G1*pow(1 - z + pow(z,2),2))) + 
         72*(-1 + z)*pow(CA,2)*pow(1 + z + pow(z,2),2)*S2(z)))/(36.*CA*(-1 + z)*(1 + z)*pow(z,2));
}

double PPgq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return -(f*(-60*beta0 + 202*CA + 60*beta0*z - 258*CA*z + 45*CF*z + 216*CA*z*log(z) - 36*CF*z*log(z) + 6*CA*pow(Pi,2) - 6*CA*z*pow(Pi,2) - 
        48*beta0*pow(z,2) + 102*CA*pow(z,2) + 63*CF*pow(z,2) + 90*CA*log(z)*pow(z,2) - 63*CF*log(z)*pow(z,2) + 
        3*CA*pow(Pi,2)*pow(z,2) + 18*log(1 - z)*(-2*beta0 + 2*beta0*z - beta0*pow(z,2) - 2*CA*pow(z,2) + 
           2*CA*log(z)*(2 - 2*z + pow(z,2)) + CF*(6 - 6*z + 5*pow(z,2))) - 88*CA*pow(z,3) + 48*CA*log(z)*pow(z,3) - 
        18*(CA - CF)*(2 - 2*z + pow(z,2))*pow(log(1 - z),2) - 36*CA*z*pow(log(z),2) + 18*CF*z*pow(log(z),2) - 
        18*CA*pow(z,2)*pow(log(z),2) - 9*CF*pow(z,2)*pow(log(z),2) + 18*CA*(2 + 2*z + pow(z,2))*S2(z)))/(18.*pow(z,2));
}

//-----------------------------------------------

double P0qiqi(double z, double x, double f, double f0) { return CF*Pqq(z,x,f,f0)+CF*3/2.*f0; }
double P0qg(double z, double x, double f, double f0) { return TF*Pqg(z,x,f,f0); }
double P0gg(double z, double x, double f, double f0) { return CA*Pgg(z,x,f,f0)+beta0/2.*f0; }
double P0gq(double z, double x, double f, double f0) { return CF*Pgq(z,x,f,f0); }

double P1qiqi(double z, double x, double f, double f0) { return CF*(PPqqV(z,x,f,f0) + PPqqS(z,x,f,f0)); }
double P1qiqbi(double z, double x, double f, double f0) { return CF*(PPqqbV(z,x,f,f0) + PPqqS(z,x,f,f0)); }
double P1qiqj(double z, double x, double f, double f0) { return CF*PPqqS(z,x,f,f0); }
double P1qg(double z, double x, double f, double f0) { return TF*PPqg(z,x,f,f0); }
double P1gg(double z, double x, double f, double f0) { return CA*PPgg(z,x,f,f0); }
double P1gq(double z, double x, double f, double f0) { return CF*PPgq(z,x,f,f0); }

double I1qiqi(double z, double x, double f, double f0) {
	if (z<x) return -(CF*f0*pow(Pi,2))/6.;
	return (-2*CF*f*(-1 + z))/z - (CF*f0*pow(Pi,2))/6.;
}
double I1qg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return 4*f*TF*(1 - z);
}
double I1gg(double z, double x, double f, double f0) {
	return -(CA*f0*pow(Pi,2))/6.;
}
double I1gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return 2*CF*f;
}
double Ii1gg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (-4*CA*f*(1 - z))/pow(z,2);
}
double Ii1gq(double z, double x, double f, double f0) { 
	if (z<x) return 0;
	return (-4*CF*f*(1 - z))/pow(z,2);
}

double I2qiqi(double z, double x, double f, double f0) {
	if (z<x) return (CF*f0*(-2624*nf*TF + 2016*nf*TF*Zeta3 + 96*(56*nf*TF + CA*(-202 + 189*Zeta3))*log(1 - x) + 360*nf*TF*pow(Pi,2) + 
       9*CF*pow(Pi,4) + 2*CA*(4856 - 2772*Zeta3 - 603*pow(Pi,2) + 18*pow(Pi,4))))/648.;
	return (CF*f0*(-2624*nf*TF + 2016*nf*TF*Zeta3 + 96*(56*nf*TF + CA*(-202 + 189*Zeta3))*log(1 - x) + 360*nf*TF*pow(Pi,2) + 
        9*CF*pow(Pi,4) + 2*CA*(4856 - 2772*Zeta3 - 603*pow(Pi,2) + 18*pow(Pi,4))))/648. - 
   (CF*(-8*f0*(56*nf*TF + CA*(-202 + 189*Zeta3))*pow(z,2) + 
        f*(688*TF - 16*CA*z - 1188*CF*z - 1948*TF*z + 152*nf*TF*z + 324*CA*z*Zeta3 + 1080*CF*z*Zeta3 - 
           216*CA*z*Li3(1 - z) + 216*CF*z*Li3(1 - z) + 432*CA*z*Li3(z) - 1080*CF*z*Li3(z) - 348*CA*z*log(z) + 
           540*CF*z*log(z) + 504*TF*z*log(z) + 120*nf*TF*z*log(z) + 288*TF*log(1 - z)*log(z) + 
           216*CA*z*log(1 - z)*log(z) - 648*CF*z*log(1 - z)*log(z) - 720*TF*z*log(1 - z)*log(z) - 48*TF*pow(Pi,2) - 
           54*CA*z*pow(Pi,2) + 54*CF*z*pow(Pi,2) + 120*TF*z*pow(Pi,2) - 1584*CA*pow(z,2) + 2376*CF*pow(z,2) + 
           2376*TF*pow(z,2) + 144*nf*TF*pow(z,2) + 108*CA*log(1 - z)*pow(z,2) - 108*CF*log(1 - z)*pow(z,2) + 
           432*CA*log(z)*pow(z,2) - 1404*CF*log(z)*pow(z,2) - 1224*TF*log(z)*pow(z,2) - 
           432*CA*log(1 - z)*log(z)*pow(z,2) + 1296*CF*log(1 - z)*log(z)*pow(z,2) + 
           864*TF*log(1 - z)*log(z)*pow(z,2) + 108*CA*pow(Pi,2)*pow(z,2) - 108*CF*pow(Pi,2)*pow(z,2) - 
           144*TF*pow(Pi,2)*pow(z,2) + 216*(CA - CF)*z*Li2(1 - z)*log(1 - z)*(1 + pow(z,2)) + 
           72*Li2(z)*(-3*(CA - 3*CF)*z*log(z)*(1 + pow(z,2)) + 
              pow(-1 + z,2)*(3*(CA - 2*CF)*z + TF*(4 - 2*z + 4*pow(z,2)))) - 16*CA*pow(z,3) - 1188*CF*pow(z,3) - 
           1660*TF*pow(z,3) + 152*nf*TF*pow(z,3) + 324*CA*Zeta3*pow(z,3) + 1080*CF*Zeta3*pow(z,3) - 
           216*CA*Li3(1 - z)*pow(z,3) + 216*CF*Li3(1 - z)*pow(z,3) + 432*CA*Li3(z)*pow(z,3) - 
           1080*CF*Li3(z)*pow(z,3) - 108*CA*log(1 - z)*pow(z,3) + 108*CF*log(1 - z)*pow(z,3) - 
           996*CA*log(z)*pow(z,3) + 1728*CF*log(z)*pow(z,3) + 1488*TF*log(z)*pow(z,3) + 120*nf*TF*log(z)*pow(z,3) + 
           216*CA*log(1 - z)*log(z)*pow(z,3) - 648*CF*log(1 - z)*log(z)*pow(z,3) - 720*TF*log(1 - z)*log(z)*pow(z,3) - 
           54*CA*pow(Pi,2)*pow(z,3) + 54*CF*pow(Pi,2)*pow(z,3) + 120*TF*pow(Pi,2)*pow(z,3) + 544*TF*pow(z,4) - 
           768*TF*log(z)*pow(z,4) + 288*TF*log(1 - z)*log(z)*pow(z,4) - 48*TF*pow(Pi,2)*pow(z,4) + 
           108*CF*z*log(z)*pow(log(1 - z),2) + 108*CF*log(z)*pow(z,3)*pow(log(1 - z),2) - 99*CA*z*pow(log(z),2) + 
           162*CF*z*pow(log(z),2) - 54*TF*z*pow(log(z),2) + 36*nf*TF*z*pow(log(z),2) + 
           108*CF*z*log(1 - z)*pow(log(z),2) - 108*CA*pow(z,2)*pow(log(z),2) + 108*CF*pow(z,2)*pow(log(z),2) + 
           9*CA*pow(z,3)*pow(log(z),2) - 108*CF*pow(z,3)*pow(log(z),2) - 90*TF*pow(z,3)*pow(log(z),2) + 
           36*nf*TF*pow(z,3)*pow(log(z),2) + 108*CF*log(1 - z)*pow(z,3)*pow(log(z),2) + 
           144*TF*pow(z,4)*pow(log(z),2) - 18*CA*z*pow(log(z),3) + 18*CF*z*pow(log(z),3) + 36*TF*z*pow(log(z),3) - 
           18*CA*pow(z,3)*pow(log(z),3) - 18*CF*pow(z,3)*pow(log(z),3) - 36*TF*pow(z,3)*pow(log(z),3))))/(54.*(-1 + z)*pow(z,2));
}

double I2qiqbi(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*(344*TF - 405*CA*z + 810*CF*z - 286*TF*z + 162*CA*z*Zeta3 - 324*CF*z*Zeta3 - 324*CA*z*Li3(-z) + 648*CF*z*Li3(-z) - 
       216*CA*z*Li3(z) + 432*CF*z*Li3(z) - 216*CA*z*Li3(1/(1 + z)) + 432*CF*z*Li3(1/(1 + z)) - 81*CA*z*log(z) + 
       162*CF*z*log(z) + 252*TF*z*log(z) + 144*TF*log(1 - z)*log(z) - 108*CA*z*log(1 - z)*log(z) + 216*CF*z*log(1 - z)*log(z) - 
       72*TF*z*log(1 - z)*log(z) + 108*CA*z*log(z)*log(1 + z) - 216*CF*z*log(z)*log(1 + z) - 24*TF*pow(Pi,2) + 
       27*CA*z*pow(Pi,2) - 54*CF*z*pow(Pi,2) + 12*TF*z*pow(Pi,2) - 18*CA*z*log(1 + z)*pow(Pi,2) + 36*CF*z*log(1 + z)*pow(Pi,2) - 
       72*TF*pow(z,2) - 378*CA*log(z)*pow(z,2) + 756*CF*log(z)*pow(z,2) - 108*TF*log(z)*pow(z,2) + 
       216*CA*log(z)*log(1 + z)*pow(z,2) - 432*CF*log(z)*log(1 + z)*pow(z,2) + 18*CA*pow(Pi,2)*pow(z,2) - 
       36*CF*pow(Pi,2)*pow(z,2) - 36*Li2(z)*(-3*(CA - 2*CF)*z*log(z)*(1 + pow(z,2)) + 
          (-1 + pow(z,2))*(4*TF - 3*CA*z + 6*CF*z - 2*TF*z + 4*TF*pow(z,2))) + 405*CA*pow(z,3) - 810*CF*pow(z,3) + 
       286*TF*pow(z,3) + 162*CA*Zeta3*pow(z,3) - 324*CF*Zeta3*pow(z,3) - 324*CA*Li3(-z)*pow(z,3) + 648*CF*Li3(-z)*pow(z,3) - 
       216*CA*Li3(z)*pow(z,3) + 432*CF*Li3(z)*pow(z,3) - 216*CA*Li3(1/(1 + z))*pow(z,3) + 432*CF*Li3(1/(1 + z))*pow(z,3) - 
       297*CA*log(z)*pow(z,3) + 594*CF*log(z)*pow(z,3) + 24*TF*log(z)*pow(z,3) + 108*CA*log(1 - z)*log(z)*pow(z,3) - 
       216*CF*log(1 - z)*log(z)*pow(z,3) + 72*TF*log(1 - z)*log(z)*pow(z,3) + 108*CA*log(z)*log(1 + z)*pow(z,3) - 
       216*CF*log(z)*log(1 + z)*pow(z,3) - 9*CA*pow(Pi,2)*pow(z,3) + 18*CF*pow(Pi,2)*pow(z,3) - 12*TF*pow(Pi,2)*pow(z,3) - 
       18*CA*log(1 + z)*pow(Pi,2)*pow(z,3) + 36*CF*log(1 + z)*pow(Pi,2)*pow(z,3) - 272*TF*pow(z,4) + 384*TF*log(z)*pow(z,4) - 
       144*TF*log(1 - z)*log(z)*pow(z,4) + 24*TF*pow(Pi,2)*pow(z,4) + 
       108*(CA - 2*CF)*z*Li2(-z)*(log(z)*(1 + pow(z,2)) + pow(1 + z,2)) - 27*TF*z*pow(log(z),2) - 
       54*CA*z*log(1 + z)*pow(log(z),2) + 108*CF*z*log(1 + z)*pow(log(z),2) - 54*TF*pow(z,2)*pow(log(z),2) - 
       99*TF*pow(z,3)*pow(log(z),2) - 54*CA*log(1 + z)*pow(z,3)*pow(log(z),2) + 108*CF*log(1 + z)*pow(z,3)*pow(log(z),2) - 
       72*TF*pow(z,4)*pow(log(z),2) + 9*CA*z*pow(log(z),3) - 18*CF*z*pow(log(z),3) + 18*TF*z*pow(log(z),3) + 
       36*TF*pow(z,2)*pow(log(z),3) + 9*CA*pow(z,3)*pow(log(z),3) - 18*CF*pow(z,3)*pow(log(z),3) + 
       18*TF*pow(z,3)*pow(log(z),3) + 36*CA*z*pow(log(1 + z),3) - 72*CF*z*pow(log(1 + z),3) + 36*CA*pow(z,3)*pow(log(1 + z),3) - 
       72*CF*pow(z,3)*pow(log(1 + z),3)))/(27.*(1 + z)*pow(z,2));
}

double I2qiqj(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*TF*(2*(-1 + z)*(-172 + 143*z - 136*pow(z,2) + 6*pow(Pi,2)*(2 - z + 2*pow(z,2))) - 
       72*Li2(z)*(-2 + 3*z - 3*pow(z,2) + 2*pow(z,3)) - 
       12*log(z)*(z*(-21 + 30*z - 32*pow(z,2)) + 6*log(1 - z)*(-2 + 3*z - 3*pow(z,2) + 2*pow(z,3))) - 
       9*z*(3 + 3*z + 8*pow(z,2))*pow(log(z),2) + 18*z*(1 + z)*pow(log(z),3)))/(27.*pow(z,2));
}

double I2qg(double z, double x, double f, double f0) { 
	if (z<x) return 0;
	return (f*TF*(688*CA - 1260*CA*z - 702*CF*z - 648*CA*z*Zeta3 + 1728*CF*z*Zeta3 + 288*CA*Li2(z) - 432*CA*z*Li2(z) + 
       216*CA*z*Li3(1 - z) - 216*CF*z*Li3(1 - z) + 648*CA*z*Li3(-z) - 216*CF*z*Li3(z) + 432*CA*z*Li3(1/(1 + z)) + 
       504*CA*z*log(z) + 432*CF*z*log(z) + 216*CF*z*Li2(z)*log(z) + 288*CA*log(1 - z)*log(z) - 432*CA*z*log(1 - z)*log(z) - 
       48*CA*pow(Pi,2) + 72*CA*z*pow(Pi,2) + 36*CA*z*log(1 + z)*pow(Pi,2) + 1548*CA*pow(z,2) + 4050*CF*pow(z,2) - 
       3456*CF*Zeta3*pow(z,2) + 1728*CA*Li2(z)*pow(z,2) - 432*CA*Li3(1 - z)*pow(z,2) + 432*CF*Li3(1 - z)*pow(z,2) + 
       1296*CA*Li3(-z)*pow(z,2) + 1728*CA*Li3(z)*pow(z,2) + 432*CF*Li3(z)*pow(z,2) + 864*CA*Li3(1/(1 + z))*pow(z,2) + 
       324*CA*log(1 - z)*pow(z,2) - 324*CF*log(1 - z)*pow(z,2) - 720*CA*log(z)*pow(z,2) + 810*CF*log(z)*pow(z,2) - 
       864*CA*Li2(z)*log(z)*pow(z,2) - 432*CF*Li2(z)*log(z)*pow(z,2) + 1728*CA*log(1 - z)*log(z)*pow(z,2) - 
       432*CF*log(1 - z)*log(z)*pow(z,2) + 432*CA*log(z)*log(1 + z)*pow(z,2) - 216*CA*pow(Pi,2)*pow(z,2) - 
       108*CF*pow(Pi,2)*pow(z,2) + 72*CA*log(1 + z)*pow(Pi,2)*pow(z,2) - 
       216*(CA - CF)*z*Li2(1 - z)*log(1 - z)*(1 - 2*z + 2*pow(z,2)) - 
       216*CA*z*Li2(-z)*(-2*z*(1 + z) + log(z)*(1 + 2*z + 2*pow(z,2))) - 1192*CA*pow(z,3) - 3888*CF*pow(z,3) - 
       1296*CA*Zeta3*pow(z,3) + 3456*CF*Zeta3*pow(z,3) - 1584*CA*Li2(z)*pow(z,3) + 432*CA*Li3(1 - z)*pow(z,3) - 
       432*CF*Li3(1 - z)*pow(z,3) + 1296*CA*Li3(-z)*pow(z,3) - 432*CF*Li3(z)*pow(z,3) + 864*CA*Li3(1/(1 + z))*pow(z,3) - 
       432*CA*log(1 - z)*pow(z,3) + 432*CF*log(1 - z)*pow(z,3) + 1632*CA*log(z)*pow(z,3) - 432*CF*log(z)*pow(z,3) + 
       432*CF*Li2(z)*log(z)*pow(z,3) - 1584*CA*log(1 - z)*log(z)*pow(z,3) + 432*CF*log(1 - z)*log(z)*pow(z,3) + 
       432*CA*log(z)*log(1 + z)*pow(z,3) + 264*CA*pow(Pi,2)*pow(z,3) + 108*CF*pow(Pi,2)*pow(z,3) + 
       72*CA*log(1 + z)*pow(Pi,2)*pow(z,3) + 108*CF*z*log(z)*pow(log(1 - z),2) - 216*CA*pow(z,2)*pow(log(1 - z),2) + 
       216*CF*pow(z,2)*pow(log(1 - z),2) - 216*CF*log(z)*pow(z,2)*pow(log(1 - z),2) + 216*CA*pow(z,3)*pow(log(1 - z),2) - 
       216*CF*pow(z,3)*pow(log(1 - z),2) + 216*CF*log(z)*pow(z,3)*pow(log(1 - z),2) + 36*CA*z*pow(log(1 - z),3) - 
       36*CF*z*pow(log(1 - z),3) - 72*CA*pow(z,2)*pow(log(1 - z),3) + 72*CF*pow(z,2)*pow(log(1 - z),3) + 
       72*CA*pow(z,3)*pow(log(1 - z),3) - 72*CF*pow(z,3)*pow(log(1 - z),3) - 54*CA*z*pow(log(z),2) + 27*CF*z*pow(log(z),2) + 
       108*CF*z*log(1 - z)*pow(log(z),2) + 108*CA*z*log(1 + z)*pow(log(z),2) + 216*CA*pow(z,2)*pow(log(z),2) + 
       324*CF*pow(z,2)*pow(log(z),2) - 216*CF*log(1 - z)*pow(z,2)*pow(log(z),2) + 216*CA*log(1 + z)*pow(z,2)*pow(log(z),2) - 
       792*CA*pow(z,3)*pow(log(z),2) - 216*CF*pow(z,3)*pow(log(z),2) + 216*CF*log(1 - z)*pow(z,3)*pow(log(z),2) + 
       216*CA*log(1 + z)*pow(z,3)*pow(log(z),2) + 36*CA*z*pow(log(z),3) - 18*CF*z*pow(log(z),3) + 72*CA*pow(z,2)*pow(log(z),3) + 
       36*CF*pow(z,2)*pow(log(z),3) - 72*CF*pow(z,3)*pow(log(z),3) - 72*CA*z*pow(log(1 + z),3) - 
       144*CA*pow(z,2)*pow(log(1 + z),3) - 144*CA*pow(z,3)*pow(log(1 + z),3)))/(54.*pow(z,2));
}

double I2gg(double z, double x, double f, double f0) {
	if (z<x) return (CA*f0*(96*(56*nf*TF + CA*(-202 + 189*Zeta3))*log(1 - x) + 8*nf*TF*(-328 + 252*Zeta3 + 45*pow(Pi,2)) + 
       CA*(9712 - 5544*Zeta3 - 1206*pow(Pi,2) + 45*pow(Pi,4))))/648.;
	return (CA*f0*(96*(56*nf*TF + CA*(-202 + 189*Zeta3))*log(1 - x) + 8*nf*TF*(-328 + 252*Zeta3 + 45*pow(Pi,2)) + 
        CA*(9712 - 5544*Zeta3 - 1206*pow(Pi,2) + 45*pow(Pi,4))))/648. - 
   (-4*CA*f0*(1 + z)*(56*nf*TF + CA*(-202 + 189*Zeta3))*pow(z,2) + 
      f*(484*CA*nf*TF - 72*CF*nf*TF - 440*CA*nf*TF*z + 1728*CF*nf*TF*z + 156*CA*nf*TF*z*log(z) + 648*CF*nf*TF*z*log(z) - 
         3160*pow(CA,2) + 2896*z*pow(CA,2) + 1728*Zeta3*pow(CA,2) - 1080*z*Zeta3*pow(CA,2) - 648*Li3(-z)*pow(CA,2) - 
         648*z*Li3(-z)*pow(CA,2) - 1080*Li3(z)*pow(CA,2) + 216*z*Li3(z)*pow(CA,2) - 432*Li3(1/(1 + z))*pow(CA,2) - 
         432*z*Li3(1/(1 + z))*pow(CA,2) - 2103*z*log(z)*pow(CA,2) - 792*log(1 - z)*log(z)*pow(CA,2) + 
         864*z*log(1 - z)*log(z)*pow(CA,2) + 132*pow(CA,2)*pow(Pi,2) - 144*z*pow(CA,2)*pow(Pi,2) - 
         36*log(1 + z)*pow(CA,2)*pow(Pi,2) - 36*z*log(1 + z)*pow(CA,2)*pow(Pi,2) + 180*CA*nf*TF*pow(z,2) - 
         1656*CF*nf*TF*pow(z,2) - 36*CA*nf*TF*log(1 - z)*pow(z,2) + 120*CA*nf*TF*log(z)*pow(z,2) + 
         648*CF*nf*TF*log(z)*pow(z,2) - 688*pow(CA,2)*pow(z,2) + 1728*Zeta3*pow(CA,2)*pow(z,2) - 
         648*Li3(-z)*pow(CA,2)*pow(z,2) - 1080*Li3(z)*pow(CA,2)*pow(z,2) - 432*Li3(1/(1 + z))*pow(CA,2)*pow(z,2) + 
         18*log(1 - z)*pow(CA,2)*pow(z,2) - 447*log(z)*pow(CA,2)*pow(z,2) - 72*log(1 - z)*log(z)*pow(CA,2)*pow(z,2) + 
         12*pow(CA,2)*pow(Pi,2)*pow(z,2) - 36*log(1 + z)*pow(CA,2)*pow(Pi,2)*pow(z,2) + 108*CA*nf*TF*pow(z,3) - 
         1656*CF*nf*TF*pow(z,3) - 156*CA*nf*TF*log(z)*pow(z,3) - 648*CF*nf*TF*log(z)*pow(z,3) - 364*pow(CA,2)*pow(z,3) + 
         1080*Zeta3*pow(CA,2)*pow(z,3) + 648*Li3(-z)*pow(CA,2)*pow(z,3) - 216*Li3(z)*pow(CA,2)*pow(z,3) + 
         432*Li3(1/(1 + z))*pow(CA,2)*pow(z,3) + 495*log(z)*pow(CA,2)*pow(z,3) - 72*log(1 - z)*log(z)*pow(CA,2)*pow(z,3) + 
         12*pow(CA,2)*pow(Pi,2)*pow(z,3) + 36*log(1 + z)*pow(CA,2)*pow(Pi,2)*pow(z,3) - 440*CA*nf*TF*pow(z,4) + 
         1728*CF*nf*TF*pow(z,4) + 36*CA*nf*TF*log(1 - z)*pow(z,4) - 120*CA*nf*TF*log(z)*pow(z,4) - 
         648*CF*nf*TF*log(z)*pow(z,4) + 3040*pow(CA,2)*pow(z,4) - 1728*Zeta3*pow(CA,2)*pow(z,4) + 
         648*Li3(-z)*pow(CA,2)*pow(z,4) + 1080*Li3(z)*pow(CA,2)*pow(z,4) + 432*Li3(1/(1 + z))*pow(CA,2)*pow(z,4) - 
         18*log(1 - z)*pow(CA,2)*pow(z,4) + 447*log(z)*pow(CA,2)*pow(z,4) + 864*log(1 - z)*log(z)*pow(CA,2)*pow(z,4) - 
         144*pow(CA,2)*pow(Pi,2)*pow(z,4) + 36*log(1 + z)*pow(CA,2)*pow(Pi,2)*pow(z,4) + 556*CA*nf*TF*pow(z,5) - 
         72*CF*nf*TF*pow(z,5) - 3340*pow(CA,2)*pow(z,5) + 1080*Zeta3*pow(CA,2)*pow(z,5) + 648*Li3(-z)*pow(CA,2)*pow(z,5) - 
         216*Li3(z)*pow(CA,2)*pow(z,5) + 432*Li3(1/(1 + z))*pow(CA,2)*pow(z,5) + 1608*log(z)*pow(CA,2)*pow(z,5) - 
         792*log(1 - z)*log(z)*pow(CA,2)*pow(z,5) + 132*pow(CA,2)*pow(Pi,2)*pow(z,5) + 
         36*log(1 + z)*pow(CA,2)*pow(Pi,2)*pow(z,5) + 
         72*Li2(z)*pow(CA,2)*(-(pow(-1 + z,2)*(11 + 10*z + 10*pow(z,2) + 11*pow(z,3))) + 
            3*log(z)*(3 - z + 3*pow(z,2) + pow(z,3) - 3*pow(z,4) + pow(z,5))) + 108*log(z)*pow(CA,2)*pow(log(1 - z),2) - 
         108*z*log(z)*pow(CA,2)*pow(log(1 - z),2) + 108*log(z)*pow(CA,2)*pow(z,2)*pow(log(1 - z),2) + 
         108*log(z)*pow(CA,2)*pow(z,3)*pow(log(1 - z),2) - 108*log(z)*pow(CA,2)*pow(z,4)*pow(log(1 - z),2) + 
         108*log(z)*pow(CA,2)*pow(z,5)*pow(log(1 - z),2) + 36*CA*nf*TF*z*pow(log(z),2) + 162*CF*nf*TF*z*pow(log(z),2) + 
         225*z*pow(CA,2)*pow(log(z),2) + 108*log(1 - z)*pow(CA,2)*pow(log(z),2) - 108*z*log(1 - z)*pow(CA,2)*pow(log(z),2) - 
         108*log(1 + z)*pow(CA,2)*pow(log(z),2) - 108*z*log(1 + z)*pow(CA,2)*pow(log(z),2) + 
         36*CA*nf*TF*pow(z,2)*pow(log(z),2) + 54*CF*nf*TF*pow(z,2)*pow(log(z),2) - 99*pow(CA,2)*pow(z,2)*pow(log(z),2) + 
         108*log(1 - z)*pow(CA,2)*pow(z,2)*pow(log(z),2) - 108*log(1 + z)*pow(CA,2)*pow(z,2)*pow(log(z),2) - 
         36*CA*nf*TF*pow(z,3)*pow(log(z),2) - 162*CF*nf*TF*pow(z,3)*pow(log(z),2) + 171*pow(CA,2)*pow(z,3)*pow(log(z),2) + 
         108*log(1 - z)*pow(CA,2)*pow(z,3)*pow(log(z),2) + 108*log(1 + z)*pow(CA,2)*pow(z,3)*pow(log(z),2) - 
         36*CA*nf*TF*pow(z,4)*pow(log(z),2) - 54*CF*nf*TF*pow(z,4)*pow(log(z),2) + 99*pow(CA,2)*pow(z,4)*pow(log(z),2) - 
         108*log(1 - z)*pow(CA,2)*pow(z,4)*pow(log(z),2) + 108*log(1 + z)*pow(CA,2)*pow(z,4)*pow(log(z),2) - 
         396*pow(CA,2)*pow(z,5)*pow(log(z),2) + 108*log(1 - z)*pow(CA,2)*pow(z,5)*pow(log(z),2) + 
         108*log(1 + z)*pow(CA,2)*pow(z,5)*pow(log(z),2) + 36*CF*nf*TF*z*pow(log(z),3) - 36*z*pow(CA,2)*pow(log(z),3) + 
         36*CF*nf*TF*pow(z,2)*pow(log(z),3) - 72*pow(CA,2)*pow(z,2)*pow(log(z),3) - 36*CF*nf*TF*pow(z,3)*pow(log(z),3) + 
         36*pow(CA,2)*pow(z,3)*pow(log(z),3) - 36*CF*nf*TF*pow(z,4)*pow(log(z),3) + 72*pow(CA,2)*pow(z,4)*pow(log(z),3) - 
         36*pow(CA,2)*pow(z,5)*pow(log(z),3) + 72*pow(CA,2)*pow(log(1 + z),3) + 72*z*pow(CA,2)*pow(log(1 + z),3) + 
         72*pow(CA,2)*pow(z,2)*pow(log(1 + z),3) - 72*pow(CA,2)*pow(z,3)*pow(log(1 + z),3) - 
         72*pow(CA,2)*pow(z,4)*pow(log(1 + z),3) - 72*pow(CA,2)*pow(z,5)*pow(log(1 + z),3) - 
         216*(-1 + z)*Li2(-z)*log(z)*pow(CA,2)*pow(1 + z + pow(z,2),2)))/(27.*(-1 + pow(z,2))*pow(z,2));
}

double I2gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*(-6320*CA + 896*nf*TF + 6328*CA*z + 540*CF*z - 896*nf*TF*z + 3456*CA*Zeta3 - 2160*CA*z*Zeta3 - 1296*CA*Li3(-z) - 
       1296*CA*z*Li3(-z) - 2160*CA*Li3(z) + 432*CA*z*Li3(z) - 864*CA*Li3(1/(1 + z)) - 864*CA*z*Li3(1/(1 + z)) - 
       1824*CA*log(1 - z) + 1728*CF*log(1 - z) + 480*nf*TF*log(1 - z) + 1824*CA*z*log(1 - z) - 1728*CF*z*log(1 - z) - 
       480*nf*TF*z*log(1 - z) - 2988*CA*z*log(z) - 810*CF*z*log(z) - 1584*CA*log(1 - z)*log(z) + 1728*CA*z*log(1 - z)*log(z) + 
       264*CA*pow(Pi,2) - 288*CA*z*pow(Pi,2) - 72*CA*log(1 + z)*pow(Pi,2) - 72*CA*z*log(1 + z)*pow(Pi,2) - 2144*CA*pow(z,2) - 
       54*CF*pow(z,2) + 208*nf*TF*pow(z,2) + 1728*CA*Zeta3*pow(z,2) - 648*CA*Li3(-z)*pow(z,2) - 1080*CA*Li3(z)*pow(z,2) - 
       432*CA*Li3(1/(1 + z))*pow(z,2) - 516*CA*log(1 - z)*pow(z,2) + 540*CF*log(1 - z)*pow(z,2) + 96*nf*TF*log(1 - z)*pow(z,2) + 
       72*CA*log(z)*pow(z,2) + 270*CF*log(z)*pow(z,2) - 648*CA*log(1 - z)*log(z)*pow(z,2) + 216*CA*log(z)*log(1 + z)*pow(z,2) + 
       54*CA*pow(Pi,2)*pow(z,2) - 36*CA*log(1 + z)*pow(Pi,2)*pow(z,2) + 
       216*CA*Li2(-z)*(pow(z,2) + log(z)*(2 + 2*z + pow(z,2))) + 1216*CA*pow(z,3) - 1056*CA*log(z)*pow(z,3) + 
       288*CA*log(1 - z)*log(z)*pow(z,3) - 48*CA*pow(Pi,2)*pow(z,3) + 
       72*CA*Li2(z)*(-22 + 24*z - 6*pow(z,2) + 3*log(z)*(6 - 2*z + 3*pow(z,2)) + 4*pow(z,3)) - 396*CA*pow(log(1 - z),2) + 
       324*CF*pow(log(1 - z),2) + 144*nf*TF*pow(log(1 - z),2) + 396*CA*z*pow(log(1 - z),2) - 324*CF*z*pow(log(1 - z),2) - 
       144*nf*TF*z*pow(log(1 - z),2) + 216*CA*log(z)*pow(log(1 - z),2) - 216*CA*z*log(z)*pow(log(1 - z),2) - 
       90*CA*pow(z,2)*pow(log(1 - z),2) + 54*CF*pow(z,2)*pow(log(1 - z),2) + 72*nf*TF*pow(z,2)*pow(log(1 - z),2) + 
       108*CA*log(z)*pow(z,2)*pow(log(1 - z),2) - 72*CA*pow(log(1 - z),3) + 72*CF*pow(log(1 - z),3) + 
       72*CA*z*pow(log(1 - z),3) - 72*CF*z*pow(log(1 - z),3) - 36*CA*pow(z,2)*pow(log(1 - z),3) + 
       36*CF*pow(z,2)*pow(log(1 - z),3) + 648*CA*z*pow(log(z),2) - 108*CF*z*pow(log(z),2) + 216*CA*log(1 - z)*pow(log(z),2) - 
       216*CA*z*log(1 - z)*pow(log(z),2) - 216*CA*log(1 + z)*pow(log(z),2) - 216*CA*z*log(1 + z)*pow(log(z),2) + 
       162*CA*pow(z,2)*pow(log(z),2) - 81*CF*pow(z,2)*pow(log(z),2) + 108*CA*log(1 - z)*pow(z,2)*pow(log(z),2) - 
       108*CA*log(1 + z)*pow(z,2)*pow(log(z),2) + 144*CA*pow(z,3)*pow(log(z),2) - 72*CA*z*pow(log(z),3) + 
       36*CF*z*pow(log(z),3) - 36*CA*pow(z,2)*pow(log(z),3) - 18*CF*pow(z,2)*pow(log(z),3) + 144*CA*pow(log(1 + z),3) + 
       144*CA*z*pow(log(1 + z),3) + 72*CA*pow(z,2)*pow(log(1 + z),3)))/(54.*pow(z,2));
}

//-----------------------------------------------

double P0qiqiP0qiqi(double z, double x, double f, double f0) {
	if (z<x) return f0*pow(CF,2)*(2.25 + 6*log(1 - x) - (2*pow(Pi,2))/3. + 4*pow(log(1 - x),2));
	return -((pow(CF,2)*(-2*f0*z*(3 + 4*log(1 - z)) + f*(1 + 4*z + pow(z,2) + 4*log(1 - z)*(1 + pow(z,2)) - log(z)*(1 + 3*pow(z,2)))))/
      ((-1 + z)*z)) + f0*pow(CF,2)*(2.25 + 6*log(1 - x) - (2*pow(Pi,2))/3. + 4*pow(log(1 - x),2));
}

double P0qgP0gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*TF*(4 + 3*z + 6*z*(1 + z)*log(z) - 3*pow(z,2) - 4*pow(z,3)))/(3.*pow(z,2));
}

double P0qgP0gg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (f*TF*(8*CA + 3*beta0*z + 6*CA*z + 12*CA*z*(1 + 4*z)*log(z) - 6*beta0*pow(z,2) + 48*CA*pow(z,2) + 
       12*CA*z*log(1 - z)*(1 - 2*z + 2*pow(z,2)) + 6*beta0*pow(z,3) - 62*CA*pow(z,3)))/(6.*pow(z,2)); 
}

double P0qiqiP0qg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*TF*(-1 + 4*z + log(z)*(-2 + 4*z - 8*pow(z,2)) + log(1 - z)*(4 - 8*z + 8*pow(z,2))))/(2.*z);
}

double P0ggP0gg(double z, double x, double f, double f0) {
	if (z<x) return f0*(2*beta0*CA*log(1 - x) + pow(beta0,2)/4. - (2*pow(CA,2)*(pow(Pi,2) - 6*pow(log(1 - x),2)))/3.);
	return f0*(2*beta0*CA*log(1 - x) + pow(beta0,2)/4. - (2*pow(CA,2)*(pow(Pi,2) - 6*pow(log(1 - x),2)))/3.) - 
   (2*CA*(-3*f0*(beta0 + 4*CA*log(1 - z))*pow(z,2) + 
        f*(3*beta0 - 22*CA - 6*beta0*z + 40*CA*z + 9*beta0*pow(z,2) - 36*CA*pow(z,2) - 6*beta0*pow(z,3) + 40*CA*pow(z,3) + 
           3*beta0*pow(z,4) - 22*CA*pow(z,4) - 6*CA*log(z)*(1 + 3*pow(z,2) - 4*pow(z,3) + pow(z,4)) + 
           12*CA*log(1 - z)*pow(1 - z + pow(z,2),2))))/(3.*(-1 + z)*pow(z,2));
}

double P0gqP0qg(double z, double x, double f, double f0) {
	return 2*nf*P0qgP0gq(z,x,f,f0);
}

double P0gqP0qiqi(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (f*pow(CF,2)*(-(z*(-4 + z + 2*(-2 + z)*log(z))) + 4*log(1 - z)*(2 - 2*z + pow(z,2))))/(2.*pow(z,2));
}

double P0ggP0gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*(6*beta0 - 62*CA - 6*beta0*z + 48*CA*z + 3*beta0*pow(z,2) + 6*CA*pow(z,2) + 12*CA*log(1 - z)*(2 - 2*z + pow(z,2)) - 
       24*CA*log(z)*(1 + z + pow(z,2)) + 8*CA*pow(z,3)))/(6.*pow(z,2));
}

double I1qiqiP0qiqi(double z, double x, double f, double f0) {
	if (z<x) return -(f0*(3 + 4*log(1 - x))*pow(CF,2)*pow(Pi,2))/12.;
	return -(f0*(3 + 4*log(1 - x))*pow(CF,2)*pow(Pi,2))/12. + (pow(CF,2)*
      (-2*f0*z*pow(Pi,2) + f*(6 - 12*z + pow(Pi,2) - 24*log(1 - z)*pow(-1 + z,2) + 12*log(z)*pow(-1 + z,2) + 6*pow(z,2) + 
           pow(Pi,2)*pow(z,2))))/(6.*(-1 + z)*z);
}

double I1qgP0gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (4*CF*f*TF*(-3/z + 2*z - 3*log(z) + 1/pow(z,2)))/3.;
}

double I1qgP0gg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (-2*f*TF*(12*CA*(-1 + z)*log(1 - z)*pow(z,2) + 24*CA*log(z)*pow(z,2) - 
       (-1 + z)*(-3*beta0*pow(z,2) + CA*(-2 + 4*z + 34*pow(z,2)))))/(3.*pow(z,2));
}

double I1qiqiP0qg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return (CF*f*TF*(-12*(1 + 2*z)*log(z) + pow(Pi,2)*(-1 + 2*z - 2*pow(z,2)) + 12*(-2 + z + pow(z,2))))/(6.*z);
}

double I1ggP0gg(double z, double x, double f, double f0) {
	if (z<x) return -(CA*f0*(beta0 + 4*CA*log(1 - x))*pow(Pi,2))/12.;
	return -(CA*f0*(beta0 + 4*CA*log(1 - x))*pow(Pi,2))/12. + (pow(CA,2)*pow(Pi,2)*(-(f0*pow(z,2)) + f*pow(1 - z + pow(z,2),2)))/
    (3.*(-1 + z)*pow(z,2));
}

double I1gqP0qg(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return 2*nf*CF*f*TF*(2 + 2/z - 4*z + 4*log(z));
}

double I1gqP0qiqi(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return f*(1 + 2/z + 4*log(1 - z) - 2*log(z))*pow(CF,2);
}

double I1ggP0gq(double z, double x, double f, double f0) {
	if (z<x) return 0;
	return -(CA*CF*f*pow(Pi,2)*(2 - 2*z + pow(z,2)))/(6.*pow(z,2));
}

//=================================================================================================

double Cqq0qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	return PDF(i1,x1,mu)*PDF(i2,x2,mu);
}

double Cqq1qg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fq0=PDF(i1,x1,mu);
	double fg=PDF(0,x2/z2,mu);
	double fg0=PDF(0,x2,mu);
	return fq0*(I1qg(z2,x2,fg,fg0) + 2*(L + Lmu)*P0qg(z2,x2,fg,fg0));
}

double Cqq1qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fqa=PDF(i1,x1/z1,mu);
	double fqa0=PDF(i1,x1,mu);
	double fqb=PDF(i2,x2/z2,mu);
	double fqb0=PDF(i2,x2,mu);
	return (fqa0*fqb0*(4*cH1qq + L*(2*CF*d1 + 2*gH0qq - CF*G0*L)))/2. + fqb0*I1qiqi(z1,x1,fqa,fqa0) + fqa0*I1qiqi(z2,x2,fqb,fqb0) + 
               2*(L + Lmu)*(fqb0*P0qiqi(z1,x1,fqa,fqa0) + fqa0*P0qiqi(z2,x2,fqb,fqb0));
}

double Cqq2qiqj(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fq0=PDF(i1,x1,mu);
	double fqj=0;
	double fqj0=0;
	for (int id=1; id<=Nf; id++) { if (id!=i2&&id!=-i2) { fqj += PDF(id,x2/z2,mu)+PDF(-id,x2/z2,mu); fqj0 += PDF(id,x2,mu)+PDF(-id,x2,mu); } }
	return fq0*(I2qiqj(z2,x2,fqj,fqj0) + 2*(L + Lmu)*(I1qgP0gq(z2,x2,fqj,fqj0) + 2*P1qiqj(z2,x2,fqj,fqj0)) + 
     2*P0qgP0gq(z2,x2,fqj,fqj0)*pow(L + Lmu,2));
}

double Cqq2qiqbi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fq0=PDF(i1,x1,mu);
	double fqb=PDF(-i2,x2/z2,mu);
	double fqb0=PDF(-i2,x2,mu);
	return fq0*(I2qiqbi(z2,x2,fqb,fqb0) + 2*(L + Lmu)*(I1qgP0gq(z2,x2,fqb,fqb0) + 2*P1qiqbi(z2,x2,fqb,fqb0)) + 
     2*P0qgP0gq(z2,x2,fqb,fqb0)*pow(L + Lmu,2));
}

double Cqq2qg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fq=PDF(i1,x1/z1,mu);
	double fq0=PDF(i1,x1,mu);
	double fg=PDF(0,x2/z2,mu);
	double fg0=PDF(0,x2,mu);
	return (4*I1qg(z2,x2,fg,fg0)*I1qiqi(z1,x1,fq,fq0) + 2*fq0*(4*cH1qq - Lmu*(2*gH0qq + CF*G0*Lmu))*
      (I1qg(z2,x2,fg,fg0) + 2*(L + Lmu)*P0qg(z2,x2,fg,fg0)) + 
     2*(L + Lmu)*(4*I1qiqi(z1,x1,fq,fq0)*P0qg(z2,x2,fg,fg0) + I1qg(z2,x2,fg,fg0)*(fq0*gH0qq + 4*P0qiqi(z1,x1,fq,fq0))) + 
     4*CF*fq0*(L*(d1 - G0*(L + Lmu))*I1qg(z2,x2,fg,fg0) + 
        2*P0qg(z2,x2,fg,fg0)*(d1*L*(L + Lmu) - G0*(4*Zeta3 + 2*Lmu*pow(L,2) + pow(L,3) + L*pow(Lmu,2)))) + 
     (CF*fq0*G0*I1qg(z2,x2,fg,fg0) + 4*P0qg(z2,x2,fg,fg0)*(fq0*gH0qq + 4*P0qiqi(z1,x1,fq,fq0)))*pow(L + Lmu,2) + 
     fq0*(4*I2qg(z2,x2,fg,fg0) - 2*(L + Lmu)*((2*beta0 - gH0qq)*I1qg(z2,x2,fg,fg0) - 
           4*(I1qgP0gg(z2,x2,fg,fg0) + I1qiqiP0qg(z2,x2,fg,fg0) + 2*P1qg(z2,x2,fg,fg0))) + 
        (CF*G0*I1qg(z2,x2,fg,fg0) - 4*(beta0 - gH0qq)*P0qg(z2,x2,fg,fg0) + 
           8*(P0qgP0gg(z2,x2,fg,fg0) + P0qiqiP0qg(z2,x2,fg,fg0)))*pow(L + Lmu,2) - 
        2*CF*G0*P0qg(z2,x2,fg,fg0)*(-4*Zeta3 - pow(L + Lmu,3))) - 2*CF*fq0*G0*P0qg(z2,x2,fg,fg0)*(-4*Zeta3 - pow(L + Lmu,3)))/4.;
}

double Cqq2gg(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fga=PDF(0,x1/z1,mu);
	double fga0=PDF(0,x1,mu);
	double fgb=PDF(0,x2/z2,mu);
	double fgb0=PDF(0,x2,mu);
	return (I1qg(z1,x1,fga,fga0) + 2*(L + Lmu)*P0qg(z1,x1,fga,fga0))*(I1qg(z2,x2,fgb,fgb0) + 2*(L + Lmu)*P0qg(z2,x2,fgb,fgb0));
}

double Cqq2qiqi(int i1, int i2, double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fqa=PDF(i1,x1/z1,mu);
	double fqa0=PDF(i1,x1,mu);
	double fqb=PDF(i2,x2/z2,mu);
	double fqb0=PDF(i2,x2,mu);
	return 2*cH2qq*fqa0*fqb0 + 2*CF*cH1qq*d1*fqa0*fqb0*L + CF*d2*fqa0*fqb0*L + 2*cH1qq*fqa0*fqb0*gH0qq*L + fqa0*fqb0*gH1qq*L - 
   2*beta0*cH1qq*fqa0*fqb0*Lmu - beta0*CF*d1*fqa0*fqb0*L*Lmu - beta0*fqa0*fqb0*gH0qq*L*Lmu + 
   (4*beta0*CF*fqa0*fqb0*G0*Zeta3)/3. - 2*CF*fqa0*fqb0*G0*gH0qq*Zeta3 + 2*fqb0*(L + Lmu)*I1qgP0gq(z1,x1,fqa,fqa0) + 
   2*fqa0*(L + Lmu)*I1qgP0gq(z2,x2,fqb,fqb0) + 2*cH1qq*fqb0*I1qiqi(z1,x1,fqa,fqa0) - beta0*fqb0*L*I1qiqi(z1,x1,fqa,fqa0) + 
   CF*d1*fqb0*L*I1qiqi(z1,x1,fqa,fqa0) + fqb0*gH0qq*L*I1qiqi(z1,x1,fqa,fqa0) - beta0*fqb0*Lmu*I1qiqi(z1,x1,fqa,fqa0) + 
   2*cH1qq*fqa0*I1qiqi(z2,x2,fqb,fqb0) - beta0*fqa0*L*I1qiqi(z2,x2,fqb,fqb0) + CF*d1*fqa0*L*I1qiqi(z2,x2,fqb,fqb0) + 
   fqa0*gH0qq*L*I1qiqi(z2,x2,fqb,fqb0) - beta0*fqa0*Lmu*I1qiqi(z2,x2,fqb,fqb0) + I1qiqi(z1,x1,fqa,fqa0)*I1qiqi(z2,x2,fqb,fqb0) + 
   2*fqb0*L*I1qiqiP0qiqi(z1,x1,fqa,fqa0) + 2*fqb0*Lmu*I1qiqiP0qiqi(z1,x1,fqa,fqa0) + 2*fqa0*L*I1qiqiP0qiqi(z2,x2,fqb,fqb0) + 
   2*fqa0*Lmu*I1qiqiP0qiqi(z2,x2,fqb,fqb0) + fqb0*I2qiqi(z1,x1,fqa,fqa0) + fqa0*I2qiqi(z2,x2,fqb,fqb0) + 
   4*fqb0*L*Lmu*P0qgP0gq(z1,x1,fqa,fqa0) + 4*fqa0*L*Lmu*P0qgP0gq(z2,x2,fqb,fqb0) + 4*cH1qq*fqb0*L*P0qiqi(z1,x1,fqa,fqa0) + 
   4*cH1qq*fqb0*Lmu*P0qiqi(z1,x1,fqa,fqa0) - 2*beta0*fqb0*L*Lmu*P0qiqi(z1,x1,fqa,fqa0) + 
   2*CF*d1*fqb0*L*Lmu*P0qiqi(z1,x1,fqa,fqa0) + 2*fqb0*gH0qq*L*Lmu*P0qiqi(z1,x1,fqa,fqa0) - 
   4*CF*fqb0*G0*Zeta3*P0qiqi(z1,x1,fqa,fqa0) + 2*L*I1qiqi(z2,x2,fqb,fqb0)*P0qiqi(z1,x1,fqa,fqa0) + 
   2*Lmu*I1qiqi(z2,x2,fqb,fqb0)*P0qiqi(z1,x1,fqa,fqa0) + 4*cH1qq*fqa0*L*P0qiqi(z2,x2,fqb,fqb0) + 
   4*cH1qq*fqa0*Lmu*P0qiqi(z2,x2,fqb,fqb0) - 2*beta0*fqa0*L*Lmu*P0qiqi(z2,x2,fqb,fqb0) + 
   2*CF*d1*fqa0*L*Lmu*P0qiqi(z2,x2,fqb,fqb0) + 2*fqa0*gH0qq*L*Lmu*P0qiqi(z2,x2,fqb,fqb0) - 
   4*CF*fqa0*G0*Zeta3*P0qiqi(z2,x2,fqb,fqb0) + 2*L*I1qiqi(z1,x1,fqa,fqa0)*P0qiqi(z2,x2,fqb,fqb0) + 
   2*Lmu*I1qiqi(z1,x1,fqa,fqa0)*P0qiqi(z2,x2,fqb,fqb0) + 8*L*Lmu*P0qiqi(z1,x1,fqa,fqa0)*P0qiqi(z2,x2,fqb,fqb0) + 
   4*fqb0*L*Lmu*P0qiqiP0qiqi(z1,x1,fqa,fqa0) + 4*fqa0*L*Lmu*P0qiqiP0qiqi(z2,x2,fqb,fqb0) + 4*fqb0*L*P1qiqi(z1,x1,fqa,fqa0) + 
   4*fqb0*Lmu*P1qiqi(z1,x1,fqa,fqa0) + 4*fqa0*L*P1qiqi(z2,x2,fqb,fqb0) + 4*fqa0*Lmu*P1qiqi(z2,x2,fqb,fqb0) - 
   2*d1*fqa0*fqb0*G0*Zeta3*pow(CF,2) + 2*fqa0*fqb0*L*Zeta3*pow(CF,2)*pow(G0,2) - beta0*CF*d1*fqa0*fqb0*pow(L,2) - 
   CF*cH1qq*fqa0*fqb0*G0*pow(L,2) - (CF*fqa0*fqb0*G1*pow(L,2))/2. - (beta0*fqa0*fqb0*gH0qq*pow(L,2))/2. + 
   CF*d1*fqa0*fqb0*gH0qq*pow(L,2) + (beta0*CF*fqa0*fqb0*G0*Lmu*pow(L,2))/2. - (CF*fqb0*G0*I1qiqi(z1,x1,fqa,fqa0)*pow(L,2))/2. - 
   (CF*fqa0*G0*I1qiqi(z2,x2,fqb,fqb0)*pow(L,2))/2. + 2*fqb0*P0qgP0gq(z1,x1,fqa,fqa0)*pow(L,2) + 
   2*fqa0*P0qgP0gq(z2,x2,fqb,fqb0)*pow(L,2) - beta0*fqb0*P0qiqi(z1,x1,fqa,fqa0)*pow(L,2) + 
   2*CF*d1*fqb0*P0qiqi(z1,x1,fqa,fqa0)*pow(L,2) + 2*fqb0*gH0qq*P0qiqi(z1,x1,fqa,fqa0)*pow(L,2) - 
   CF*fqb0*G0*Lmu*P0qiqi(z1,x1,fqa,fqa0)*pow(L,2) - beta0*fqa0*P0qiqi(z2,x2,fqb,fqb0)*pow(L,2) + 
   2*CF*d1*fqa0*P0qiqi(z2,x2,fqb,fqb0)*pow(L,2) + 2*fqa0*gH0qq*P0qiqi(z2,x2,fqb,fqb0)*pow(L,2) - 
   CF*fqa0*G0*Lmu*P0qiqi(z2,x2,fqb,fqb0)*pow(L,2) + 4*P0qiqi(z1,x1,fqa,fqa0)*P0qiqi(z2,x2,fqb,fqb0)*pow(L,2) + 
   2*fqb0*P0qiqiP0qiqi(z1,x1,fqa,fqa0)*pow(L,2) + 2*fqa0*P0qiqiP0qiqi(z2,x2,fqb,fqb0)*pow(L,2) + 
   (fqa0*fqb0*pow(CF,2)*pow(d1,2)*pow(L,2))/2. + (fqa0*fqb0*pow(gH0qq,2)*pow(L,2))/2. + (beta0*CF*fqa0*fqb0*G0*pow(L,3))/3. - 
   (CF*fqa0*fqb0*G0*gH0qq*pow(L,3))/2. - CF*fqb0*G0*P0qiqi(z1,x1,fqa,fqa0)*pow(L,3) - 
   CF*fqa0*G0*P0qiqi(z2,x2,fqb,fqb0)*pow(L,3) - (d1*fqa0*fqb0*G0*pow(CF,2)*pow(L,3))/2. + 
   (fqa0*fqb0*pow(CF,2)*pow(G0,2)*pow(L,4))/8. + 2*fqb0*P0qgP0gq(z1,x1,fqa,fqa0)*pow(Lmu,2) + 
   2*fqa0*P0qgP0gq(z2,x2,fqb,fqb0)*pow(Lmu,2) - beta0*fqb0*P0qiqi(z1,x1,fqa,fqa0)*pow(Lmu,2) - 
   beta0*fqa0*P0qiqi(z2,x2,fqb,fqb0)*pow(Lmu,2) + 4*P0qiqi(z1,x1,fqa,fqa0)*P0qiqi(z2,x2,fqb,fqb0)*pow(Lmu,2) + 
   2*fqb0*P0qiqiP0qiqi(z1,x1,fqa,fqa0)*pow(Lmu,2) + 2*fqa0*P0qiqiP0qiqi(z2,x2,fqb,fqb0)*pow(Lmu,2);
}

//=================================================================================================

double Cgg0gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	return PDF(0,x1,mu)*PDF(0,x2,mu);
}

double Cgg1gq(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fg0=PDF(0,x1,mu);
	double fq=0;
	double fq0=0;
	for (int id=1; id<=Nf; id++) { fq += PDF(id,x2/z2,mu)+PDF(-id,x2/z2,mu); fq0 += PDF(id,x2,mu)+PDF(-id,x2,mu); }
	return fg0*(I1gq(z2,x2,fq,fq0) + 2*(L + Lmu)*P0gq(z2,x2,fq,fq0));
}

double Cgg1gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fga=PDF(0,x1/z1,mu);
	double fga0=PDF(0,x1,mu);
	double fgb=PDF(0,x2/z2,mu);
	double fgb0=PDF(0,x2,mu);
	return (fga0*fgb0*(4*cH1gg + 2*gH0gg*L + CA*L*(2*d1 - G0*L) - 4*beta0*Lmu))/2. + fgb0*I1gg(z1,x1,fga,fga0) + fga0*I1gg(z2,x2,fgb,fgb0) + 
               2*(L + Lmu)*(fgb0*P0gg(z1,x1,fga,fga0) + fga0*P0gg(z2,x2,fgb,fgb0));
}

double Cgg2gq(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fg=PDF(0,x1/z1,mu);
	double fg0=PDF(0,x1,mu);
	double fq=0;
	double fq0=0;
	for (int id=1; id<=Nf; id++) { fq += PDF(id,x2/z2,mu)+PDF(-id,x2/z2,mu); fq0 += PDF(id,x2,mu)+PDF(-id,x2,mu); }
	return (4*CA*fg0*L*(d1 - G0*(L + Lmu))*I1gq(z2,x2,fq,fq0) + 4*I1gg(z1,x1,fg,fg0)*I1gq(z2,x2,fq,fq0) + 4*fg0*I2gq(z2,x2,fq,fq0) + 
     4*Ii1gg(z1,x1,fg,fg0)*Ii1gq(z2,x2,fq,fq0) + 8*CA*d1*fg0*L*(L + Lmu)*P0gq(z2,x2,fq,fq0) - 
     16*CA*fg0*G0*Zeta3*P0gq(z2,x2,fq,fq0) + 2*fg0*(4*cH1gg - Lmu*(4*beta0 + 2*gH0gg + CA*G0*Lmu))*
      (I1gq(z2,x2,fq,fq0) + 2*(L + Lmu)*P0gq(z2,x2,fq,fq0)) + 
     2*(L + Lmu)*(I1gq(z2,x2,fq,fq0)*(fg0*gH0gg + 4*P0gg(z1,x1,fg,fg0)) + 4*I1gg(z1,x1,fg,fg0)*P0gq(z2,x2,fq,fq0)) - 
     2*fg0*(L + Lmu)*((2*beta0 - gH0gg)*I1gq(z2,x2,fq,fq0) - 4*(I1ggP0gq(z2,x2,fq,fq0) + I1gqP0qiqi(z2,x2,fq,fq0)) - 
        8*P1gq(z2,x2,fq,fq0)) - 16*CA*fg0*G0*Lmu*P0gq(z2,x2,fq,fq0)*pow(L,2) - 8*CA*fg0*G0*P0gq(z2,x2,fq,fq0)*pow(L,3) - 
     8*CA*fg0*G0*L*P0gq(z2,x2,fq,fq0)*pow(Lmu,2) + 
     (CA*fg0*G0*I1gq(z2,x2,fq,fq0) + 4*(fg0*gH0gg + 4*P0gg(z1,x1,fg,fg0))*P0gq(z2,x2,fq,fq0))*pow(L + Lmu,2) + 
     fg0*(CA*G0*I1gq(z2,x2,fq,fq0) + 4*(-beta0 + gH0gg)*P0gq(z2,x2,fq,fq0) + 
        8*(P0ggP0gq(z2,x2,fq,fq0) + P0gqP0qiqi(z2,x2,fq,fq0)))*pow(L + Lmu,2) + 4*CA*fg0*G0*P0gq(z2,x2,fq,fq0)*pow(L + Lmu,3))/4.;
}

double Cgg2qq(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fqa=0;
	double fqa0=0;
	double fqb=0;
	double fqb0=0;
	for (int id=1; id<=Nf; id++) { fqa += PDF(id,x1/z1,mu)+PDF(-id,x1/z1,mu); fqa0 += PDF(id,x1,mu)+PDF(-id,x1,mu); }
	for (int id=1; id<=Nf; id++) { fqb += PDF(id,x2/z2,mu)+PDF(-id,x2/z2,mu); fqb0 += PDF(id,x2,mu)+PDF(-id,x2,mu); }
	return Ii1gq(z1,x1,fqa,fqa0)*Ii1gq(z2,x2,fqb,fqb0) +
               (I1gq(z1,x1,fqa,fqa0) + 2*(L + Lmu)*P0gq(z1,x1,fqa,fqa0))*(I1gq(z2,x2,fqb,fqb0) + 2*(L + Lmu)*P0gq(z2,x2,fqb,fqb0));
}

double Cgg2gg(double x1, double x2, double z1, double z2, double qT, double mu, double Q) {
	double L=2.*log(qT/Q);
	double Lmu=2.*log(Q/mu);
	double fga=PDF(0,x1/z1,mu);
	double fga0=PDF(0,x1,mu);
	double fgb=PDF(0,x2/z2,mu);
	double fgb0=PDF(0,x2,mu);
	return 2*fgb0*L*I1ggP0gg(z1,x1,fga,fga0) + 2*fgb0*Lmu*I1ggP0gg(z1,x1,fga,fga0) + 2*fga0*L*I1ggP0gg(z2,x2,fgb,fgb0) + 
   2*fga0*Lmu*I1ggP0gg(z2,x2,fgb,fgb0) + 2*fgb0*L*I1gqP0qg(z1,x1,fga,fga0) + 2*fgb0*Lmu*I1gqP0qg(z1,x1,fga,fga0) + 
   2*fga0*L*I1gqP0qg(z2,x2,fgb,fgb0) + 2*fga0*Lmu*I1gqP0qg(z2,x2,fgb,fgb0) + fgb0*I2gg(z1,x1,fga,fga0) + 
   fga0*I2gg(z2,x2,fgb,fgb0) + Ii1gg(z1,x1,fga,fga0)*Ii1gg(z2,x2,fgb,fgb0) + 4*cH1gg*fgb0*L*P0gg(z1,x1,fga,fga0) + 
   4*cH1gg*fgb0*Lmu*P0gg(z1,x1,fga,fga0) - 6*beta0*fgb0*L*Lmu*P0gg(z1,x1,fga,fga0) + 2*CA*d1*fgb0*L*Lmu*P0gg(z1,x1,fga,fga0) + 
   2*fgb0*gH0gg*L*Lmu*P0gg(z1,x1,fga,fga0) + (I1gg(z2,x2,fgb,fgb0)*
      (fga0*(4*cH1gg + 2*gH0gg*L + CA*L*(2*d1 - G0*L) - 2*beta0*(L + 3*Lmu)) + 4*(L + Lmu)*P0gg(z1,x1,fga,fga0)))/2. + 
   4*cH1gg*fga0*L*P0gg(z2,x2,fgb,fgb0) + 4*cH1gg*fga0*Lmu*P0gg(z2,x2,fgb,fgb0) - 6*beta0*fga0*L*Lmu*P0gg(z2,x2,fgb,fgb0) + 
   2*CA*d1*fga0*L*Lmu*P0gg(z2,x2,fgb,fgb0) + 2*fga0*gH0gg*L*Lmu*P0gg(z2,x2,fgb,fgb0) + 
   8*L*Lmu*P0gg(z1,x1,fga,fga0)*P0gg(z2,x2,fgb,fgb0) - 
   (2*CA*G0*Zeta3*(fga0*fgb0*(-2*beta0 + 3*CA*d1 + 3*gH0gg - 3*CA*G0*L) + 6*fgb0*P0gg(z1,x1,fga,fga0) + 
        6*fga0*P0gg(z2,x2,fgb,fgb0)))/3. + (I1gg(z1,x1,fga,fga0)*
      (fgb0*(4*cH1gg + 2*gH0gg*L + CA*L*(2*d1 - G0*L) - 2*beta0*(L + 3*Lmu)) + 2*I1gg(z2,x2,fgb,fgb0) + 
        4*(L + Lmu)*P0gg(z2,x2,fgb,fgb0)))/2. + 4*fgb0*L*Lmu*P0ggP0gg(z1,x1,fga,fga0) + 
   4*fga0*L*Lmu*P0ggP0gg(z2,x2,fgb,fgb0) + 4*fgb0*L*Lmu*P0gqP0qg(z1,x1,fga,fga0) + 4*fga0*L*Lmu*P0gqP0qg(z2,x2,fgb,fgb0) + 
   4*fgb0*L*P1gg(z1,x1,fga,fga0) + 4*fgb0*Lmu*P1gg(z1,x1,fga,fga0) + 4*fga0*L*P1gg(z2,x2,fgb,fgb0) + 
   4*fga0*Lmu*P1gg(z2,x2,fgb,fgb0) - beta0*fgb0*P0gg(z1,x1,fga,fga0)*pow(L,2) + 2*CA*d1*fgb0*P0gg(z1,x1,fga,fga0)*pow(L,2) + 
   2*fgb0*gH0gg*P0gg(z1,x1,fga,fga0)*pow(L,2) - CA*fgb0*G0*Lmu*P0gg(z1,x1,fga,fga0)*pow(L,2) - 
   beta0*fga0*P0gg(z2,x2,fgb,fgb0)*pow(L,2) + 2*CA*d1*fga0*P0gg(z2,x2,fgb,fgb0)*pow(L,2) + 
   2*fga0*gH0gg*P0gg(z2,x2,fgb,fgb0)*pow(L,2) - CA*fga0*G0*Lmu*P0gg(z2,x2,fgb,fgb0)*pow(L,2) + 
   4*P0gg(z1,x1,fga,fga0)*P0gg(z2,x2,fgb,fgb0)*pow(L,2) + 2*fgb0*P0ggP0gg(z1,x1,fga,fga0)*pow(L,2) + 
   2*fga0*P0ggP0gg(z2,x2,fgb,fgb0)*pow(L,2) + 2*fgb0*P0gqP0qg(z1,x1,fga,fga0)*pow(L,2) + 
   2*fga0*P0gqP0qg(z2,x2,fgb,fgb0)*pow(L,2) - CA*fgb0*G0*P0gg(z1,x1,fga,fga0)*pow(L,3) - 
   CA*fga0*G0*P0gg(z2,x2,fgb,fgb0)*pow(L,3) - 5*beta0*fgb0*P0gg(z1,x1,fga,fga0)*pow(Lmu,2) - 
   5*beta0*fga0*P0gg(z2,x2,fgb,fgb0)*pow(Lmu,2) + 4*P0gg(z1,x1,fga,fga0)*P0gg(z2,x2,fgb,fgb0)*pow(Lmu,2) + 
   2*fgb0*P0ggP0gg(z1,x1,fga,fga0)*pow(Lmu,2) + 2*fga0*P0ggP0gg(z2,x2,fgb,fgb0)*pow(Lmu,2) + 
   2*fgb0*P0gqP0qg(z1,x1,fga,fga0)*pow(Lmu,2) + 2*fga0*P0gqP0qg(z2,x2,fgb,fgb0)*pow(Lmu,2) + 
   (fga0*fgb0*(48*cH2gg + 4*CA*L*(6*d2 + 6*d1*gH0gg*L + 6*cH1gg*(2*d1 - G0*L) - 6*beta0*d1*(L + 3*Lmu) + 
           L*(-3*G1 + G0*(2*beta0 - 3*gH0gg)*L + 9*beta0*G0*Lmu)) + 3*pow(CA,2)*pow(L,2)*pow(-2*d1 + G0*L,2) + 
        12*(L*(2*gH1gg + gH0gg*(-beta0 + gH0gg)*L) - 2*(4*beta1 + 3*beta0*gH0gg*L)*Lmu + 4*cH1gg*(gH0gg*L - 3*beta0*Lmu) + 
           6*pow(beta0,2)*pow(Lmu,2))))/24.;
}

}
