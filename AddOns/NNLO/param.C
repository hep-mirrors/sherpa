#include "param.H"

namespace SHNNLO {

double G_F,Mh,Mt,Mb,Mc,Mw,Mz,Gw,Gz,sw2;
double e2,Qu,Qd,eZ,guL,gdL,guR,gdR;

double nf = 5.;
double sumQsq = (nf>=6.?3.:2.)*4./9.+(nf>=5.?3.:2)*1./9.; // sum of active quark charge squared in the loop
int Nf = 5; // number of active flavors for PDF

double beta0 = 11./3.*CA-4./3.*TF*nf;
double beta1 = 34./3.*pow(CA,2)-(20./3.*CA+4.*CF)*TF*nf;
double beta2 = 2857./54.*pow(CA,3)
      +(2.*pow(CF,2)-205./9.*CF*CA-1415./27.*pow(CA,2))*TF*nf
      +(44./9.*CF+158./27.*CA)*pow(TF*nf,2);

double G0 = 4;
double G1 = 4*(CA*(67./9. - pow(Pi,2)/3.) - (20.*nf*TF)/9.);
double G2 = 4*((-16*pow(nf,2)*pow(TF,2))/27. + CA*nf*TF*(-418/27. + (40*pow(Pi,2))/27. - (56*Zeta3)/3.) + 
     pow(CA,2)*(245/6. - (134*pow(Pi,2))/27. + (11*pow(Pi,4))/45. + (22*Zeta3)/3.) + CF*nf*TF*(-55/3. + 16*Zeta3));

double gS0 = 0;
double gS1 = nf*(448./27. - (8.*pow(Pi,2))/9.)*TF +
    CA*(-1616./27. + (22.*pow(Pi,2))/9. + 56.*Zeta3);
double gS2 = pow(CA,2)*(-273562/729. + (2632*Zeta3)/3. - 384*Zeta5 + (12650*pow(Pi,2))/243. - (176*Zeta3*pow(Pi,2))/9. - (176*pow(Pi,4))/45.) + 
   CF*nf*TF*(6844/27. - (1216*Zeta3)/9. - (8*pow(Pi,2))/3. - (32*pow(Pi,4))/45.) + 
   CA*nf*TF*(47368/729. - (2912*Zeta3)/27. - (5656*pow(Pi,2))/243. + (32*pow(Pi,4))/15.) + 
   pow(nf,2)*(16640/729. - (896*Zeta3)/27. + (160*pow(Pi,2))/81.)*pow(TF,2);

double d1 = 0;
double d2 = (-224*nf*TF)/27. + CA*(808/27. - 28*Zeta3);
double d3 = pow(nf,2)*pow(TF,2)*(7424/729. + (128*Zeta3)/9.) + CA*nf*TF*(-125252/729. + (824*pow(Pi,2))/243. - (4*pow(Pi,4))/27. + (1808*Zeta3)/27.) + 
   CF*nf*TF*(-3422/27. + (16*pow(Pi,4))/45. + (608*Zeta3)/9.) + 
   pow(CA,2)*(297029/729. - (3196*pow(Pi,2))/243. - (77*pow(Pi,4))/135. - (12328*Zeta3)/27. + (88*pow(Pi,2)*Zeta3)/9. + 192*Zeta5);

double gH0qq = -6*CF;
double gH1qq = CF*(CA*(-961/27. + 52*Zeta3 - (11*pow(Pi,2))/3.) + nf*TF*(260/27. + (4*pow(Pi,2))/3.) + CF*(-3 - 48*Zeta3 + 4*pow(Pi,2)));
double gH2qq = pow(CF,3)*(-29 - 136*Zeta3 + 480*Zeta5 - 6*pow(Pi,2) + (32*Zeta3*pow(Pi,2))/3. - (16*pow(Pi,4))/5.) + 
   nf*TF*pow(CF,2)*(5906/27. + (1024*Zeta3)/9. - (52*pow(Pi,2))/9. - (56*pow(Pi,4))/27.) + 
   CF*pow(CA,2)*(-139345/1458. + (7052*Zeta3)/9. - 272*Zeta5 - (7163*pow(Pi,2))/243. - (88*Zeta3*pow(Pi,2))/9. - (83*pow(Pi,4))/45.) + 
   CA*CF*nf*TF*(-34636/729. - (3856*Zeta3)/27. + (5188*pow(Pi,2))/243. + (44*pow(Pi,4))/45.) + 
   CA*pow(CF,2)*(-151/2. - (1688*Zeta3)/3. - 240*Zeta5 + (410*pow(Pi,2))/9. - (16*Zeta3*pow(Pi,2))/3. + (494*pow(Pi,4))/135.) + 
   CF*pow(nf,2)*(19336/729. - (64*Zeta3)/27. - (80*pow(Pi,2))/27.)*pow(TF,2);
double gH0gg = -2*beta0;
double gH1gg = -4*beta1 - 8*CF*nf*TF + CA*nf*TF*(-208/27. - (4*pow(Pi,2))/9.) + pow(CA,2)*(160/27. + 4*Zeta3 + (11*pow(Pi,2))/9.);
double gH2gg = -6*beta2 + 8*nf*TF*pow(CF,2) + pow(CA,3)*(37045/729. - 32*Zeta5 + Zeta3*(244/3. - (40*pow(Pi,2))/9.) + (6109*pow(Pi,2))/243. - (319*pow(Pi,4))/135.) + 
   CA*CF*nf*TF*(1178/27. - (608*Zeta3)/9. - (4*pow(Pi,2))/3. - (16*pow(Pi,4))/45.) + 
   nf*TF*pow(CA,2)*(-167800/729. + (1424*Zeta3)/27. - (2396*pow(Pi,2))/243. + (164*pow(Pi,4))/135.) + (176*CF*pow(nf,2)*pow(TF,2))/9. + 
   CA*pow(nf,2)*(24520/729. - (448*Zeta3)/27. + (80*pow(Pi,2))/81.)*pow(TF,2) ;

double gB0q = -gH0qq - (CF*gS0)/2.;
double gB1q = -gH1qq - (CF*gS1)/2.;
double gB2q = -gH2qq - (CF*gS2)/2.;
double gB0g = -gH0gg - (CA*gS0)/2.;
double gB1g = -gH1gg - (CA*gS1)/2.;
double gB2g = -gH2gg - (CA*gS2)/2.;

double cH1qq = CF*(-8 + (7*pow(Pi,2))/6.);
double cH2qq = CF*(CF*(511/8. - (83*pow(Pi,2))/6. + (67*pow(Pi,4))/60. - 30*Zeta3) +
                      nf*TF*(4085/162. - (91*pow(Pi,2))/27. + (4*Zeta3)/9.) +
                      CA*(-51157/648. + (1061*pow(Pi,2))/108. - (4*pow(Pi,4))/45. + (313*Zeta3)/9.));
double cH3qq = 0.;
double cH1gg = (7*CA*pow(Pi,2))/6.;
double cH2gg = CF*nf*TF*(-67./3. + 16*Zeta3) + CA*nf*TF*(-1832./81. - (92*Zeta3)/9. - (25*pow(Pi,2))/9.) +
   pow(CA,2)*(5105./162. - (143*Zeta3)/9. + (335*pow(Pi,2))/36. + (37*pow(Pi,4))/36.);
double cH3gg = 0.;

double cS1 = pow(Pi,2)/3.;
double cS2 = nf*TF*(80./81. + (74.*pow(Pi,2))/27. - (232.*Zeta3)/9.) +
    CA*(-2140./81. - (335.*pow(Pi,2))/54. + (22.*pow(Pi,4))/45. + (638.*Zeta3)/9.);

//IR subtraction constants in Eq. (6.11) of 1503.04812
double deltaqT1 = 10./3.*Zeta3*beta0/2. + (-1214./81.+67./18.*Zeta2)*CA + (164./81.-5./9.*Zeta2)*nf;
double I1hat = Pi*Pi*CF/6.;
double I2hat = 1./72.*CF*(pow(Pi,4)*CF + 36*deltaqT1 + 12*Pi*Pi*G1/8. + 48*beta0/2.*Zeta3);

void updateparam() {
	// couplings
	e2 = 8.*G_F/sqrt(2.)*Mw*Mw*sw2;
	Qu = 2./3.;
	Qd = -1./3.;
	eZ = (1.-sw2)*sqrt(e2/sw2/(1.-sw2));
	guL = eZ/2./(1-sw2)*(1.0/2-2.0/3*sw2);
	gdL = eZ/2./(1-sw2)*(-1.0/2+1.0/3*sw2);
	guR = -eZ/2./(1-sw2)*2.0/3*sw2;
	gdR = eZ/2./(1-sw2)*1.0/3*sw2;
	// other constants
	sumQsq = (nf>=6.?3.:2.)*4./9.+(nf>=5.?3.:2)*1./9.; // sum of active quark charge squared in the loop
	Nf = nf; // number of active flavors for PDF
	beta0 = 11./3.*CA-4./3.*TF*nf;
	beta1 = 34./3.*pow(CA,2)-(20./3.*CA+4.*CF)*TF*nf;
	beta2 = 2857./54.*pow(CA,3)
      +(2.*pow(CF,2)-205./9.*CF*CA-1415./27.*pow(CA,2))*TF*nf
      +(44./9.*CF+158./27.*CA)*pow(TF*nf,2);
	G0 = 4;
	G1 = 4*(CA*(67./9. - pow(Pi,2)/3.) - (20.*nf*TF)/9.);
	G2 = 4*((-16*pow(nf,2)*pow(TF,2))/27. + CA*nf*TF*(-418/27. + (40*pow(Pi,2))/27. - (56*Zeta3)/3.) +
     pow(CA,2)*(245/6. - (134*pow(Pi,2))/27. + (11*pow(Pi,4))/45. + (22*Zeta3)/3.) + CF*nf*TF*(-55/3. + 16*Zeta3));
	gS0 = 0;
	gS1 = nf*(448./27. - (8.*pow(Pi,2))/9.)*TF +
    CA*(-1616./27. + (22.*pow(Pi,2))/9. + 56.*Zeta3);
	gS2 = pow(CA,2)*(-273562/729. + (2632*Zeta3)/3. - 384*Zeta5 + (12650*pow(Pi,2))/243. - (176*Zeta3*pow(Pi,2))/9. - (176*pow(Pi,4))/45.) +
   CF*nf*TF*(6844/27. - (1216*Zeta3)/9. - (8*pow(Pi,2))/3. - (32*pow(Pi,4))/45.) +
   CA*nf*TF*(47368/729. - (2912*Zeta3)/27. - (5656*pow(Pi,2))/243. + (32*pow(Pi,4))/15.) +
   pow(nf,2)*(16640/729. - (896*Zeta3)/27. + (160*pow(Pi,2))/81.)*pow(TF,2);
	d1 = 0;
	d2 = (-224*nf*TF)/27. + CA*(808/27. - 28*Zeta3);
	d3 = pow(nf,2)*pow(TF,2)*(7424/729. + (128*Zeta3)/9.) + CA*nf*TF*(-125252/729. + (824*pow(Pi,2))/243. - (4*pow(Pi,4))/27. + (1808*Zeta3)/27.) +
   CF*nf*TF*(-3422/27. + (16*pow(Pi,4))/45. + (608*Zeta3)/9.) +
   pow(CA,2)*(297029/729. - (3196*pow(Pi,2))/243. - (77*pow(Pi,4))/135. - (12328*Zeta3)/27. + (88*pow(Pi,2)*Zeta3)/9. + 192*Zeta5);
	gH0qq = -6*CF;
	gH1qq = CF*(CA*(-961/27. + 52*Zeta3 - (11*pow(Pi,2))/3.) + nf*TF*(260/27. + (4*pow(Pi,2))/3.) + CF*(-3 - 48*Zeta3 + 4*pow(Pi,2)));
	gH2qq = pow(CF,3)*(-29 - 136*Zeta3 + 480*Zeta5 - 6*pow(Pi,2) + (32*Zeta3*pow(Pi,2))/3. - (16*pow(Pi,4))/5.) +
   nf*TF*pow(CF,2)*(5906/27. + (1024*Zeta3)/9. - (52*pow(Pi,2))/9. - (56*pow(Pi,4))/27.) +
   CF*pow(CA,2)*(-139345/1458. + (7052*Zeta3)/9. - 272*Zeta5 - (7163*pow(Pi,2))/243. - (88*Zeta3*pow(Pi,2))/9. - (83*pow(Pi,4))/45.) +
   CA*CF*nf*TF*(-34636/729. - (3856*Zeta3)/27. + (5188*pow(Pi,2))/243. + (44*pow(Pi,4))/45.) +
   CA*pow(CF,2)*(-151/2. - (1688*Zeta3)/3. - 240*Zeta5 + (410*pow(Pi,2))/9. - (16*Zeta3*pow(Pi,2))/3. + (494*pow(Pi,4))/135.) +
   CF*pow(nf,2)*(19336/729. - (64*Zeta3)/27. - (80*pow(Pi,2))/27.)*pow(TF,2);
	gH0gg = -2*beta0;
	gH1gg = -4*beta1 - 8*CF*nf*TF + CA*nf*TF*(-208/27. - (4*pow(Pi,2))/9.) + pow(CA,2)*(-160/27. + 4*Zeta3 + (11*pow(Pi,2))/9.);
	gH2gg = -6*beta2 + 8*nf*TF*pow(CF,2) + pow(CA,3)*(37045/729. - 32*Zeta5 + Zeta3*(244/3. - (40*pow(Pi,2))/9.) + (6109*pow(Pi,2))/243. - (319*pow(Pi,4))/135.) +
   CA*CF*nf*TF*(1178/27. - (608*Zeta3)/9. - (4*pow(Pi,2))/3. - (16*pow(Pi,4))/45.) +
   nf*TF*pow(CA,2)*(-167800/729. + (1424*Zeta3)/27. - (2396*pow(Pi,2))/243. + (164*pow(Pi,4))/135.) + (176*CF*pow(nf,2)*pow(TF,2))/9. +
   CA*pow(nf,2)*(24520/729. - (448*Zeta3)/27. + (80*pow(Pi,2))/81.)*pow(TF,2) ;
	gB0q = -gH0qq - (CF*gS0)/2.;
	gB1q = -gH1qq - (CF*gS1)/2.;
	gB2q = -gH2qq - (CF*gS2)/2.;
	gB0g = -gH0gg - (CA*gS0)/2.;
	gB1g = -gH1gg - (CA*gS1)/2.;
	gB2g = -gH2gg - (CA*gS2)/2.;
	cH1qq = CF*(-8 + (7*pow(Pi,2))/6.);
	cH2qq = CF*(CF*(511/8. - (83*pow(Pi,2))/6. + (67*pow(Pi,4))/60. - 30*Zeta3) +
                      nf*TF*(4085/162. - (91*pow(Pi,2))/27. + (4*Zeta3)/9.) +
                      CA*(-51157/648. + (1061*pow(Pi,2))/108. - (4*pow(Pi,4))/45. + (313*Zeta3)/9.));
	cH3qq = 0.;
	cH1gg = (7*CA*pow(Pi,2))/6.;
	cH2gg = CF*nf*TF*(-67./3. + 16*Zeta3) + CA*nf*TF*(-1832./81. - (92*Zeta3)/9. - (25*pow(Pi,2))/9.) +
   pow(CA,2)*(5105./162. - (143*Zeta3)/9. + (335*pow(Pi,2))/36. + (37*pow(Pi,4))/36.);
	cH3gg = 0.;
	cS1 = pow(Pi,2)/3.;
	cS2 = nf*TF*(80./81. + (74.*pow(Pi,2))/27. - (232.*Zeta3)/9.) +
    CA*(-2140./81. - (335.*pow(Pi,2))/54. + (22.*pow(Pi,4))/45. + (638.*Zeta3)/9.);
	deltaqT1 = 10./3.*Zeta3*beta0/2. + (-1214./81.+67./18.*Zeta2)*CA + (164./81.-5./9.*Zeta2)*nf;
	I1hat = Pi*Pi*CF/6.;
	I2hat = 1./72.*CF*(pow(Pi,4)*CF + 36*deltaqT1 + 12*Pi*Pi*G1/8. + 48*beta0/2.*Zeta3);
};

}
