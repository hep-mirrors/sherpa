#ifndef NNLO__param_H
#define NNLO__param_H

#include "ATOOLS/Math/MyComplex.H"
#include <math.h>
#include <iomanip>
#include <stdlib.h> // optional
#include "ATOOLS/Org/Message.H"

namespace SHNNLO {

static double Pi = 3.14159265358979323846;
static double Zeta2 = 1.6449340668482264;
static double Zeta3 = 1.2020569031595942;
static double Zeta4 = 1.0823232337111382;
static double Zeta5 = 1.0369277551433699;
static double CA = 3.;
static double CF = 4./3.;
static double TF = 1./2.;

void updateparam();

extern double G_F,Mh,Mt,Mb,Mc,Mw,Mz,Gw,Gz,sw2;
extern double e2,Qu,Qd,eZ,guL,gdL,guR,gdR;

extern double nf,sumQsq;
extern int Nf;

extern double beta0,beta1,beta2;
extern double G0,G1,G2;
extern double d1,d2,d3;
extern double gS0,gS1,gS2;

extern double gB0q,gB1q,gB2q,gB0g,gB1g,gB2g;
extern double gH0qq,gH1qq,gH2qq,cH1qq,cH2qq,cH3qq;
extern double gH0gg,gH1gg,gH2gg,cH1gg,cH2gg,cH3gg;

extern double cS1,cS2;
extern double deltaqT1,I1hat,I2hat;

}

#endif
