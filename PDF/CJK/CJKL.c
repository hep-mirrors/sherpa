/* CJKL.f -- translated by f2c (version 20200916).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

Extern struct {
    integer flav;
} flav_;

#define flav_1 flav_

Extern struct {
    doublereal mc, mb;
} mass_;

#define mass_1 mass_

/* Table of constant values */

static integer c__1 = 1;
static integer c__2 = 2;
static integer c__3 = 3;
static doublereal c_b12 = 1.;

/* ******************************************************************** */
/* ***   LO parametrization of the parton densities in the real    **** */
/* ***             photon for the CJKL fit based on                **** */
/* ***                                                             **** */
/* ***        F.Cornet, P.Jankowski, M.Krawczyk and A.Lorca        **** */
/* ***      "A New 5 Flavour LO Analysis and Parametrization       **** */
/* ***         of Parton Distributions in the Real Photon"         **** */
/* ***               Phys. Rev. D68: 014010, 2003                  **** */
/* ***                      hep-ph/0212160                         **** */
/* ******************************************************************** */
/* ***                                                             **** */
/* ***   valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2      **** */
/* ***                                                             **** */
/* ***       x      - Bjorken x variable                           **** */
/* ***       Q2     - square of momentum scale (in GeV**2)         **** */
/* ***   XPDF(-5:5) - matrix containing x*f(x,Q2)/alfa             **** */
/* ***                                                             **** */
/* ***   PDF =   -5 ,  -4 ,  -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5          **** */
/* ***         b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b          **** */
/* ***                                                             **** */
/* ***   heavy quarks masses:                                      **** */
/* ***             M_charm=1.3, M_beauty=4.3 GeV                   **** */
/* ***                                                             **** */
/* ***   Lambda_QCD values for active N_f falvours:                **** */
/* ***   N_f=   3      4      5                                    **** */
/* ***        0.314  0.280  0.221 GeV                              **** */
/* ***                                                             **** */
/* ***  LIPARTS: only light partons (gl,d,u,s)                     **** */
/* ***  PARTONS: (gl,d,u,s,c,b) parton densities                   **** */
/* ***  PARTF2:  (gl,d,u,s,c,b) parton densities                   **** */
/* ***           & structure function F2 according to Eq. (22)     **** */
/* ***                                                             **** */
/* ******************************************************************** */
/* ***  Evolution, parametrization, checks performed and programs  **** */
/* ***  written by                                                 **** */
/* ***     Pawel Jankowski, (pjank@fuw.edu.pl)                     **** */
/* ***                      Institute of Theoretical Physics,      **** */
/* ***                      Warsaw University                      **** */
/* ***                      ul. Hoza 69                            **** */
/* ***                      00-681 Warsaw                          **** */
/* ***                      Poland                                 **** */
/* ***                                                             **** */
/* ***  Last changes - 09 DEC 2002                                 **** */
/* ******************************************************************** */
/*<       SUBROUTINE LIPARTS(x,Q2,XPDF) >*/
/* Subroutine */ int liparts_(doublereal *x, doublereal *q2, doublereal *xpdf)
{
    static integer i__;
    static doublereal f2, part[11];
    extern /* Subroutine */ int param_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension XPDF(-3:3),PART(-5:5) >*/
/*<       integer I >*/
/* ***  x * par / alfa  **** */
/*<       call PARAM(1,x,Q2,PART,F2) >*/
    /* Parameter adjustments */
    xpdf -= -3;

    /* Function Body */
    param_(&c__1, x, q2, part, &f2);
/*<       do 1 I=-3,3 >*/
    for (i__ = -3; i__ <= 3; ++i__) {
/*<  1       XPDF(I)=PART(I) >*/
/* L1: */
	xpdf[i__] = part[i__ + 5];
    }
/*<       END >*/
    return 0;
} /* liparts_ */

/*<       SUBROUTINE PARTONS(x,Q2,XPDF) >*/
/* Subroutine */ int partons_(doublereal *x, doublereal *q2, doublereal *xpdf)
{
    static doublereal f2;
    extern /* Subroutine */ int param_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension XPDF(-5:5) >*/
/* ***  x * par / alfa  **** */
/*<       call PARAM(2,x,Q2,XPDF,F2) >*/
    /* Parameter adjustments */
    xpdf -= -5;

    /* Function Body */
    param_(&c__2, x, q2, &xpdf[-5], &f2);
/*<       END >*/
    return 0;
} /* partons_ */

/*<       SUBROUTINE PARTF2(x,Q2,XPDF,F2) >*/
/* Subroutine */ int partf2_(doublereal *x, doublereal *q2, doublereal *xpdf, 
	doublereal *f2)
{
    extern /* Subroutine */ int param_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension XPDF(-5:5) >*/
/* ***  x * par / alfa  **** */
/* ***  F2 / alfa       **** */
/*<       call PARAM(3,x,Q2,XPDF,F2) >*/
    /* Parameter adjustments */
    xpdf -= -5;

    /* Function Body */
    param_(&c__3, x, q2, &xpdf[-5], f2);
/*<       END >*/
    return 0;
} /* partf2_ */

/* ********************************* */
/*<       SUBROUTINE PARAM(OPT,x,Q2,XPDF,F2) >*/
/* Subroutine */ int param_(integer *opt, doublereal *x, doublereal *q2, 
	doublereal *xpdf, doublereal *f2)
{
    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static doublereal f2beauty;
    static integer i__;
    static logical t;
    extern /* Subroutine */ int pointlike_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal dn, f2bresover, f2cresover, up, eb2, ec2, eb4, ec4, f2c,
	     f2b;
    extern /* Subroutine */ int hadronlike_(integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static doublereal mb2, mc2, chm, mbq, mcq, als, bot, glu, eps, res[2], 
	    str, lam3, chhd, beta, glhd, bthd, chpl, dnpl, glpl, vlhd, btpl, 
	    sthd;
    static integer step;
    static doublereal uppl, beta2;
    extern doublereal alfas_(doublereal *, integer *, doublereal *);
    static doublereal zetab, zetac, f2alfa;
    extern /* Subroutine */ int intxr_(doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, logical *);
    static doublereal f2bres, f2cres, f2charm, f2bover, f2cover;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<        >*/
/*<       dimension XPDF(-5:5),res(2) >*/
/*<       logical t >*/
/*<       integer OPT,I,step,flav >*/
/*<       common /flav/ flav >*/
/*<       common /mass/ mc,mb >*/
/*<       mc = 1.3d0 >*/
    /* Parameter adjustments */
    xpdf -= -5;

    /* Function Body */
    mass_1.mc = 1.3;
/*<       mb = 4.3d0 >*/
    mass_1.mb = 4.3;
/*<       Lam3 = 0.314d0 >*/
    lam3 = .314;
/*<       call POINTLIKE(OPT,x,Q2,glpl,uppl,dnpl,chpl,btpl) >*/
    pointlike_(opt, x, q2, &glpl, &uppl, &dnpl, &chpl, &btpl);
/*<       call HADRONLIKE(OPT,x,Q2,glhd,vlhd,sthd,chhd,bthd) >*/
    hadronlike_(opt, x, q2, &glhd, &vlhd, &sthd, &chhd, &bthd);
/* ***  x * par / alfa  **** */
/*<       XPDF(0) = glhd + glpl >*/
    xpdf[0] = glhd + glpl;
/*<       XPDF(2) = 0.5d0*vlhd + sthd + uppl >*/
    xpdf[2] = vlhd * .5 + sthd + uppl;
/*<       XPDF(1) = 0.5d0*vlhd + sthd + dnpl >*/
    xpdf[1] = vlhd * .5 + sthd + dnpl;
/*<       XPDF(3) = sthd + dnpl >*/
    xpdf[3] = sthd + dnpl;
/*<       XPDF(4) = chhd + chpl >*/
    xpdf[4] = chhd + chpl;
/*<       XPDF(5) = bthd + btpl >*/
    xpdf[5] = bthd + btpl;
/*<       do 1 I=1,5 >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<  1       XPDF(-I) = XPDF(I) >*/
/* L1: */
	xpdf[-i__] = xpdf[i__];
    }
/*<       do 3 I=-5,5 >*/
    for (i__ = -5; i__ <= 5; ++i__) {
/*<  3       if (XPDF(I).LT.0.d0) XPDF(I)=0.d0 >*/
/* L3: */
	if (xpdf[i__] < 0.) {
	    xpdf[i__] = 0.;
	}
    }
/*<       F2alfa = 0.d0 >*/
    f2alfa = 0.;
/*<       if (OPT.LE.2) goto 2 >*/
    if (*opt <= 2) {
	goto L2;
    }
/*<       glu = alfa*XPDF(0)/x >*/
    glu = xpdf[0] * .00729735308 / *x;
/*<       up  = alfa*XPDF(2)/x >*/
    up = xpdf[2] * .00729735308 / *x;
/*<       dn  = alfa*XPDF(1)/x >*/
    dn = xpdf[1] * .00729735308 / *x;
/*<       str = alfa*XPDF(3)/x >*/
    str = xpdf[3] * .00729735308 / *x;
/*<       chm = alfa*XPDF(4)/x >*/
    chm = xpdf[4] * .00729735308 / *x;
/*<       bot = alfa*XPDF(5)/x >*/
    bot = xpdf[5] * .00729735308 / *x;
/*<       mc2 = mc*mc >*/
    mc2 = mass_1.mc * mass_1.mc;
/*<       ec2 = e2*e2 >*/
    ec2 = .44444444444444442;
/*<       ec4 = ec2*ec2 >*/
    ec4 = ec2 * ec2;
/*<       mb2 = mb*mb >*/
    mb2 = mass_1.mb * mass_1.mb;
/*<       eb2 = e1*e1 >*/
    eb2 = .1111111111111111;
/*<       eb4 = eb2*eb2 >*/
    eb4 = eb2 * eb2;
/*<       zetac = x*(1.d0+4.d0*mc2/Q2) >*/
    zetac = *x * (mc2 * 4. / *q2 + 1.);
/*<       zetab = x*(1.d0+4.d0*mb2/Q2) >*/
    zetab = *x * (mb2 * 4. / *q2 + 1.);
/*<       if (zetac.GE.1.d0) chm = 0.d0 >*/
    if (zetac >= 1.) {
	chm = 0.;
    }
/*<       if (zetab.GE.1.d0) bot = 0.d0 >*/
    if (zetab >= 1.) {
	bot = 0.;
    }
/* ********************************************************************* */
/* ***   The Bethe-Heitler cross-section for gamma*gamma -> ccbar   **** */
/* ********************************************************************* */
/*<       beta2 = 1.d0-4.d0*mc2*x/((1.d0-x)*Q2) >*/
    beta2 = 1. - mc2 * 4. * *x / ((1. - *x) * *q2);
/*<       if ((x.LT.1.d0).AND.(beta2.GT.0.d0)) then >*/
    if (*x < 1. && beta2 > 0.) {
/*<          beta = DSQRT(beta2) >*/
	beta = sqrt(beta2);
/*<          mcQ = 4.d0*mc2/Q2 >*/
	mcq = mc2 * 4. / *q2;
/*<        >*/
	f2c = *x * 3. * ec4 * .00729735308 / 3.1415926535897932 * (beta * (*x 
		* 8. * (1. - *x) - 1. - *x * (1. - *x) * mcq) + (*x * *x + ((
		float)1. - *x) * ((float)1. - *x) + *x * ((float)1. - *x * (
		float)3.) * mcq - *x * *x * mcq * mcq / (float)2.) * log((
		beta + (float)1.) / ((float)1. - beta)));
/*<       else >*/
    } else {
/*<          F2c = 0.d0 >*/
	f2c = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/* ***      Substraction of the overlaping Bethe-Heitler term       **** */
/* ********************************************************************* */
/*<       if (zetac.LT.1.d0) then >*/
    if (zetac < 1.) {
/*<        >*/
	f2cover = zetac * 3. * ec4 * .00729735308 / 3.1415926535897932 * (
		zetac * zetac + (1. - zetac) * (1. - zetac)) * log(*q2 / mc2);
/*<       else >*/
    } else {
/*<          F2cover = 0.d0 >*/
	f2cover = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/* ***   The Bethe-Heitler cross-section for gamma*gamma -> bbbar   **** */
/* ********************************************************************* */
/*<       beta2 = 1.d0-4.d0*mb2*x/((1.d0-x)*Q2) >*/
    beta2 = 1. - mb2 * 4. * *x / ((1. - *x) * *q2);
/*<       if ((x.LT.1.d0).AND.(beta2.GT.0.d0)) then >*/
    if (*x < 1. && beta2 > 0.) {
/*<          beta = DSQRT(beta2) >*/
	beta = sqrt(beta2);
/*<          mbQ = 4.d0*mb2/Q2 >*/
	mbq = mb2 * 4. / *q2;
/*<        >*/
	f2b = *x * 3. * eb4 * .00729735308 / 3.1415926535897932 * (beta * (*x 
		* 8. * (1. - *x) - 1. - *x * (1. - *x) * mbq) + (*x * *x + ((
		float)1. - *x) * ((float)1. - *x) + *x * ((float)1. - *x * (
		float)3.) * mbq - *x * *x * mbq * mbq / (float)2.) * log((
		beta + (float)1.) / ((float)1. - beta)));
/*<       else >*/
    } else {
/*<          F2b = 0.d0 >*/
	f2b = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/* ***      Substraction of the overlaping Bethe-Heitler term       **** */
/* ********************************************************************* */
/*<       if (zetab.LT.1.d0) then >*/
    if (zetab < 1.) {
/*<        >*/
	f2bover = zetab * 3. * eb4 * .00729735308 / 3.1415926535897932 * (
		zetab * zetab + (1. - zetab) * (1. - zetab)) * log(*q2 / mb2);
/*<       else >*/
    } else {
/*<          F2bover = 0.d0 >*/
	f2bover = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/* ***       1:   CHARM - 'resolved' = gamma*G -> ccbar             **** */
/* ***       2:   Substraction of the overlapping term              **** */
/* ********************************************************************* */
/*<       if (zetac.LT.1.d0) then >*/
    if (zetac < 1.) {
/*<         step = 5 >*/
	step = 5;
/*<         eps = 1.d-6 >*/
	eps = 1e-6;
/*<         flav = 1 >*/
	flav_1.flav = 1;
/*<         call intxr(x,step,zetac,Q2,zetac,1.d0,step,eps,res,t) >*/
	intxr_(x, &step, &zetac, q2, &zetac, &c_b12, &step, &eps, res, &t);
/*<         ALS = alfas(Q2,3,Lam3) >*/
	als = alfas_(q2, &c__3, &lam3);
/*<         F2cres = ec2 * ALS/(2.d0*Pi) * res(1) >*/
	f2cres = ec2 * als / 6.2831853071795862 * res[0];
/*<         F2cresover = zetac*ec2 * ALS/Pi * res(2) * DLOG(Q2/mc2) >*/
	f2cresover = zetac * ec2 * als / 3.1415926535897932 * res[1] * log(*
		q2 / mc2);
/*<       else >*/
    } else {
/*<         F2cres = 0.d0 >*/
	f2cres = 0.;
/*<         F2cresover = 0.d0 >*/
	f2cresover = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/* ***       1:   BEAUTY - 'resolved' = gamma*G -> bbbar            **** */
/* ***       2:   Substraction of the overlapping term              **** */
/* ********************************************************************* */
/*<       if (zetab.LT.1.d0) then >*/
    if (zetab < 1.) {
/*<         step = 5 >*/
	step = 5;
/*<         eps = 1.d-6 >*/
	eps = 1e-6;
/*<         flav = 2 >*/
	flav_1.flav = 2;
/*<         call intxr(x,step,zetab,Q2,zetab,1.d0,step,eps,res,t) >*/
	intxr_(x, &step, &zetab, q2, &zetab, &c_b12, &step, &eps, res, &t);
/*<         ALS = alfas(Q2,3,Lam3) >*/
	als = alfas_(q2, &c__3, &lam3);
/*<         F2bres = eb2 * ALS/(2.d0*Pi) * res(1) >*/
	f2bres = eb2 * als / 6.2831853071795862 * res[0];
/*<         F2bresover = zetab*eb2 * ALS/Pi * res(2) * DLOG(Q2/mb2) >*/
	f2bresover = zetab * eb2 * als / 3.1415926535897932 * res[1] * log(*
		q2 / mb2);
/*<       else >*/
    } else {
/*<         F2bres = 0.d0 >*/
	f2bres = 0.;
/*<         F2bresover = 0.d0 >*/
	f2bresover = 0.;
/*<       endif >*/
    }
/* ********************************************************************* */
/*<       F2charm = 2.d0*x*ec2*chm + F2c - F2cover + F2cres - F2cresover >*/
    f2charm = *x * 2. * ec2 * chm + f2c - f2cover + f2cres - f2cresover;
/*<       if (F2charm.LT.0.d0) F2charm = 0.d0 >*/
    if (f2charm < 0.) {
	f2charm = 0.;
    }
/*<       F2beauty = 2.d0*x*eb2*bot + F2b - F2bover + F2bres - F2bresover >*/
    f2beauty = *x * 2. * eb2 * bot + f2b - f2bover + f2bres - f2bresover;
/*<       if (F2beauty.LT.0.d0) F2beauty = 0.d0 >*/
    if (f2beauty < 0.) {
	f2beauty = 0.;
    }
/*<       F2 = 2.d0*x*(ec2*up + eb2*dn + eb2*str) + F2charm + F2beauty >*/
    *f2 = *x * 2. * (ec2 * up + eb2 * dn + eb2 * str) + f2charm + f2beauty;
/*<       F2 = F2/alfa >*/
    *f2 /= .00729735308;
/*<  2    continue >*/
L2:
/*<       end >*/
    return 0;
} /* param_ */

/* ******************************************************************** */
/* ******************************************************************** */
/*<       SUBROUTINE POINTLIKE(OPT,x,Q2,glpl,uppl,dnpl,chpl,btpl) >*/
/* Subroutine */ int pointlike_(integer *opt, doublereal *x, doublereal *q2, 
	doublereal *glpl, doublereal *uppl, doublereal *dnpl, doublereal *
	chpl, doublereal *btpl)
{
    /* Initialized data */

    static doublereal pargl[19] = { -.43865,2.7174,.36752,.086893,.010556,
	    -.099005,1.0648,3.6717,2.1944,.236795,-.19994,-.34992,.049525,
	    .3483,.14342,2.5071,1.9358,-.11849,.028124 };
    static doublereal parup[19] = { -1.0711,3.132,.69243,-.058266,.0097377,
	    -.0068345,.22297,6.4289,1.7302,.8794,2.6878,.20506,-.10617,.15211,
	    .013567,2.2802,.76997,-.110241,-.040252 };
    static doublereal pardn[19] = { -1.1357,3.1187,.6629,.098814,-.092892,
	    -.006614,-.31385,6.4671,1.6996,11.777,-11.124,-.0673,.049949,
	    .020427,-.0037558,2.2834,.84262,.03476,-.20135 };
    static doublereal parch1[21] = { 2.9808,28.682,2.4863,-.18826,.18508,
	    -.0014153,-.48961,.20911,2.7644,-7.6307,394.58,.13565,-.11764,
	    -.01151,.1881,-2.8544,.93717,5.6807,-541.82,14.256,200.82 };
    static doublereal parch2[21] = { -1.8095,7.9399,.041563,-.54831,.19484,
	    -.39046,.12717,8.7191,4.2616,-.30307,7.2383,.33412,.041562,.37194,
	    .05928,3.0194,.73993,.2943,-1.5995,0.,0. };
    static doublereal parbt1[21] = { 2.2849,6.0408,-.11577,-.26971,.27033,
	    .0022862,.30807,14.812,1.7148,3.814,2.2292,.17942,-.18358,
	    -.0016837,-.1049,-1.2977,2.3532,-1.0514,20.194,.0061059,.053734 };
    static doublereal parbt2[21] = { -5.0607,16.59,.8719,-.7279,-.62903,
	    -2.4467,.56575,1.4687,1.1706,-.084651,9.6036,.36549,.56817,1.6783,
	    -.1912,9.6071,.99674,-.083206,-3.4864,0.,0. };

    extern doublereal pl_(doublereal *, doublereal *, doublereal *), bpl_(
	    doublereal *, doublereal *, doublereal *), cpl_(doublereal *, 
	    doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension pargl(19),parup(19),pardn(19) >*/
/*<       dimension parch1(21),parch2(21),parbt1(21),parbt2(21) >*/
/*<       integer OPT >*/
/* ***  gluon       **** */
/*<        >*/
/* ***  up          **** */
/*<        >*/
/* ***  down = str  **** */
/*<        >*/
/* ***  charm  **** */
/*<        >*/
/*<        >*/
/* ***  bottom  **** */
/*<        >*/
/*<        >*/
/*<       glpl = pl(x,Q2,pargl) >*/
    *glpl = pl_(x, q2, pargl);
/*<       uppl = pl(x,Q2,parup) >*/
    *uppl = pl_(x, q2, parup);
/*<       dnpl = pl(x,Q2,pardn) >*/
    *dnpl = pl_(x, q2, pardn);
/*<       if (OPT.EQ.1) goto 1 >*/
    if (*opt == 1) {
	goto L1;
    }
/*<       if (Q2.LE.10.d0) then >*/
    if (*q2 <= 10.) {
/*<          chpl = cpl(x,Q2,parch1) >*/
	*chpl = cpl_(x, q2, parch1);
/*<       else >*/
    } else {
/*<          chpl = cpl(x,Q2,parch2) >*/
	*chpl = cpl_(x, q2, parch2);
/*<       endif >*/
    }
/*<       if (Q2.LE.100.d0) then >*/
    if (*q2 <= 100.) {
/*<          btpl = bpl(x,Q2,parbt1) >*/
	*btpl = bpl_(x, q2, parbt1);
/*<       else >*/
    } else {
/*<          btpl = bpl(x,Q2,parbt2) >*/
	*btpl = bpl_(x, q2, parbt2);
/*<       endif >*/
    }
/*<  1    continue >*/
L1:
/*<       END >*/
    return 0;
} /* pointlike_ */

/* ********************* */
/* ********************* */
/*<       SUBROUTINE HADRONLIKE(OPT,x,Q2,glhd,vlhd,sthd,chhd,bthd) >*/
/* Subroutine */ int hadronlike_(integer *opt, doublereal *x, doublereal *q2, 
	doublereal *glhd, doublereal *vlhd, doublereal *sthd, doublereal *
	chhd, doublereal *bthd)
{
    /* Initialized data */

    static doublereal pargl[18] = { .59945,1.1285,-.19898,1.9942,-1.9848,
	    -.34948,1.0012,1.2287,4.923,.21294,.57414,-1.8306,1.4136,.47058,
	    .99767,2.4447,.18526,2.745 };
    static doublereal parvl[10] = { 1.0898,.78391,.42654,-1.6576,.96155,
	    .38087,-.06872,-1.2128,1.7075,1.8441 };
    static doublereal parst[14] = { .7166,.72289,.60478,4.2106,4.1494,4.5179,
	    5.2812,1.0497,-.21562,.03616,-.85835,.34866,1.9219,-.152 };
    static doublereal parch1[17] = { 5.6729,1.6248,-2586.4,2695.,1.5146,
	    -3.9185,3.6126,1.4575,-.70433,1910.1,-1688.2,3.1028,11.738,
	    -1.0291,0.,0.,0. };
    static doublereal parch2[17] = { -1.647,-.78809,-2.0561,2.1266,3.0301,
	    4.1282,.89599,.72738,.90278,.75576,.66383,-1.7499,1.6929,1.2761,
	    -.15061,-.26292,1.6466 };
    static doublereal parbt1[16] = { -10.21,.82278,-99.613,492.61,3.3917,
	    5.6829,-2.0137,-2.2296,.081818,171.25,-420.45,.084256,-.23571,
	    4.6955,0.,0. };
    static doublereal parbt2[16] = { 2.4198,-.98933,-2.1109,9.0196,3.6455,
	    4.6196,.66454,.40703,.42366,1.2711,-3.6082,-4.1353,2.4212,1.1109,
	    .15817,2.3615 };

    extern doublereal chm_(doublereal *, doublereal *, doublereal *), val_(
	    doublereal *, doublereal *, doublereal *), bot_(doublereal *, 
	    doublereal *, doublereal *), glu_(doublereal *, doublereal *, 
	    doublereal *), str_(doublereal *, doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension pargl(18),parvl(10),parst(14) >*/
/*<       dimension parch1(17),parch2(17),parbt1(16),parbt2(16) >*/
/*<       integer OPT >*/
/* ***  gluon       **** */
/*<        >*/
/* ***  valence     **** */
/*<        >*/
/* ***  strange     **** */
/*<        >*/
/* ***  charm  **** */
/*<        >*/
/*<        >*/
/* ***  bottom  **** */
/*<        >*/
/*<        >*/
/*<       glhd = glu(x,Q2,pargl) >*/
    *glhd = glu_(x, q2, pargl);
/*<       vlhd = val(x,Q2,parvl) >*/
    *vlhd = val_(x, q2, parvl);
/*<       sthd = str(x,Q2,parst) >*/
    *sthd = str_(x, q2, parst);
/*<       if (OPT.EQ.1) goto 1 >*/
    if (*opt == 1) {
	goto L1;
    }
/*<       if (Q2.LE.10.d0) then >*/
    if (*q2 <= 10.) {
/*<          chhd = chm(x,Q2,parch1) >*/
	*chhd = chm_(x, q2, parch1);
/*<       else >*/
    } else {
/*<          chhd = chm(x,Q2,parch2) >*/
	*chhd = chm_(x, q2, parch2);
/*<       endif >*/
    }
/*<       if (Q2.LE.100.d0) then >*/
    if (*q2 <= 100.) {
/*<          bthd = bot(x,Q2,parbt1) >*/
	*bthd = bot_(x, q2, parbt1);
/*<       else >*/
    } else {
/*<          bthd = bot(x,Q2,parbt2) >*/
	*bthd = bot_(x, q2, parbt2);
/*<       endif >*/
    }
/*<  1    continue >*/
L1:
/*<       END >*/
    return 0;
} /* hadronlike_ */

/* ******************************************************************** */
/*<       SUBROUTINE GLUON(x,Q2,glun) >*/
/* Subroutine */ int gluon_(doublereal *x, doublereal *q2, doublereal *glun)
{
    /* Initialized data */

    static doublereal plglu[19] = { -.4387,2.717,.3675,.08689,.01056,-.099,
	    1.065,3.672,2.194,.2368,-.1999,-.3499,.04953,.3483,.1434,2.507,
	    1.936,-.1185,.02812 };
    static doublereal hadglu[18] = { .599449,1.12849,-.19898,1.99417,-1.9848,
	    -.349479,1.0012,1.22871,4.92304,.212938,.57414,-1.8306,1.41359,
	    .470584,.997665,2.4447,.185258,2.745 };

    extern doublereal pl_(doublereal *, doublereal *, doublereal *), glu_(
	    doublereal *, doublereal *, doublereal *);
    static doublereal glhd, glpl;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       Parameter (alfa=7.29735308D-3) >*/
/*<       dimension plglu(19),hadglu(18) >*/
/* ***  point-like       **** */
/*<        >*/
/* ***  hadron-like      **** */
/*<        >*/
/*<       glpl = pl(x,Q2,plglu) >*/
    glpl = pl_(x, q2, plglu);
/*<       glhd = glu(x,Q2,hadglu) >*/
    glhd = glu_(x, q2, hadglu);
/*<       glun = alfa*(glpl + glhd)/x >*/
    *glun = (glpl + glhd) * .00729735308 / *x;
/*<       END >*/
    return 0;
} /* gluon_ */

/* ******************************************************************* */
/* ***                      POINT-LIKE                            **** */
/* ******************************************************************* */
/*<       DOUBLE PRECISION FUNCTION pl(x,Q2,PAR) >*/
doublereal pl_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__, e, s, as, ep, bs, dlg, lam2, alfa, beta,
	     alfap;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(19) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa  = PAR(1) >*/
    alfa = par[1];
/*<       alfap = PAR(2) >*/
    alfap = par[2];
/*<       beta  = PAR(3) >*/
    beta = par[3];
/*<       A =  PAR(4)  + PAR(12)*s >*/
    a = par[4] + par[12] * s;
/*<       B =  PAR(5)  + PAR(13)*s >*/
    b = par[5] + par[13] * s;
/*<       C =  PAR(6)  + PAR(14)*s >*/
    c__ = par[6] + par[14] * s;
/*<       D =  PAR(7)  + PAR(15)*s >*/
    d__ = par[7] + par[15] * s;
/*<       E =  PAR(8)  + PAR(16)*s >*/
    e = par[8] + par[16] * s;
/*<       EP = PAR(9)  + PAR(17)*s >*/
    ep = par[9] + par[17] * s;
/*<       AS = PAR(10) + PAR(18)*s >*/
    as = par[10] + par[18] * s;
/*<       BS = PAR(11) + PAR(19)*s >*/
    bs = par[11] + par[19] * s;
/*<       pl = s**alfa*x**AS*(A+B*DSQRT(x)+C*x**BS) >*/
    ret_val = pow_dd(&s, &alfa) * pow_dd(x, &as) * (a + b * sqrt(*x) + c__ * 
	    pow_dd(x, &bs));
/*<       pl = pl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg)) >*/
    ret_val += pow_dd(&s, &alfap) * exp(-e + sqrt(ep * pow_dd(&s, &beta) * 
	    dlg));
/*<       pl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*pl*(1.d0-x)**D >*/
    d__1 = 1. - *x;
    ret_val = log(*q2 / lam2) * .71619724391352901 * ret_val * pow_dd(&d__1, &
	    d__);
/*<       END >*/
    return ret_val;
} /* pl_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION cpl(x,Q2,PAR) >*/
doublereal cpl_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__, e, s, y, s2, as, ep, bs, dlg, lam2, 
	    alfa, beta, alfap;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(21) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       s2 = s*s >*/
    s2 = s * s;
/*<       cpl = 0.d0 >*/
    ret_val = 0.;
/*<       y = x + 1.d0 - Q2/(Q2+6.76d0) >*/
    y = *x + 1. - *q2 / (*q2 + 6.76);
/*<       if (y.GE.1.d0) goto 10 >*/
    if (y >= 1.) {
	goto L10;
    }
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa  = PAR(1) >*/
    alfa = par[1];
/*<       alfap = PAR(2) >*/
    alfap = par[2];
/*<       beta  = PAR(3) >*/
    beta = par[3];
/*<       A =  PAR(4)  + PAR(12)*s >*/
    a = par[4] + par[12] * s;
/*<       B =  PAR(5)  + PAR(13)*s >*/
    b = par[5] + par[13] * s;
/*<       C =  PAR(6)  + PAR(14)*s >*/
    c__ = par[6] + par[14] * s;
/*<       D =  PAR(7)  + PAR(15)*s >*/
    d__ = par[7] + par[15] * s;
/*<       E =  PAR(8)  + PAR(16)*s + PAR(20)*s2 >*/
    e = par[8] + par[16] * s + par[20] * s2;
/*<       EP = PAR(9)  + PAR(17)*s >*/
    ep = par[9] + par[17] * s;
/*<       AS = PAR(10) + PAR(18)*s >*/
    as = par[10] + par[18] * s;
/*<       BS = PAR(11) + PAR(19)*s + PAR(21)*s2 >*/
    bs = par[11] + par[19] * s + par[21] * s2;
/*<       cpl = s**alfa*y**AS*(A+B*DSQRT(y)+C*y**BS) >*/
    ret_val = pow_dd(&s, &alfa) * pow_dd(&y, &as) * (a + b * sqrt(y) + c__ * 
	    pow_dd(&y, &bs));
/*<       cpl = cpl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg)) >*/
    ret_val += pow_dd(&s, &alfap) * exp(-e + sqrt(ep * pow_dd(&s, &beta) * 
	    dlg));
/*<       cpl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*cpl*(1.d0-y)**D >*/
    d__1 = 1. - y;
    ret_val = log(*q2 / lam2) * .71619724391352901 * ret_val * pow_dd(&d__1, &
	    d__);
/*<  10   continue >*/
L10:
/*<       END >*/
    return ret_val;
} /* cpl_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION bpl(x,Q2,PAR) >*/
doublereal bpl_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal a, b, c__, d__, e, s, y, s2, as, ep, bs, ds, dlg, lam2, 
	    alfa, beta, alfap;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(21) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       s2 = s*s >*/
    s2 = s * s;
/*<       ds = DSQRT(s) >*/
    ds = sqrt(s);
/*<       bpl = 0.d0 >*/
    ret_val = 0.;
/*<       y = x + 1.d0 - Q2/(Q2+73.96d0) >*/
    y = *x + 1. - *q2 / (*q2 + 73.96);
/*<       if (y.GE.1.d0) goto 10 >*/
    if (y >= 1.) {
	goto L10;
    }
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa  = PAR(1) >*/
    alfa = par[1];
/*<       alfap = PAR(2) >*/
    alfap = par[2];
/*<       beta  = PAR(3) >*/
    beta = par[3];
/*<       A =  PAR(4)  + PAR(12)*s >*/
    a = par[4] + par[12] * s;
/*<       B =  PAR(5)  + PAR(13)*s + PAR(20)*s2 >*/
    b = par[5] + par[13] * s + par[20] * s2;
/*<       C =  PAR(6)  + PAR(14)*s >*/
    c__ = par[6] + par[14] * s;
/*<       D =  PAR(7)  + PAR(15)*s >*/
    d__ = par[7] + par[15] * s;
/*<       E =  PAR(8)  + PAR(16)*s >*/
    e = par[8] + par[16] * s;
/*<       EP = PAR(9)  + PAR(17)*s + PAR(21)*ds >*/
    ep = par[9] + par[17] * s + par[21] * ds;
/*<       AS = PAR(10) + PAR(18)*s >*/
    as = par[10] + par[18] * s;
/*<       BS = PAR(11) + PAR(19)*s >*/
    bs = par[11] + par[19] * s;
/*<       bpl = s**alfa*y**AS*(A+B*DSQRT(y)+C*y**BS) >*/
    ret_val = pow_dd(&s, &alfa) * pow_dd(&y, &as) * (a + b * sqrt(y) + c__ * 
	    pow_dd(&y, &bs));
/*<       bpl = bpl + s**alfap*DEXP(-E + DSQRT(EP*s**beta*dlg)) >*/
    ret_val += pow_dd(&s, &alfap) * exp(-e + sqrt(ep * pow_dd(&s, &beta) * 
	    dlg));
/*<       bpl = 9.d0/(4.d0*Pi)*DLOG(Q2/Lam2)*bpl*(1.d0-y)**D >*/
    d__1 = 1. - y;
    ret_val = log(*q2 / lam2) * .71619724391352901 * ret_val * pow_dd(&d__1, &
	    d__);
/*<  10   continue >*/
L10:
/*<       END >*/
    return ret_val;
} /* bpl_ */

/* ******************************************************************* */
/* ***                     HADRON-LIKE                            **** */
/* ******************************************************************* */
/*<       DOUBLE PRECISION FUNCTION val(x,Q2,PAR) >*/
doublereal val_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *);

    /* Local variables */
    static doublereal c__, d__, s, ac, bc, as, lam2;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(10) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       AC = PAR(1) + PAR(6)*s >*/
    ac = par[1] + par[6] * s;
/*<       AS = PAR(2) + PAR(7)*s >*/
    as = par[2] + par[7] * s;
/*<       BC = PAR(3) + PAR(8)*s >*/
    bc = par[3] + par[8] * s;
/*<       C  = PAR(4) + PAR(9)*s >*/
    c__ = par[4] + par[9] * s;
/*<       D  = PAR(5) + PAR(10)*s >*/
    d__ = par[5] + par[10] * s;
/*<       val = AC*x**AS*(1.d0+BC*DSQRT(x)+C*x) >*/
    ret_val = ac * pow_dd(x, &as) * (bc * sqrt(*x) + 1. + c__ * *x);
/*<       val = val*(1.d0-x)**D >*/
    d__1 = 1. - *x;
    ret_val *= pow_dd(&d__1, &d__);
/*<       END >*/
    return ret_val;
} /* val_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION glu(x,Q2,PAR) >*/
doublereal glu_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal c__, d__, e, s, ac, bc, as, bs, ep, dlg, lam2, alfa, 
	    beta;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(18) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa = PAR(1) >*/
    alfa = par[1];
/*<       beta = PAR(2) >*/
    beta = par[2];
/*<       AC = PAR(3)  + PAR(11)*s >*/
    ac = par[3] + par[11] * s;
/*<       BC = PAR(4)  + PAR(12)*s >*/
    bc = par[4] + par[12] * s;
/*<       C  = PAR(5)  + PAR(13)*s >*/
    c__ = par[5] + par[13] * s;
/*<       AS = PAR(6)  + PAR(14)*s >*/
    as = par[6] + par[14] * s;
/*<       BS = PAR(7)  + PAR(15)*s >*/
    bs = par[7] + par[15] * s;
/*<       E  = PAR(8)  + PAR(16)*s >*/
    e = par[8] + par[16] * s;
/*<       EP = PAR(9)  + PAR(17)*s >*/
    ep = par[9] + par[17] * s;
/*<       D  = PAR(10) + PAR(18)*s >*/
    d__ = par[10] + par[18] * s;
/*<       glu = x**AS*(AC+BC*DSQRT(x)+C*x) >*/
    ret_val = pow_dd(x, &as) * (ac + bc * sqrt(*x) + c__ * *x);
/*<       glu = glu + s**alfa*DEXP(-E+DSQRT(EP*s**beta*dlg)) >*/
    ret_val += pow_dd(&s, &alfa) * exp(-e + sqrt(ep * pow_dd(&s, &beta) * dlg)
	    );
/*<       glu = glu*(1.d0-x)**D >*/
    d__1 = 1. - *x;
    ret_val *= pow_dd(&d__1, &d__);
/*<       END >*/
    return ret_val;
} /* glu_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION str(x,Q2,PAR) >*/
doublereal str_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal d__, e, s, ac, bc, as, ep, dlg, lam2, alfa, beta;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(14) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa = PAR(1) >*/
    alfa = par[1];
/*<       AS = PAR(2)  + PAR(9)*s >*/
    as = par[2] + par[9] * s;
/*<       AC = PAR(3)  + PAR(10)*s >*/
    ac = par[3] + par[10] * s;
/*<       BC = PAR(4)  + PAR(11)*s >*/
    bc = par[4] + par[11] * s;
/*<       D  = PAR(5)  + PAR(12)*s >*/
    d__ = par[5] + par[12] * s;
/*<       E  = PAR(6)  + PAR(13)*s >*/
    e = par[6] + par[13] * s;
/*<       EP = PAR(7)  + PAR(14)*s >*/
    ep = par[7] + par[14] * s;
/*<       beta = PAR(8) >*/
    beta = par[8];
/*<       str = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(x)+BC*x) >*/
    ret_val = pow_dd(&s, &alfa) / pow_dd(&dlg, &as) * (ac * sqrt(*x) + 1. + 
	    bc * *x);
/*<       str = str*(1.d0-x)**D*DEXP(-E+DSQRT(EP*s**beta*dlg)) >*/
    d__1 = 1. - *x;
    ret_val = ret_val * pow_dd(&d__1, &d__) * exp(-e + sqrt(ep * pow_dd(&s, &
	    beta) * dlg));
/*<       END >*/
    return ret_val;
} /* str_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION chm(x,Q2,PAR) >*/
doublereal chm_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal d__, e, s, y, s2, ac, bc, as, ep, dlg, lam2, alfa, beta;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(17) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       s2 = s*s >*/
    s2 = s * s;
/*<       chm = 0.d0 >*/
    ret_val = 0.;
/*<       y = x + 1.d0 - Q2/(Q2+6.76d0) >*/
    y = *x + 1. - *q2 / (*q2 + 6.76);
/*<       if (y.GE.1.d0) goto 10 >*/
    if (y >= 1.) {
	goto L10;
    }
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa = PAR(1) >*/
    alfa = par[1];
/*<       AS = PAR(2)  + PAR(9)*s >*/
    as = par[2] + par[9] * s;
/*<       AC = PAR(3)  + PAR(10)*s >*/
    ac = par[3] + par[10] * s;
/*<       BC = PAR(4)  + PAR(11)*s >*/
    bc = par[4] + par[11] * s;
/*<       D  = PAR(5)  + PAR(12)*s + PAR(17)*s2 >*/
    d__ = par[5] + par[12] * s + par[17] * s2;
/*<       E  = PAR(6)  + PAR(13)*s + PAR(16)*s2 >*/
    e = par[6] + par[13] * s + par[16] * s2;
/*<       EP = PAR(7)  + PAR(14)*s + PAR(15)*s2 >*/
    ep = par[7] + par[14] * s + par[15] * s2;
/*<       beta = PAR(8) >*/
    beta = par[8];
/*<       chm = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(y)+BC*y) >*/
    ret_val = pow_dd(&s, &alfa) / pow_dd(&dlg, &as) * (ac * sqrt(y) + 1. + bc 
	    * y);
/*<       chm = chm*(1.d0-y)**D*DEXP(-E+EP*DSQRT(s**beta*dlg)) >*/
    d__1 = 1. - y;
    ret_val = ret_val * pow_dd(&d__1, &d__) * exp(-e + ep * sqrt(pow_dd(&s, &
	    beta) * dlg));
/*<  10   continue >*/
L10:
/*<       END >*/
    return ret_val;
} /* chm_ */

/* ********* */
/* ********* */
/*<       DOUBLE PRECISION FUNCTION bot(x,Q2,PAR) >*/
doublereal bot_(doublereal *x, doublereal *q2, doublereal *par)
{
    /* System generated locals */
    doublereal ret_val, d__1;

    /* Builtin functions */
    double log(doublereal), sqrt(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static doublereal d__, e, s, y, s2, ac, bc, as, ep, dlg, lam2, alfa, beta;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension PAR(16) >*/
/*<       Parameter (Pi=3.1415926535897932d0) >*/
/*<       Lam2 = 0.221d0*0.221d0 >*/
    /* Parameter adjustments */
    --par;

    /* Function Body */
    lam2 = .048841000000000002;
/*<       s = DLOG(DLOG(Q2/Lam2)/DLOG(0.25d0/Lam2)) >*/
    s = log(log(*q2 / lam2) / log(.25 / lam2));
/*<       s2 = s*s >*/
    s2 = s * s;
/*<       bot = 0.d0 >*/
    ret_val = 0.;
/*<       y = x + 1.d0 - Q2/(Q2+73.96d0) >*/
    y = *x + 1. - *q2 / (*q2 + 73.96);
/*<       if (y.GE.1.d0) goto 10 >*/
    if (y >= 1.) {
	goto L10;
    }
/*<       dlg = DLOG(1.d0/x) >*/
    dlg = log(1. / *x);
/*<       alfa = PAR(1) >*/
    alfa = par[1];
/*<       AS = PAR(2)  + PAR(9)*s  + PAR(15)*s2 >*/
    as = par[2] + par[9] * s + par[15] * s2;
/*<       AC = PAR(3)  + PAR(10)*s >*/
    ac = par[3] + par[10] * s;
/*<       BC = PAR(4)  + PAR(11)*s >*/
    bc = par[4] + par[11] * s;
/*<       D  = PAR(5)  + PAR(12)*s + PAR(16)*s2 >*/
    d__ = par[5] + par[12] * s + par[16] * s2;
/*<       E  = PAR(6)  + PAR(13)*s >*/
    e = par[6] + par[13] * s;
/*<       EP = PAR(7)  + PAR(14)*s >*/
    ep = par[7] + par[14] * s;
/*<       beta = PAR(8) >*/
    beta = par[8];
/*<       bot = s**alfa/(dlg**AS)*(1.d0+AC*DSQRT(y)+BC*y) >*/
    ret_val = pow_dd(&s, &alfa) / pow_dd(&dlg, &as) * (ac * sqrt(y) + 1. + bc 
	    * y);
/*<       bot = bot*(1.d0-y)**D*DEXP(-E+EP*DSQRT(s**beta*dlg)) >*/
    d__1 = 1. - y;
    ret_val = ret_val * pow_dd(&d__1, &d__) * exp(-e + ep * sqrt(pow_dd(&s, &
	    beta) * dlg));
/*<  10   continue >*/
L10:
/*<       END >*/
    return ret_val;
} /* bot_ */

/* ************************************************************************ */
/* ***                    Running alpha strong                         **** */
/* ************************************************************************ */
/*<       DOUBLE PRECISION FUNCTION alfas(Q2,Nf,LAM) >*/
doublereal alfas_(doublereal *q2, integer *nf, doublereal *lam)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double log(doublereal);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       PARAMETER(PI=3.1415926535898d0) >*/
/*<       INTEGER Nf >*/
/*<       alfas = 12.0d0*PI/((33.0d0-2.0d0*Nf)*DLOG(Q2/(LAM*LAM)) ) >*/
    ret_val = 37.699111843077603 / ((33. - *nf * 2.) * log(*q2 / (*lam * *lam)
	    ));
/*<       END >*/
    return ret_val;
} /* alfas_ */

/* *************************************************************************** */
/*<       subroutine intxr(x0,step,x,Q2,a,b,n,eps,result,t) >*/
/* Subroutine */ int intxr_(doublereal *x0, integer *step, doublereal *x, 
	doublereal *q2, doublereal *a, doublereal *b, integer *n, doublereal *
	eps, doublereal *result, logical *t)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Local variables */
    static doublereal c__;
    static integer i__;
    static doublereal aa[60], bb[60], ep, ra1[2], ra2[2], rb1[2], rb2[2];
    static integer ind;
    static doublereal roz, eps1[2], eps2[2];
    static integer iadr[60];
    static doublereal rozn, reslt2[2];
    extern /* Subroutine */ int gauscxr_(doublereal *, integer *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/*<       implicit double precision (a-z) >*/
/*<       integer *4 iadr(60) >*/
/*<       integer n,I,ind,step >*/
/*<        >*/
/*<       logical t >*/
/*<       t=.true. >*/
    /* Parameter adjustments */
    --result;

    /* Function Body */
    *t = TRUE_;
/*<       DO 110 I=1,2 >*/
    for (i__ = 1; i__ <= 2; ++i__) {
/*<          result(I)=0.d0 >*/
	result[i__] = 0.;
/*< 110      reslt2(I)=0.d0 >*/
/* L110: */
	reslt2[i__ - 1] = 0.;
    }
/*<       ind=1 >*/
    ind = 1;
/*<       iadr(1)=-1 >*/
    iadr[0] = -1;
/*<       aa(1)=a >*/
    aa[0] = *a;
/*<       bb(1)=b >*/
    bb[0] = *b;
/*< 1     c=(aa(ind)+bb(ind))/2.0d0 >*/
L1:
    c__ = (aa[ind - 1] + bb[ind - 1]) / 2.;
/*<       call gauscxr(x0,step,x,Q2,aa(ind),c,ra1,ra2) >*/
    gauscxr_(x0, step, x, q2, &aa[ind - 1], &c__, ra1, ra2);
/*<       call gauscxr(x0,step,x,Q2,c,bb(ind),rb1,rb2) >*/
    gauscxr_(x0, step, x, q2, &c__, &bb[ind - 1], rb1, rb2);
/*<       DO 200 I=1,2 >*/
    for (i__ = 1; i__ <= 2; ++i__) {
/*<         eps1(I)=dabs(ra1(I)-ra2(I))/(dabs(ra1(I)+result(I))+1.0d-300) >*/
	eps1[i__ - 1] = (d__2 = ra1[i__ - 1] - ra2[i__ - 1], abs(d__2)) / ((
		d__1 = ra1[i__ - 1] + result[i__], abs(d__1)) + 1e-300);
/*< 200     eps2(I)=dabs(rb1(I)-rb2(I))/(dabs(rb1(I)+result(I))+1.0d-300) >*/
/* L200: */
	eps2[i__ - 1] = (d__2 = rb1[i__ - 1] - rb2[i__ - 1], abs(d__2)) / ((
		d__1 = rb1[i__ - 1] + result[i__], abs(d__1)) + 1e-300);
    }
/*<       rozn = eps1(1) - eps2(1) >*/
    rozn = eps1[0] - eps2[0];
/*<       DO 300 I=2,2 >*/
    for (i__ = 2; i__ <= 2; ++i__) {
/*<          roz = eps1(I)-eps2(I) >*/
	roz = eps1[i__ - 1] - eps2[i__ - 1];
/*<          if (roz.GT.rozn) rozn = roz >*/
	if (roz > rozn) {
	    rozn = roz;
	}
/*< 300   CONTINUE >*/
/* L300: */
    }
/*<       if(rozn) 10,10,20 >*/
    if (rozn <= 0.) {
	goto L10;
    } else {
	goto L20;
    }
/*< 10    rozn = eps1(1) - eps >*/
L10:
    rozn = eps1[0] - *eps;
/*<       DO 400 I=2,2 >*/
    for (i__ = 2; i__ <= 2; ++i__) {
/*<          roz = eps1(I) - eps >*/
	roz = eps1[i__ - 1] - *eps;
/*<          if (roz.GT.rozn) rozn = roz >*/
	if (roz > rozn) {
	    rozn = roz;
	}
/*< 400   CONTINUE >*/
/* L400: */
    }
/*<       if(rozn) 12,12,11 >*/
    if (rozn <= 0.) {
	goto L12;
    } else {
	goto L11;
    }
/*< 11    if(ind-n) 13,15,15 >*/
L11:
    if (ind - *n >= 0) {
	goto L15;
    } else {
	goto L13;
    }
/*< 15    t=.false. >*/
L15:
    *t = FALSE_;
/*< 12    DO 500 I=1,2 >*/
L12:
    for (i__ = 1; i__ <= 2; ++i__) {
/*<          result(I)=result(I)+ra1(I) >*/
	result[i__] += ra1[i__ - 1];
/*< 500      reslt2(I)=reslt2(I)+ra2(I) >*/
/* L500: */
	reslt2[i__ - 1] += ra2[i__ - 1];
    }
/*<       iadr(ind)=iadr(ind)+100 >*/
    iadr[ind - 1] += 100;
/*<       if(iadr(ind)-150) 20,20,30 >*/
    if (iadr[ind - 1] - 150 <= 0) {
	goto L20;
    } else {
	goto L30;
    }
/*< 13    ind=ind+1 >*/
L13:
    ++ind;
/*<       iadr(ind)=0 >*/
    iadr[ind - 1] = 0;
/*<       aa(ind)=aa(ind-1) >*/
    aa[ind - 1] = aa[ind - 2];
/*<       bb(ind)=(aa(ind-1)+bb(ind-1))/2. >*/
    bb[ind - 1] = (aa[ind - 2] + bb[ind - 2]) / (float)2.;
/*<       go to 1 >*/
    goto L1;
/*< 14    iadr(ind)=iadr(ind)+100 >*/
L14:
    iadr[ind - 1] += 100;
/*<       if(iadr(ind)-150) 23,23,30 >*/
    if (iadr[ind - 1] - 150 <= 0) {
	goto L23;
    } else {
	goto L30;
    }
/*< 20    rozn = eps2(1) - eps >*/
L20:
    rozn = eps2[0] - *eps;
/*<       DO 600 I=2,2 >*/
    for (i__ = 2; i__ <= 2; ++i__) {
/*<          roz = eps2(I) - eps >*/
	roz = eps2[i__ - 1] - *eps;
/*<          if (roz.GT.rozn) rozn = roz >*/
	if (roz > rozn) {
	    rozn = roz;
	}
/*< 600   CONTINUE >*/
/* L600: */
    }
/*<       if(rozn) 22,22,21 >*/
    if (rozn <= 0.) {
	goto L22;
    } else {
	goto L21;
    }
/*< 21    if(ind-n) 23,25,25 >*/
L21:
    if (ind - *n >= 0) {
	goto L25;
    } else {
	goto L23;
    }
/*< 25    t =.false. >*/
L25:
    *t = FALSE_;
/*< 22    DO 700 I=1,2 >*/
L22:
    for (i__ = 1; i__ <= 2; ++i__) {
/*<          result(I)=result(I)+rb1(I) >*/
	result[i__] += rb1[i__ - 1];
/*< 700      reslt2(I)=reslt2(I)+rb2(I) >*/
/* L700: */
	reslt2[i__ - 1] += rb2[i__ - 1];
    }
/*<       iadr(ind)=iadr(ind)+100 >*/
    iadr[ind - 1] += 100;
/*<       if(iadr(ind)-150) 10,10,30 >*/
    if (iadr[ind - 1] - 150 <= 0) {
	goto L10;
    } else {
	goto L30;
    }
/*< 23    ind=ind+1 >*/
L23:
    ++ind;
/*<       iadr(ind)=1 >*/
    iadr[ind - 1] = 1;
/*<       aa(ind)=(aa(ind-1)+bb(ind-1))/2. >*/
    aa[ind - 1] = (aa[ind - 2] + bb[ind - 2]) / (float)2.;
/*<       bb(ind)=bb(ind-1) >*/
    bb[ind - 1] = bb[ind - 2];
/*<       go to 1 >*/
    goto L1;
/*< 24    iadr(ind)=iadr(ind)+100 >*/
L24:
    iadr[ind - 1] += 100;
/*<       if(iadr(ind)-150) 13,13,30 >*/
    if (iadr[ind - 1] - 150 <= 0) {
	goto L13;
    } else {
	goto L30;
    }
/*< 30    ind=ind-1 >*/
L30:
    --ind;
/*<       if(iadr(ind+1)-200) 100,14,24 >*/
    if ((i__1 = iadr[ind] - 200) < 0) {
	goto L100;
    } else if (i__1 == 0) {
	goto L14;
    } else {
	goto L24;
    }
/*< 100   eps = dabs(result(1)-reslt2(1))/(dabs(result(1))+1.d-300) >*/
L100:
    *eps = (d__1 = result[1] - reslt2[0], abs(d__1)) / (abs(result[1]) + 
	    1e-300);
/*<       DO 800 I=2,2 >*/
    for (i__ = 2; i__ <= 2; ++i__) {
/*<          ep = dabs(result(I)-reslt2(I))/(dabs(result(I)+1.d-300)) >*/
	ep = (d__1 = result[i__] - reslt2[i__ - 1], abs(d__1)) / (d__2 = 
		result[i__] + 1e-300, abs(d__2));
/*<          if (ep.GT.eps) eps = ep >*/
	if (ep > *eps) {
	    *eps = ep;
	}
/*< 800   CONTINUE >*/
/* L800: */
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* intxr_ */

/* *************************************************************************** */
/*<       subroutine gauscxr(x0,step,x,Q2,a,b,gauskr,gaus) >*/
/* Subroutine */ int gauscxr_(doublereal *x0, integer *step, doublereal *x, 
	doublereal *q2, doublereal *a, doublereal *b, doublereal *gauskr, 
	doublereal *gaus)
{
    /* Initialized data */

    static doublereal g[24]	/* was [3][8] */ = { .993379875881716,0.,
	    .0178223833207104,.960289856497536,.101228536290376,
	    .0494393950021394,.894120906847456,0.,.0824822989313584,
	    .796666477413626,.222381034453374,.11164637082684,
	    .672354070945158,0.,.136263109255172,.525532409916329,
	    .313706645877887,.156652606168188,.360701097928132,0.,
	    .172070608555211,.18343464249565,.362683783378362,
	    .181400025068035 };
    static doublereal g39 = .184446405744692;

    /* System generated locals */
    doublereal d__1;

    /* Local variables */
    static doublereal c__[2], d__, e;
    static integer i__, l;
    static doublereal y, z1, z2, z3;
    extern doublereal pqg_(doublereal *), fun_(doublereal *, doublereal *);
    static doublereal glu1, glu2, res1[2], res2[2], res3[2], glu3;
    extern /* Subroutine */ int gluon_(doublereal *, doublereal *, doublereal 
	    *);

/*<       implicit double precision(a-z) >*/
/*<       parameter (nvar = 12) >*/
/*<       integer l,I,step >*/
/*<        >*/
/*<       Parameter (Pi=3.1415926535897932d0, alfa=7.29735308D-3) >*/
/*<        >*/
    /* Parameter adjustments */
    --gaus;
    --gauskr;

    /* Function Body */
/*<       data g39/1.844464 0574 4692d-1/ >*/
/*<       do 10 I=1,2 >*/
    for (i__ = 1; i__ <= 2; ++i__) {
/*<         gaus(I)=0.d0 >*/
	gaus[i__] = 0.;
/*<  10     gauskr(I)=0.d0 >*/
/* L10: */
	gauskr[i__] = 0.;
    }
/*<       d=(b-a)/2. >*/
    d__ = (*b - *a) / (float)2.;
/*<       e=(b+a)/2. >*/
    e = (*b + *a) / (float)2.;
/*<       z3 = e >*/
    z3 = e;
/*<       call GLUON(z3,Q2,glu3) >*/
    gluon_(&z3, q2, &glu3);
/*<       res3(1) = glu3*fun(x0/z3,Q2) >*/
    d__1 = *x0 / z3;
    res3[0] = glu3 * fun_(&d__1, q2);
/*<       res3(2) = glu3/z3*Pqg(x/z3) >*/
    d__1 = *x / z3;
    res3[1] = glu3 / z3 * pqg_(&d__1);
/*<       do 100 l=1,8 >*/
    for (l = 1; l <= 8; ++l) {
/*<       y=d*g(1,l) >*/
	y = d__ * g[l * 3 - 3];
/*<       z1 = e+y >*/
	z1 = e + y;
/*<       z2 = e-y >*/
	z2 = e - y;
/*<       call GLUON(z1,Q2,glu1) >*/
	gluon_(&z1, q2, &glu1);
/*<       res1(1) = glu1*fun(x0/z1,Q2) >*/
	d__1 = *x0 / z1;
	res1[0] = glu1 * fun_(&d__1, q2);
/*<       res1(2) = glu1/z1*Pqg(x/z1) >*/
	d__1 = *x / z1;
	res1[1] = glu1 / z1 * pqg_(&d__1);
/*<       call GLUON(z2,Q2,glu2) >*/
	gluon_(&z2, q2, &glu2);
/*<       res2(1) = glu2*fun(x0/z2,Q2) >*/
	d__1 = *x0 / z2;
	res2[0] = glu2 * fun_(&d__1, q2);
/*<       res2(2) = glu2/z2*Pqg(x/z2) >*/
	d__1 = *x / z2;
	res2[1] = glu2 / z2 * pqg_(&d__1);
/*<       do 20 I=1,2 >*/
	for (i__ = 1; i__ <= 2; ++i__) {
/*<          c(I) = res1(I) + res2(I) >*/
	    c__[i__ - 1] = res1[i__ - 1] + res2[i__ - 1];
/*<          gaus(I) = gaus(I)+c(I)*g(2,l) >*/
	    gaus[i__] += c__[i__ - 1] * g[l * 3 - 2];
/*<  20      gauskr(I) = gauskr(I)+c(I)*g(3,l) >*/
/* L20: */
	    gauskr[i__] += c__[i__ - 1] * g[l * 3 - 1];
	}
/*< 100   continue >*/
/* L100: */
    }
/*<       do 30 I=1,2 >*/
    for (i__ = 1; i__ <= 2; ++i__) {
/*<          gaus(I) = d*gaus(I) >*/
	gaus[i__] = d__ * gaus[i__];
/*<  30      gauskr(I) = d*(gauskr(I)+g39*res3(I)) >*/
/* L30: */
	gauskr[i__] = d__ * (gauskr[i__] + g39 * res3[i__ - 1]);
    }
/*<       return >*/
    return 0;
/*<       end >*/
} /* gauscxr_ */

/* *************************************************************************** */
/*<       DOUBLE PRECISION FUNCTION fun(X,Q2) >*/
doublereal fun_(doublereal *x, doublereal *q2)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static doublereal mb2, mc2, mbq, mcq, beta, beta2;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       integer flav >*/
/*<       common /flav/ flav >*/
/*<       common /mass/ mc,mb >*/
/*<       if (flav.EQ.1) then >*/
    if (flav_1.flav == 1) {
/*<         mc2 = mc*mc >*/
	mc2 = mass_1.mc * mass_1.mc;
/*<         beta2 = 1.d0-4.d0*mc2*X / ((1.d0-X)*Q2) >*/
	beta2 = 1. - mc2 * 4. * *x / ((1. - *x) * *q2);
/*<         if (beta2.GT.0.d0) then >*/
	if (beta2 > 0.) {
/*<            mcQ = 4.d0*mc2/Q2 >*/
	    mcq = mc2 * 4. / *q2;
/*<            beta = DSQRT(beta2) >*/
	    beta = sqrt(beta2);
/*<        >*/
	    ret_val = *x * (beta * (*x * 8. * (1. - *x) - 1. - *x * (1. - *x) 
		    * mcq) + (*x * *x + ((float)1. - *x) * ((float)1. - *x) + 
		    *x * ((float)1. - *x * (float)3.) * mcq - *x * *x * mcq * 
		    mcq / (float)2.) * log((beta + (float)1.) / ((float)1. - 
		    beta)));
/*<         else >*/
	} else {
/*<            fun = 0.d0 >*/
	    ret_val = 0.;
/*<         endif >*/
	}
/*<       elseif (flav.EQ.2) then >*/
    } else if (flav_1.flav == 2) {
/*<         mb2 = mb*mb >*/
	mb2 = mass_1.mb * mass_1.mb;
/*<         beta2 = 1.d0-4.d0*mb2*X / ((1.d0-X)*Q2) >*/
	beta2 = 1. - mb2 * 4. * *x / ((1. - *x) * *q2);
/*<         if (beta2.GT.0.d0) then >*/
	if (beta2 > 0.) {
/*<            mbQ = 4.d0*mb2/Q2 >*/
	    mbq = mb2 * 4. / *q2;
/*<            beta = DSQRT(beta2) >*/
	    beta = sqrt(beta2);
/*<        >*/
	    ret_val = *x * (beta * (*x * 8. * (1. - *x) - 1. - *x * (1. - *x) 
		    * mbq) + (*x * *x + ((float)1. - *x) * ((float)1. - *x) + 
		    *x * ((float)1. - *x * (float)3.) * mbq - *x * *x * mbq * 
		    mbq / (float)2.) * log((beta + (float)1.) / ((float)1. - 
		    beta)));
/*<         else >*/
	} else {
/*<            fun = 0.d0 >*/
	    ret_val = 0.;
/*<         endif >*/
	}
/*<       endif >*/
    }
/*<       END >*/
    return ret_val;
} /* fun_ */

/* *************************************************************************** */
/*<       DOUBLE PRECISION FUNCTION Pqg(X) >*/
doublereal pqg_(doublereal *x)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       Pqg = 1.d0/2.d0*( X*X + (1.d0-X)*(1.d0-X) ) >*/
    ret_val = (*x * *x + (1. - *x) * (1. - *x)) * .5;
/*<       END >*/
    return ret_val;
} /* pqg_ */

#ifdef __cplusplus
	}
#endif
