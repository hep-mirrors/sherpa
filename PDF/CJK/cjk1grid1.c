/* cjk1grid1.f -- translated by f2c (version 20200916).
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
    doublereal lam3;
} lam_;

#define lam_1 lam_

Extern struct {
    doublereal mc, mb;
} mass_;

#define mass_1 mass_

Extern struct {
    integer flav;
} flavv1_;

#define flavv1_1 flavv1_

Extern struct {
    integer ist;
} istv1_;

#define istv1_1 istv1_

Extern struct {
    integer iread;
} ireadv1_;

#define ireadv1_1 ireadv1_

Extern struct {
    doublereal gl[16524]	/* was [9][34][54] */, dn[16524]	/* 
	    was [9][34][54] */, up[16524]	/* was [9][34][54] */, st[
	    16524]	/* was [9][34][54] */, ch[24300]	/* was [9][50]
	    [54] */, bt[24300]	/* was [9][50][54] */, chx[10800]	/* 
	    was [9][50][24] */, btx[10800]	/* was [9][50][24] */;
} partv1_;

#define partv1_1 partv1_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;

/* ******************************************************************** */
/* ***     The same parametrization as in cjk1grid.f but for       **** */
/* ***    calculations according to the "ACOT2(chi)" approach      **** */
/* ***        described in "application.ps" available at           **** */
/* ***          http//www.fuw.edu.pl/~pjank/param.html             **** */
/* ***            (part of P.Jankowski's PhD thesis)               **** */
/* ***                                                             **** */
/* ***   LO parametrization of the parton densities in the real    **** */
/* ***             photon for the CJK 1 fit based on               **** */
/* ***                                                             **** */
/* ***            F.Cornet, P.Jankowski and M.Krawczyk             **** */
/* ***     "CJK - Improved LO Parton Distributions in the Real     **** */
/* ***      Photon and Their Experimental Uncertainties"           **** */
/* ***          Nucl. Phys. Proc. Suppl. 126: 28, 2004             **** */
/* ***                     hep-ph/0310029                          **** */
/* ***                                                             **** */
/* ***     with additional parametrizations of the test parton     **** */
/* ***  densities which allow for calculation of uncertainties of  **** */
/* ***            any observable X(parton densities).              **** */
/* ***     To read about the method used to obtain them see        **** */
/* ***                                                             **** */
/* ***                       P.Jankowski                           **** */
/* *** "Uncertainties of the CJK 5 Flavour LO Parton Distributions **** */
/* ***                   in the Real Photon"                       **** */
/* ***                     hep-ph/0312056                          **** */
/* ******************************************************************** */
/* ***                                                             **** */
/* ***   valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2      **** */
/* ***                                                             **** */
/* ***       x      - Bjorken x variable                           **** */
/* ***       xc     - chi_c calculated for the given process       **** */
/* ***       xb     - chi_b calculated for the given process       **** */
/* ***       Q2     - square of momentum scale (in GeV**2)         **** */
/* ***   XPDF(-5:5) - matrix containing x*f(x,Q2)/alfa             **** */
/* ***                                                             **** */
/* ***   PDF =   -5,  -4,    -3 ,  -2 ,  -1 ,0 ,1,2,3,4,5          **** */
/* ***         b_bar,c_bar,s_bar,u_bar,d_bar,gl,d,u,s,c,b          **** */
/* ***                                                             **** */
/* ***  All antiquark and corresponding quark densities are equal. **** */
/* ***                                                             **** */
/* ***   heavy-quark masses:                                       **** */
/* ***             M_charm=1.3, M_beauty=4.3 GeV                   **** */
/* ***                                                             **** */
/* ***   Lambda_QCD values for active N_f falvours:                **** */
/* ***   N_f =   3      4      5                                   **** */
/* ***         0.138  0.115  0.084 GeV                             **** */
/* ***                                                             **** */
/* ***  Grid parametrization utilizing the bicubic interpolation   **** */
/* ***  in the Hermite polynomials basis.                          **** */
/* ***                                                             **** */
/* ***  To use it one must add in ones main program:               **** */
/* ***        INTEGER IREAD                                        **** */
/* ***        common /IREADV1/ IREAD                               **** */
/* ***        IREAD = 0                                            **** */
/* ***  This allows for fast multiple use of the parametrization.  **** */
/* ***                                                             **** */
/* ***                                                             **** */
/* ***   IOPT = 1 : light partons (gl,up,dn,str)                   **** */
/* ***        = 2 : all partons   (gl,up,dn,str,chm,bot)           **** */
/* ***                                                             **** */
/* ***   ISET = 0 : bestfit S0 - "cjk1best.dat"                    **** */
/* ***                                                             **** */
/* ***   Other ISET values (other files of grid data) can be       **** */
/* ***   used to calculate uncertainties of any observable X(pdf)  **** */
/* ***   depending on photon parton densities.                     **** */
/* ***   One needs to use the master equation:                     **** */
/* ***                                                             **** */
/* ***      DX = T/10* (Sum_I [X(SI+)-X(SI-)]^2)^1/2    I=1..4     **** */
/* ***                    D\chi^2 = T^2                            **** */
/* ***   with                                                      **** */
/* ***        DX (T)      - uncertainty of X                       **** */
/* ***        D\chi^2 (T) - corresponding displacement from the    **** */
/* ***                      minimal (best) fit of parton densities **** */
/* ***                                                             **** */
/* ***   ISET = 1 : S1+  - "cjk1set1pl.dat"                        **** */
/* ***   ISET = 2 : S1-  - "cjk1set1mn.dat"                        **** */
/* ***   ISET = 3 : S2+  - "cjk1set2pl.dat"                        **** */
/* ***   ISET = 4 : S2-  - "cjk1set2mn.dat"                        **** */
/* ***   ISET = 5 : S3+  - "cjk1set3pl.dat"                        **** */
/* ***   ISET = 6 : S3-  - "cjk1set3mn.dat"                        **** */
/* ***   ISET = 7 : S4+  - "cjk1set4pl.dat"                        **** */
/* ***   ISET = 8 : S4-  - "cjk1set4mn.dat"                        **** */
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
/* ***  Last changes - 26 May 2004                                 **** */
/* ******************************************************************** */
/*<       SUBROUTINE CJK1GRID(ISET,IOPT,X,XC,XB,Q2,XPDF,F2alfa) >*/
/* Subroutine */ int cjk1grid_(integer *iset, integer *iopt, doublereal *x, 
	doublereal *xc, doublereal *xb, doublereal *q2, doublereal *xpdf, 
	doublereal *f2alfa)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, nf;
    static doublereal xx, mb2, mc2, qq2, mu2, xdn, xxb, xxc, xup, xchm, xbot, 
	    xglu, xstr;
    extern /* Subroutine */ int gridv1_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };


/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       INTEGER I,flav,Nf,step,IOPT,ISET,IST,IREAD,IROPT >*/
/*<        >*/
/*<       logical t >*/
/*<       dimension XPDF(-5:5),resres(2) >*/
/*<       common /Lam/ Lam3 >*/
/*<       common /mass/ mc,mb >*/
/*<       common /flavv1/ flav >*/
/*<       common /ISTV1/ IST >*/
/*<       IST = ISET >*/
    /* Parameter adjustments */
    xpdf -= -5;

    /* Function Body */
    istv1_1.ist = *iset;
/*<       XX = X >*/
    xx = *x;
/*<       QQ2 = Q2 >*/
    qq2 = *q2;
/*<       if ((XX.LE.1.d-5).OR.(XX.GT.1.d0)) then >*/
    if (xx <= 1e-5 || xx > 1.) {
/*<          print *,'X out of range: ',XX >*/
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, "X out of range: ", (ftnlen)16);
	do_lio(&c__5, &c__1, (char *)&xx, (ftnlen)sizeof(doublereal));
	e_wsle();
/*<          stop >*/
	s_stop("", (ftnlen)0);
/*<       endif >*/
    }
/*<       if (XX.EQ.1.d0) then >*/
    if (xx == 1.) {
/*<          do 10 I=1,5 >*/
	for (i__ = 1; i__ <= 5; ++i__) {
/*<  10         XPDF(I) = 0.d0 >*/
/* L10: */
	    xpdf[i__] = 0.;
	}
/*<          goto 1000 >*/
	goto L1000;
/*<       endif >*/
    }
/*<       if ((QQ2.LE.5.d-1).OR.(QQ2.GE.5.d5)) then >*/
    if (qq2 <= .5 || qq2 >= 5e5) {
/*<          print *,'Q2 out of range: ',QQ2 >*/
	s_wsle(&io___5);
	do_lio(&c__9, &c__1, "Q2 out of range: ", (ftnlen)17);
	do_lio(&c__5, &c__1, (char *)&qq2, (ftnlen)sizeof(doublereal));
	e_wsle();
/*<          stop >*/
	s_stop("", (ftnlen)0);
/*<       endif >*/
    }
/*<       mc = 1.3d0 >*/
    mass_1.mc = 1.3;
/*<       mc2 = mc*mc >*/
    mc2 = mass_1.mc * mass_1.mc;
/*<       mb = 4.3d0 >*/
    mass_1.mb = 4.3;
/*<       mb2 = mb*mb >*/
    mb2 = mass_1.mb * mass_1.mb;
/*<       XXC = XC/(1.d0+4.d0*mc2/Q2) >*/
    xxc = *xc / (mc2 * 4. / *q2 + 1.);
/*<       XXB = XB/(1.d0+4.d0*mb2/Q2) >*/
    xxb = *xb / (mb2 * 4. / *q2 + 1.);
/*<       MU2 = 0.25d0 >*/
    mu2 = .25;
/*<       if (QQ2.LT.mc2) then >*/
    if (qq2 < mc2) {
/*<          Nf = 3 >*/
	nf = 3;
/*<       elseif (QQ2.LT.mb2) then >*/
    } else if (qq2 < mb2) {
/*<          Nf = 4 >*/
	nf = 4;
/*<       else >*/
    } else {
/*<          Nf = 5 >*/
	nf = 5;
/*<       endif >*/
    }
/*<       if (IOPT.EQ.1) then >*/
    if (*iopt == 1) {
/*<          call GRIDV1(1,XX,XXC,XXB,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT) >*/
	gridv1_(&c__1, &xx, &xxc, &xxb, &qq2, &xglu, &xdn, &xup, &xstr, &xchm,
		 &xbot);
/*<          XPDF(0) = xglu/alfa >*/
	xpdf[0] = xglu / .00729735308;
/*<          XPDF(1) = xdn/alfa >*/
	xpdf[1] = xdn / .00729735308;
/*<          XPDF(2) = xup/alfa >*/
	xpdf[2] = xup / .00729735308;
/*<          XPDF(3) = xstr/alfa >*/
	xpdf[3] = xstr / .00729735308;
/*<          XPDF(4) = 0.d0 >*/
	xpdf[4] = 0.;
/*<          XPDF(5) = 0.d0 >*/
	xpdf[5] = 0.;
/*<          F2alfa  = 0.d0 >*/
	*f2alfa = 0.;
/*<          goto 1000 >*/
	goto L1000;
/*<       endif >*/
    }
/*<       if (IOPT.EQ.2) then >*/
    if (*iopt == 2) {
/*<          call GRIDV1(2,XX,XXC,XXB,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT) >*/
	gridv1_(&c__2, &xx, &xxc, &xxb, &qq2, &xglu, &xdn, &xup, &xstr, &xchm,
		 &xbot);
/*<          XPDF(0) = xglu/alfa >*/
	xpdf[0] = xglu / .00729735308;
/*<          XPDF(1) = xdn/alfa >*/
	xpdf[1] = xdn / .00729735308;
/*<          XPDF(2) = xup/alfa >*/
	xpdf[2] = xup / .00729735308;
/*<          XPDF(3) = xstr/alfa >*/
	xpdf[3] = xstr / .00729735308;
/*<          XPDF(4) = xchm/alfa >*/
	xpdf[4] = xchm / .00729735308;
/*<          XPDF(5) = xbot/alfa >*/
	xpdf[5] = xbot / .00729735308;
/*<          F2alfa  = 0.d0 >*/
	*f2alfa = 0.;
/*<       endif >*/
    }
/*<  1000 continue >*/
L1000:
/*<       do 20 I=1,5 >*/
    for (i__ = 1; i__ <= 5; ++i__) {
/*<  20      XPDF(-I) = XPDF(I) >*/
/* L20: */
	xpdf[-i__] = xpdf[i__];
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* cjk1grid_ */

/* ******************************************************************** */
/* ******************************************************************** */
/*<        >*/
/* Subroutine */ int gridv1_(integer *iopt, doublereal *xin, doublereal *xcin,
	 doublereal *xbin, doublereal *q2in, doublereal *xglu, doublereal *
	xdn, doublereal *xup, doublereal *xstr, doublereal *xchm, doublereal *
	xbot)
{
    /* Initialized data */

    static doublereal xdata[54] = { 0.,1e-5,2e-5,4e-5,6e-5,8e-5,1e-4,2e-4,
	    4e-4,6e-4,8e-4,.001,.002,.004,.006,.008,.01,.014,.02,.03,.04,.06,
	    .08,.1,.125,.15,.175,.2,.225,.25,.275,.3,.325,.35,.375,.4,.425,
	    .45,.475,.5,.525,.55,.575,.6,.65,.7,.75,.8,.85,.9,.95,.98,1.,0. };
    static doublereal q2data[34] = { 0.,.5,.75,1.,1.25,1.5,2.,2.5,3.2,4.,5.,
	    6.4,8.,10.,12.,18.,26.,40.,64.,100.,160.,240.,400.,640.,1e3,1800.,
	    3200.,5600.,1e4,1.8e4,3.2e4,5.6e4,1e5,0. };
    static doublereal q2hdata[50] = { 0.,.5,.75,1.,1.25,1.5,1.75,2.,2.25,2.5,
	    2.75,3.2,3.6,4.,4.5,5.,5.4,6.,6.4,7.2,8.,10.,12.,15.,18.,22.,26.,
	    33.,40.,52.,64.,72.,85.,100.,130.,160.,200.,240.,400.,640.,1e3,
	    1800.,3200.,5600.,1e4,1.8e4,3.2e4,5.6e4,1e5,0. };
    static integer imaxc[48] = { 21,22,23,24,25,26,27,28,29,29,30,31,32,33,34,
	    35,36,37,38,39,41,42,43,44,45,45,46,46,47,47,48,48,48,48,48,48,49,
	    49,49,49,49,49,49,49,49,49,49,49 };
    static integer imaxb[48] = { 14,15,16,17,17,18,18,18,19,19,19,20,20,20,20,
	    21,21,21,22,22,23,24,25,26,27,28,30,32,34,36,37,38,40,42,43,44,45,
	    46,47,48,48,49,49,49,49,49,49,49 };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    extern /* Subroutine */ int findxhv1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *), readtabv1_(integer *);
    static integer i__, j, k, l, m, n;
    static doublereal x[4], q2[4], mb, mc;
    static integer jh;
    static doublereal xb[4], xc[4], xx, mb2, mc2, q2h[4], qq2, chh[16]	/* 
	    was [4][4] */, dnh[16]	/* was [4][4] */, glh[16]	/* 
	    was [4][4] */, bth[16]	/* was [4][4] */, uph[16]	/* 
	    was [4][4] */, sth[16]	/* was [4][4] */, xxb, xxc;
    static integer nqq2, bord;
    static doublereal xmxb, xmxc;
    static integer nqq2h;
    extern /* Subroutine */ int fitv1_(integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static integer bordb, bordc;
    static doublereal xmaxb, xmaxc;
    extern /* Subroutine */ int findxv1_(doublereal *, doublereal *, integer *
	    ), findq2v1_(integer *, doublereal *, doublereal *, integer *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<        >*/
/*<       parameter (NX=52,NQ2=32,NQ2H=48,alfa=7.29735308D-3) >*/
/*<        >*/
/*<       common /IREADV1/ IREAD >*/
/*<       common /ISTV1/ IST >*/
/*<       common /PARTV1/ gl,dn,up,st,ch,bt,chx,btx >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<       XX = XIN >*/
    xx = *xin;
/*<       XXC = XCIN >*/
    xxc = *xcin;
/*<       XXB = XBIN >*/
    xxb = *xbin;
/*<       QQ2 = Q2IN >*/
    qq2 = *q2in;
/*<       mc = 1.3d0 >*/
    mc = 1.3;
/*<       mc2 = mc*mc >*/
    mc2 = mc * mc;
/*<       mb = 4.3d0 >*/
    mb = 4.3;
/*<       mb2 = mb*mb >*/
    mb2 = mb * mb;
/*<       xmaxc = 1.d0/(1.d0+4.d0*mc2/QQ2) >*/
    xmaxc = 1. / (mc2 * 4. / qq2 + 1.);
/*<       xmaxb = 1.d0/(1.d0+4.d0*mb2/QQ2) >*/
    xmaxb = 1. / (mb2 * 4. / qq2 + 1.);
/*<       xchm = 0.d0 >*/
    *xchm = 0.;
/*<       xbot = 0.d0 >*/
    *xbot = 0.;
/*<       bord = 1 >*/
    bord = 1;
/*<       bordc = 1 >*/
    bordc = 1;
/*<       bordb = 1 >*/
    bordb = 1;
/* **************************************************************** */
/* ***                  reading grid data                      **** */
/* **************************************************************** */
/*<       if (IREAD.EQ.0) then >*/
    if (ireadv1_1.iread == 0) {
/*<          call readtabv1(IOPT) >*/
	readtabv1_(iopt);
/*<       endif >*/
    }
/*<       IREAD = 100 >*/
    ireadv1_1.iread = 100;
/* **************************************************************** */
/* ***   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  **** */
/* ***    searching for the I such that: x(I) < XX < x(I+1)    **** */
/* ***            for the light quarks and gluon               **** */
/* **************************************************************** */
/*<       NQQ2 = NQ2 >*/
    nqq2 = 32;
/*<       call findq2v1(NQQ2,q2data,QQ2,J) >*/
    findq2v1_(&nqq2, q2data, &qq2, &j);
/*<       call findxv1(xdata,XX,I) >*/
    findxv1_(xdata, &xx, &i__);
/*<       if (I.EQ.1.OR.I.EQ.NX-1) bord = 0 >*/
    if (i__ == 1 || i__ == 51) {
	bord = 0;
    }
/*<       if (J.EQ.1.OR.J.EQ.NQ2-1) bord = 0 >*/
    if (j == 1 || j == 31) {
	bord = 0;
    }
/* **************************************************************** */
/* ***     *x1(I-1)   $x2(I)   $$xx   $x3(I+1)   *x4(I+2)      **** */
/* **************************************************************** */
/* ***              Only 3 points at borders!!!                **** */
/* **************************************************************** */
/*<       do 10 K=1,4 >*/
    for (k = 1; k <= 4; ++k) {
/*<          x(K)  = xdata(I+K-2) >*/
	x[k - 1] = xdata[i__ + k - 2];
/*<  10      q2(K) = q2data(J+K-2) >*/
/* L10: */
	q2[k - 1] = q2data[j + k - 2];
    }
/* **************************************************************** */
/* ***                  tbh(1,1) = tb(J-1,I-1)                 **** */
/* ***                  tbh(1,2) = tb(J-1,I)                   **** */
/* ***                  tbh(1,3) = tb(J-1,I+1)                 **** */
/* ***                  tbh(1,4) = tb(J-1,I+2) ...             **** */
/* **************************************************************** */
/*<       if (IOPT.EQ.4) then >*/
    if (*iopt == 4) {
/*<          do 20 K=J-1,J+2 >*/
	i__1 = j + 2;
	for (k = j - 1; k <= i__1; ++k) {
/*<             do 20 L=I-1,I+2 >*/
	    i__2 = i__ + 2;
	    for (l = i__ - 1; l <= i__2; ++l) {
/*<                M = K+2-J >*/
		m = k + 2 - j;
/*<                N = L+2-I >*/
		n = l + 2 - i__;
/*<  20            glh(M,N) = gl(IST,K,L) >*/
/* L20: */
		glh[m + (n << 2) - 5] = partv1_1.gl[istv1_1.ist + (k + l * 34)
			 * 9];
	    }
	}
/*<          call fitv1(bord,QQ2,XX,q2,x,glh,xglu) >*/
	fitv1_(&bord, &qq2, &xx, q2, x, glh, xglu);
/*<          xglu = alfa*xglu >*/
	*xglu *= .00729735308;
/*<       else >*/
    } else {
/*<          do 21 K=J-1,J+2 >*/
	i__2 = j + 2;
	for (k = j - 1; k <= i__2; ++k) {
/*<             do 21 L=I-1,I+2 >*/
	    i__1 = i__ + 2;
	    for (l = i__ - 1; l <= i__1; ++l) {
/*<                M = K+2-J >*/
		m = k + 2 - j;
/*<                N = L+2-I >*/
		n = l + 2 - i__;
/*<                glh(M,N) = gl(IST,K,L) >*/
		glh[m + (n << 2) - 5] = partv1_1.gl[istv1_1.ist + (k + l * 34)
			 * 9];
/*<                dnh(M,N) = dn(IST,K,L) >*/
		dnh[m + (n << 2) - 5] = partv1_1.dn[istv1_1.ist + (k + l * 34)
			 * 9];
/*<                uph(M,N) = up(IST,K,L) >*/
		uph[m + (n << 2) - 5] = partv1_1.up[istv1_1.ist + (k + l * 34)
			 * 9];
/*<  21            sth(M,N) = st(IST,K,L) >*/
/* L21: */
		sth[m + (n << 2) - 5] = partv1_1.st[istv1_1.ist + (k + l * 34)
			 * 9];
	    }
	}
/*<          call fitv1(bord,QQ2,XX,q2,x,glh,xglu) >*/
	fitv1_(&bord, &qq2, &xx, q2, x, glh, xglu);
/*<          call fitv1(bord,QQ2,XX,q2,x,dnh,xdn) >*/
	fitv1_(&bord, &qq2, &xx, q2, x, dnh, xdn);
/*<          call fitv1(bord,QQ2,XX,q2,x,uph,xup) >*/
	fitv1_(&bord, &qq2, &xx, q2, x, uph, xup);
/*<          call fitv1(bord,QQ2,XX,q2,x,sth,xstr) >*/
	fitv1_(&bord, &qq2, &xx, q2, x, sth, xstr);
/*<          xglu = alfa*xglu >*/
	*xglu *= .00729735308;
/*<          xdn  = alfa*xdn >*/
	*xdn *= .00729735308;
/*<          xup  = alfa*xup >*/
	*xup *= .00729735308;
/*<          xstr = alfa*xstr >*/
	*xstr *= .00729735308;
/*<       endif >*/
    }
/*<       if ((IOPT.EQ.1).OR.(IOPT.EQ.4)) goto 1000 >*/
    if (*iopt == 1 || *iopt == 4) {
	goto L1000;
    }
/* **************************************************************** */
/* ***   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  **** */
/* ***    searching for the I such that: x(I) < XX < x(I+1)    **** */
/* ***                  for the heavy quarks                   **** */
/* **************************************************************** */
/*<       NQQ2H = NQ2H >*/
    nqq2h = 48;
/*<       call findq2v1(NQQ2H,q2hdata,QQ2,JH) >*/
    findq2v1_(&nqq2h, q2hdata, &qq2, &jh);
/*<       if (JH.EQ.1.OR.JH.EQ.NQ2H-1) then >*/
    if (jh == 1 || jh == 47) {
/*<          bordc = 0 >*/
	bordc = 0;
/*<          bordb = 0 >*/
	bordb = 0;
/*<       endif >*/
    }
/*<       do 30 K=1,4 >*/
    for (k = 1; k <= 4; ++k) {
/*<  30      q2h(K) = q2hdata(JH+K-2) >*/
/* L30: */
	q2h[k - 1] = q2hdata[jh + k - 2];
    }
/*<       if (XXC.LT.xmaxc) then >*/
    if (xxc < xmaxc) {
/*<          xmxc = 1.d0/(1.d0+4.d0*mc2/q2h(2)) >*/
	xmxc = 1. / (mc2 * 4. / q2h[1] + 1.);
/*<          call findxhv1(xmxc,xdata,ch,chx,XXC,Imaxc(JH),JH,bordc,xc,chh) >*/
	findxhv1_(&xmxc, xdata, partv1_1.ch, partv1_1.chx, &xxc, &imaxc[jh - 
		1], &jh, &bordc, xc, chh);
/*<          call fitv1(bordc,QQ2,XXC,q2h,xc,chh,xchm) >*/
	fitv1_(&bordc, &qq2, &xxc, q2h, xc, chh, xchm);
/*<          xchm = XX/XXC*alfa*xchm >*/
	*xchm = xx / xxc * .00729735308 * *xchm;
/*<       endif >*/
    }
/*<       if (XXB.LT.xmaxb) then >*/
    if (xxb < xmaxb) {
/*<          xmxb = 1.d0/(1.d0+4.d0*mb2/q2h(2)) >*/
	xmxb = 1. / (mb2 * 4. / q2h[1] + 1.);
/*<          call findxhv1(xmxb,xdata,bt,btx,XXB,Imaxb(JH),JH,bordb,xb,bth) >*/
	findxhv1_(&xmxb, xdata, partv1_1.bt, partv1_1.btx, &xxb, &imaxb[jh - 
		1], &jh, &bordb, xb, bth);
/*<          call fitv1(bordb,QQ2,XXB,q2h,xb,bth,xbot) >*/
	fitv1_(&bordb, &qq2, &xxb, q2h, xb, bth, xbot);
/*<          xbot = XX/XXB*alfa*xbot >*/
	*xbot = xx / xxb * .00729735308 * *xbot;
/*<       endif >*/
    }
/*<  1000 continue >*/
L1000:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* gridv1_ */

/* **************************************************************** */
/*<       SUBROUTINE READTABV1(IOPT) >*/
/* Subroutine */ int readtabv1_(integer *iopt)
{
    /* Format strings */
    static char fmt_100[] = "(f11.7,3(2x,f10.7))";
    static char fmt_200[] = "(f10.7,2x,f10.7)";
    static char fmt_300[] = "(23(f10.7,2x),f10.7)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(), f_clos(cllist *);

    /* Local variables */
    static integer i__, j, ist;
    static char name__[15];

    /* Fortran I/O blocks */
    static cilist io___62 = { 0, 10, 0, fmt_100, 0 };
    static cilist io___63 = { 0, 10, 0, fmt_200, 0 };
    static cilist io___64 = { 0, 10, 0, fmt_300, 0 };
    static cilist io___65 = { 0, 10, 0, fmt_300, 0 };


/*<       DOUBLE PRECISION gl,dn,up,st,ch,bt,chx,btx >*/
/*<       parameter (NX=52,NQ2=32,NQ2H=48) >*/
/*<        >*/
/*<       character*15 name >*/
/*<       common /PARTV1/ gl,dn,up,st,ch,bt,chx,btx >*/
/*<       do 2000 IST = 0,8 >*/
    for (ist = 0; ist <= 8; ++ist) {
/*<          if (IST.EQ.0) name = 'cjk1best.dat' >*/
	if (ist == 0) {
	    s_copy(name__, "cjk1best.dat", (ftnlen)15, (ftnlen)12);
	}
/*<          if (IST.EQ.1) name = 'cjk1set1pl.dat' >*/
	if (ist == 1) {
	    s_copy(name__, "cjk1set1pl.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.2) name = 'cjk1set1mn.dat' >*/
	if (ist == 2) {
	    s_copy(name__, "cjk1set1mn.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.3) name = 'cjk1set2pl.dat' >*/
	if (ist == 3) {
	    s_copy(name__, "cjk1set2pl.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.4) name = 'cjk1set2mn.dat' >*/
	if (ist == 4) {
	    s_copy(name__, "cjk1set2mn.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.5) name = 'cjk1set3pl.dat' >*/
	if (ist == 5) {
	    s_copy(name__, "cjk1set3pl.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.6) name = 'cjk1set3mn.dat' >*/
	if (ist == 6) {
	    s_copy(name__, "cjk1set3mn.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.7) name = 'cjk1set4pl.dat' >*/
	if (ist == 7) {
	    s_copy(name__, "cjk1set4pl.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          if (IST.EQ.8) name = 'cjk1set4mn.dat' >*/
	if (ist == 8) {
	    s_copy(name__, "cjk1set4mn.dat", (ftnlen)15, (ftnlen)14);
	}
/*<          open(10,file=name,status='old') >*/
	o__1.oerr = 0;
	o__1.ounit = 10;
	o__1.ofnmlen = 15;
	o__1.ofnm = name__;
	o__1.orl = 0;
	o__1.osta = "old";
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
/*<          do 1 I=0,NQ2+1 >*/
	for (i__ = 0; i__ <= 33; ++i__) {
/*<             do 1 J=0,NX+1 >*/
	    for (j = 0; j <= 53; ++j) {
/*<                gl(IST,I,J) = 0.d0 >*/
		partv1_1.gl[ist + (i__ + j * 34) * 9] = 0.;
/*<                dn(IST,I,J) = 0.d0 >*/
		partv1_1.dn[ist + (i__ + j * 34) * 9] = 0.;
/*<                up(IST,I,J) = 0.d0 >*/
		partv1_1.up[ist + (i__ + j * 34) * 9] = 0.;
/*<  1             st(IST,I,J) = 0.d0 >*/
/* L1: */
		partv1_1.st[ist + (i__ + j * 34) * 9] = 0.;
	    }
	}
/*<          do 10 I=1,NQ2 >*/
	for (i__ = 1; i__ <= 32; ++i__) {
/*<             do 10 J=1,NX >*/
	    for (j = 1; j <= 52; ++j) {
/*<  10    >*/
/* L10: */
		s_rsfe(&io___62);
		do_fio(&c__1, (char *)&partv1_1.gl[ist + (i__ + j * 34) * 9], 
			(ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&partv1_1.dn[ist + (i__ + j * 34) * 9], 
			(ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&partv1_1.up[ist + (i__ + j * 34) * 9], 
			(ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&partv1_1.st[ist + (i__ + j * 34) * 9], 
			(ftnlen)sizeof(doublereal));
		e_rsfe();
	    }
	}
/*<  100          format (F11.7,3(2X,F10.7)) >*/
/*<          if ((IOPT.EQ.1).OR.(IOPT.EQ.4)) goto 1000 >*/
	if (*iopt == 1 || *iopt == 4) {
	    goto L1000;
	}
/*<          do 2 I=0,NQ2H+1 >*/
	for (i__ = 0; i__ <= 49; ++i__) {
/*<             do 2 J=0,NX+1 >*/
	    for (j = 0; j <= 53; ++j) {
/*<                ch(IST,I,J) = 0.d0 >*/
		partv1_1.ch[ist + (i__ + j * 50) * 9] = 0.;
/*<  2             bt(IST,I,J) = 0.d0 >*/
/* L2: */
		partv1_1.bt[ist + (i__ + j * 50) * 9] = 0.;
	    }
	}
/*<          do 20 I=1,NQ2H >*/
	for (i__ = 1; i__ <= 48; ++i__) {
/*<             do 20 J=1,NX >*/
	    for (j = 1; j <= 52; ++j) {
/*<  20            read(10,200) ch(IST,I,J),bt(IST,I,J) >*/
/* L20: */
		s_rsfe(&io___63);
		do_fio(&c__1, (char *)&partv1_1.ch[ist + (i__ + j * 50) * 9], 
			(ftnlen)sizeof(doublereal));
		do_fio(&c__1, (char *)&partv1_1.bt[ist + (i__ + j * 50) * 9], 
			(ftnlen)sizeof(doublereal));
		e_rsfe();
	    }
	}
/*<  200           format (F10.7,2X,F10.7) >*/
/*<          do 30 I=1,NQ2H >*/
	for (i__ = 1; i__ <= 48; ++i__) {
/*<             read(10,300) (chx(IST,I,J),J=1,24) >*/
	    s_rsfe(&io___64);
	    for (j = 1; j <= 24; ++j) {
		do_fio(&c__1, (char *)&partv1_1.chx[ist + (i__ + j * 50) * 9 
			- 450], (ftnlen)sizeof(doublereal));
	    }
	    e_rsfe();
/*<  30      continue >*/
/* L30: */
	}
/*<          do 40 I=1,NQ2H >*/
	for (i__ = 1; i__ <= 48; ++i__) {
/*<             read(10,300) (btx(IST,I,J),J=1,24) >*/
	    s_rsfe(&io___65);
	    for (j = 1; j <= 24; ++j) {
		do_fio(&c__1, (char *)&partv1_1.btx[ist + (i__ + j * 50) * 9 
			- 450], (ftnlen)sizeof(doublereal));
	    }
	    e_rsfe();
/*<  40      continue >*/
/* L40: */
	}
/*<  300     format (23(F10.7,2X),F10.7) >*/
/*<          do 31 I=0,NQ2H+1,NQ2H+1 >*/
	for (i__ = 0; i__ <= 49; i__ += 49) {
/*<             do 31 J=1,24 >*/
	    for (j = 1; j <= 24; ++j) {
/*<                chx(IST,I,J) = 0.d0 >*/
		partv1_1.chx[ist + (i__ + j * 50) * 9 - 450] = 0.;
/*<  31            btx(IST,I,J) = 0.d0 >*/
/* L31: */
		partv1_1.btx[ist + (i__ + j * 50) * 9 - 450] = 0.;
	    }
	}
/*<  1000    continue >*/
L1000:
/*<          close(10) >*/
	cl__1.cerr = 0;
	cl__1.cunit = 10;
	cl__1.csta = 0;
	f_clos(&cl__1);
/*<  2000 continue >*/
/* L2000: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* readtabv1_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDQ2V1(NQ2,q2data,QQ2,I) >*/
/* Subroutine */ int findq2v1_(integer *nq2, doublereal *q2data, doublereal *
	qq2, integer *i__)
{
    static doublereal q2;
    static integer il, iu;

/*<       DOUBLE PRECISION QQ2,Q2,q2data >*/
/*<       INTEGER I,iu,ul,NQ2 >*/
/*<       dimension q2data(0:NQ2+1) >*/
/*<       il = 1 >*/
    il = 1;
/*<       iu = NQ2 >*/
    iu = *nq2;
/*<  10   if (iu-il.GT.1) then >*/
L10:
    if (iu - il > 1) {
/*<          I = (iu+il)/2 >*/
	*i__ = (iu + il) / 2;
/*<          Q2 = q2data(I) >*/
	q2 = q2data[*i__];
/*<          if (QQ2.GE.Q2) then >*/
	if (*qq2 >= q2) {
/*<             il = I >*/
	    il = *i__;
/*<          else >*/
	} else {
/*<             iu = I >*/
	    iu = *i__;
/*<          endif >*/
	}
/*<       goto 10 >*/
	goto L10;
/*<       endif >*/
    }
/*<       I = il >*/
    *i__ = il;
/*<  100  continue >*/
/* L100: */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* findq2v1_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDXV1(xdata,XX,I) >*/
/* Subroutine */ int findxv1_(doublereal *xdata, doublereal *xx, integer *i__)
{
    static doublereal x;
    static integer il, iu;

/*<       DOUBLE PRECISION XX,x,xdata >*/
/*<       INTEGER I,il,iu >*/
/*<       parameter (NX=52) >*/
/*<       dimension xdata(0:NX+1) >*/
/*<       il = 1 >*/
    il = 1;
/*<       iu = NX >*/
    iu = 52;
/*<  10   if (iu-il.GT.1) then >*/
L10:
    if (iu - il > 1) {
/*<          I = (iu+il)/2 >*/
	*i__ = (iu + il) / 2;
/*<          x = xdata(I) >*/
	x = xdata[*i__];
/*<          if (XX.GE.x) then >*/
	if (*xx >= x) {
/*<             il = I >*/
	    il = *i__;
/*<          else >*/
	} else {
/*<             iu = I >*/
	    iu = *i__;
/*<          endif >*/
	}
/*<       goto 10 >*/
	goto L10;
/*<       endif >*/
    }
/*<       I = il >*/
    *i__ = il;
/*<  100  continue >*/
/* L100: */
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* findxv1_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDXHV1(xmax,xdata,hq,hqx,XX,Imax,J,bord,x,hqh) >*/
/* Subroutine */ int findxhv1_(doublereal *xmax, doublereal *xdata, 
	doublereal *hq, doublereal *hqx, doublereal *xx, integer *imax, 
	integer *j, integer *bord, doublereal *x, doublereal *hqh)
{
    /* Initialized data */

    static doublereal per[7] = { .93,.95,.96,.97,.99,.999,0. };

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    static integer i__, k, l, m, n, il;
    static doublereal xb;
    static integer iu;
    static doublereal xdatah[108];

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       INTEGER I,J,K,L,M,N,il,iu,Imax,bord,IST,NX,NQ2 >*/
/*<       parameter (NX=52,NQ2=48) >*/
/*<        >*/
/*<       common /ISTV1/ IST >*/
/*<       data per/0.93d0,0.95d0,0.96d0,0.97d0,0.99d0,0.999d0,0d0/ >*/
    /* Parameter adjustments */
    hqh -= 5;
    --x;
    hqx -= 450;

    /* Function Body */
/*<       do 1 I=0,Imax >*/
    i__1 = *imax;
    for (i__ = 0; i__ <= i__1; ++i__) {
/*<  1       xdatah(I) = xdata(I) >*/
/* L1: */
	xdatah[i__] = xdata[i__];
    }
/*<       do 2 I=1,7 >*/
    for (i__ = 1; i__ <= 7; ++i__) {
/*<  2       xdatah(Imax+I) = per(I)*xmax >*/
/* L2: */
	xdatah[*imax + i__] = per[i__ - 1] * *xmax;
    }
/*<       il = 1 >*/
    il = 1;
/*<       iu = Imax+6 >*/
    iu = *imax + 6;
/*<  10   if (iu-il.GT.1) then >*/
L10:
    if (iu - il > 1) {
/*<          I = (iu+il)/2 >*/
	i__ = (iu + il) / 2;
/*<          xb = xdatah(I) >*/
	xb = xdatah[i__];
/*<          if (XX.GE.xb) then >*/
	if (*xx >= xb) {
/*<             il = I >*/
	    il = i__;
/*<          else >*/
	} else {
/*<             iu = I >*/
	    iu = i__;
/*<          endif >*/
	}
/*<          goto 10 >*/
	goto L10;
/*<        endif >*/
    }
/*<        I = il >*/
    i__ = il;
/*<        K = I-Imax+1 >*/
    k = i__ - *imax + 1;
/*<        if (K.LT.0) K = 0 >*/
    if (k < 0) {
	k = 0;
    }
/* **************************************************************** */
/* ***     *x1(I-1)   $x2(I)   $$xx   $x3(I+1)   *x4(I+2)      **** */
/* **************************************************************** */
/* ***              Only 3 points at borders!!!                **** */
/* **************************************************************** */
/*<       do 20 L=1,4 >*/
    for (l = 1; l <= 4; ++l) {
/*<  20      x(L) = xdatah(I+L-2) >*/
/* L20: */
	x[l] = xdatah[i__ + l - 2];
    }
/*<       if (K.EQ.0) then >*/
    if (k == 0) {
/*<          if (I+2.LE.Imax) then >*/
	if (i__ + 2 <= *imax) {
/*<             do 30 M=J-1,J+2 >*/
	    i__1 = *j + 2;
	    for (m = *j - 1; m <= i__1; ++m) {
/*<                do 30 L=I-1,I+2 >*/
		i__2 = i__ + 2;
		for (l = i__ - 1; l <= i__2; ++l) {
/*<  30               hqh(M+2-J,L+2-I) = hq(IST,M,L) >*/
/* L30: */
		    hqh[m + 2 - *j + (l + 2 - i__ << 2)] = hq[istv1_1.ist + (
			    m + l * 50) * 9];
		}
	    }
/*<          else >*/
	} else {
/*<             do 40 M=J-1,J+2 >*/
	    i__2 = *j + 2;
	    for (m = *j - 1; m <= i__2; ++m) {
/*<                N = M+2-J >*/
		n = m + 2 - *j;
/*<                hqh(N,4) = hqx(IST,M,6*N-5) >*/
		hqh[n + 16] = hqx[istv1_1.ist + (m + (n * 6 - 5) * 50) * 9];
/*<                do 40 L=I-1,I+1 >*/
		i__1 = i__ + 1;
		for (l = i__ - 1; l <= i__1; ++l) {
/*<  40               hqh(N,L+2-I) = hq(IST,M,L) >*/
/* L40: */
		    hqh[n + (l + 2 - i__ << 2)] = hq[istv1_1.ist + (m + l * 
			    50) * 9];
		}
	    }
/*<          endif >*/
	}
/*<       elseif (K.EQ.1) then >*/
    } else if (k == 1) {
/*<          do 50 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             hqh(N,1) = hq(IST,M,Imax-1) >*/
	    hqh[n + 4] = hq[istv1_1.ist + (m + (*imax - 1) * 50) * 9];
/*<             hqh(N,2) = hq(IST,M,Imax) >*/
	    hqh[n + 8] = hq[istv1_1.ist + (m + *imax * 50) * 9];
/*<             hqh(N,3) = hqx(IST,M,6*N-5) >*/
	    hqh[n + 12] = hqx[istv1_1.ist + (m + (n * 6 - 5) * 50) * 9];
/*<  50         hqh(N,4) = hqx(IST,M,6*N-4) >*/
/* L50: */
	    hqh[n + 16] = hqx[istv1_1.ist + (m + (n * 6 - 4) * 50) * 9];
	}
/*<       elseif (K.EQ.2) then >*/
    } else if (k == 2) {
/*<          do 60 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             hqh(N,1) = hq(IST,M,Imax) >*/
	    hqh[n + 4] = hq[istv1_1.ist + (m + *imax * 50) * 9];
/*<             do 60 I=2,4 >*/
	    for (i__ = 2; i__ <= 4; ++i__) {
/*<  60            hqh(N,I) = hqx(IST,M,6*N+I-7) >*/
/* L60: */
		hqh[n + (i__ << 2)] = hqx[istv1_1.ist + (m + (n * 6 + i__ - 7)
			 * 50) * 9];
	    }
	}
/*<          if (hqh(1,3).LT.1.d-10) bord = 0 >*/
	if (hqh[13] < 1e-10) {
	    *bord = 0;
	}
/*<       elseif (K.EQ.3) then >*/
    } else if (k == 3) {
/*<          do 70 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             do 70 I=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<  70            hqh(N,I) = hqx(IST,M,6*N+I-6) >*/
/* L70: */
		hqh[n + (i__ << 2)] = hqx[istv1_1.ist + (m + (n * 6 + i__ - 6)
			 * 50) * 9];
	    }
	}
/*<          if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0 >*/
	if (hqh[13] < 1e-10 || hqh[18] < 1e-10) {
	    *bord = 0;
	}
/*<       elseif (K.EQ.4) then >*/
    } else if (k == 4) {
/*<          do 80 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             do 80 I=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<  80            hqh(N,I) = hqx(IST,M,6*N+I-5) >*/
/* L80: */
		hqh[n + (i__ << 2)] = hqx[istv1_1.ist + (m + (n * 6 + i__ - 5)
			 * 50) * 9];
	    }
	}
/*<          if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0 >*/
	if (hqh[13] < 1e-10 || hqh[18] < 1e-10) {
	    *bord = 0;
	}
/*<       elseif (K.EQ.5) then >*/
    } else if (k == 5) {
/*<          do 90 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             do 90 I=1,4 >*/
	    for (i__ = 1; i__ <= 4; ++i__) {
/*<  90            hqh(N,I) = hqx(IST,M,6*N+I-4) >*/
/* L90: */
		hqh[n + (i__ << 2)] = hqx[istv1_1.ist + (m + (n * 6 + i__ - 4)
			 * 50) * 9];
	    }
	}
/*<          if (hqh(1,3).LT.1.d-10.OR.hqh(2,4).LT.1.d-10) bord = 0 >*/
	if (hqh[13] < 1e-10 || hqh[18] < 1e-10) {
	    *bord = 0;
	}
/*<       elseif (K.EQ.6) then >*/
    } else if (k == 6) {
/*<          bord = 0 >*/
	*bord = 0;
/*<          do 100 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             hqh(N,4) = 0.d0 >*/
	    hqh[n + 16] = 0.;
/*<             do 100 I=1,3 >*/
	    for (i__ = 1; i__ <= 3; ++i__) {
/*<  100           hqh(N,I) = hqx(IST,M,6*N+I-3) >*/
/* L100: */
		hqh[n + (i__ << 2)] = hqx[istv1_1.ist + (m + (n * 6 + i__ - 3)
			 * 50) * 9];
	    }
	}
/*<       endif >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* findxhv1_ */

/* ******************************************************************** */
/* ***         Here I use the bicubic interpolation in             **** */
/* ***              the Hermit polynomials basis                   **** */
/* ******************************************************************** */
/*<       SUBROUTINE FITV1(bord,xx1,xx2,x1,x2,yg,result) >*/
/* Subroutine */ int fitv1_(integer *bord, doublereal *xx1, doublereal *xx2, 
	doublereal *x1, doublereal *x2, doublereal *yg, doublereal *result)
{
    static doublereal y[4], y1[4], y2[4], y12[4], x1l, x2l, x1u, x2u;
    extern doublereal d1fv1_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), d2fv1_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int iterv1_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension x1(4),x2(4),yg(4,4),y(4),y1(4),y2(4),y12(4) >*/
/*<       integer bord >*/
/*<       external d1fv1,d2fv1 >*/
/* **************************************************************** */
/* ***              4 *(x1l,x2u)    3 *(x1u,x2u)               **** */
/* ***              1 *(x1l,x2l)    2 *(x1u,x2l)               **** */
/* **************************************************************** */
/*<       x1l = x1(2) >*/
    /* Parameter adjustments */
    yg -= 5;
    --x2;
    --x1;

    /* Function Body */
    x1l = x1[2];
/*<       x1u = x1(3) >*/
    x1u = x1[3];
/*<       x2l = x2(2) >*/
    x2l = x2[2];
/*<       x2u = x2(3) >*/
    x2u = x2[3];
/* **************************************************************** */
/* ***      Function values, first and cross-derivatives       **** */
/* ***              at 4 corners of a grid cell                **** */
/* **************************************************************** */
/*<       y(1) = yg(2,2) >*/
    y[0] = yg[10];
/*<       y(2) = yg(3,2) >*/
    y[1] = yg[11];
/*<       y(3) = yg(3,3) >*/
    y[2] = yg[15];
/*<       y(4) = yg(2,3) >*/
    y[3] = yg[14];
/*<       if (bord.EQ.1) then >*/
    if (*bord == 1) {
/*<          y1(1) = d1fv1(x1(1),x1(3),yg(1,2),yg(3,2)) >*/
	y1[0] = d1fv1_(&x1[1], &x1[3], &yg[9], &yg[11]);
/*<          y1(2) = d1fv1(x1(2),x1(4),yg(2,2),yg(4,2)) >*/
	y1[1] = d1fv1_(&x1[2], &x1[4], &yg[10], &yg[12]);
/*<          y1(3) = d1fv1(x1(2),x1(4),yg(2,3),yg(4,3)) >*/
	y1[2] = d1fv1_(&x1[2], &x1[4], &yg[14], &yg[16]);
/*<          y1(4) = d1fv1(x1(1),x1(3),yg(1,3),yg(3,3)) >*/
	y1[3] = d1fv1_(&x1[1], &x1[3], &yg[13], &yg[15]);
/*<          y2(1) = d1fv1(x2(1),x2(3),yg(2,1),yg(2,3)) >*/
	y2[0] = d1fv1_(&x2[1], &x2[3], &yg[6], &yg[14]);
/*<          y2(2) = d1fv1(x2(1),x2(3),yg(3,1),yg(3,3)) >*/
	y2[1] = d1fv1_(&x2[1], &x2[3], &yg[7], &yg[15]);
/*<          y2(3) = d1fv1(x2(2),x2(4),yg(3,2),yg(3,4)) >*/
	y2[2] = d1fv1_(&x2[2], &x2[4], &yg[11], &yg[19]);
/*<          y2(4) = d1fv1(x2(2),x2(4),yg(2,2),yg(2,4)) >*/
	y2[3] = d1fv1_(&x2[2], &x2[4], &yg[10], &yg[18]);
/*<        >*/
	y12[0] = d2fv1_(&x1[1], &x1[3], &x2[1], &x2[3], &yg[5], &yg[13], &yg[
		7], &yg[15]);
/*<        >*/
	y12[1] = d2fv1_(&x1[2], &x1[4], &x2[1], &x2[3], &yg[6], &yg[14], &yg[
		8], &yg[16]);
/*<        >*/
	y12[2] = d2fv1_(&x1[2], &x1[4], &x2[2], &x2[4], &yg[10], &yg[18], &yg[
		12], &yg[20]);
/*<        >*/
	y12[3] = d2fv1_(&x1[1], &x1[3], &x2[2], &x2[4], &yg[9], &yg[17], &yg[
		11], &yg[19]);
/*<       else >*/
    } else {
/*<          y1(1) = d1fv1(x1(2),x1(3),yg(2,2),yg(3,2)) >*/
	y1[0] = d1fv1_(&x1[2], &x1[3], &yg[10], &yg[11]);
/*<          y1(2) = y1(1) >*/
	y1[1] = y1[0];
/*<          y1(3) = d1fv1(x1(2),x1(3),yg(2,3),yg(3,3)) >*/
	y1[2] = d1fv1_(&x1[2], &x1[3], &yg[14], &yg[15]);
/*<          y1(4) = y1(3) >*/
	y1[3] = y1[2];
/*<          y2(1) = d1fv1(x2(2),x2(3),yg(2,2),yg(2,3)) >*/
	y2[0] = d1fv1_(&x2[2], &x2[3], &yg[10], &yg[14]);
/*<          y2(2) = d1fv1(x2(2),x2(3),yg(3,2),yg(3,3)) >*/
	y2[1] = d1fv1_(&x2[2], &x2[3], &yg[11], &yg[15]);
/*<          y2(3) = y2(2) >*/
	y2[2] = y2[1];
/*<          y2(4) = y2(1) >*/
	y2[3] = y2[0];
/*<        >*/
	y12[0] = d2fv1_(&x1[2], &x1[3], &x2[2], &x2[3], &yg[10], &yg[14], &yg[
		11], &yg[15]);
/*<          y12(2) = y12(1) >*/
	y12[1] = y12[0];
/*<          y12(3) = y12(1) >*/
	y12[2] = y12[0];
/*<          y12(4) = y12(1) >*/
	y12[3] = y12[0];
/*<       endif >*/
    }
/*<       call iterv1(y,y1,y2,y12,x1l,x1u,x2l,x2u,xx1,xx2,result) >*/
    iterv1_(y, y1, y2, y12, &x1l, &x1u, &x2l, &x2u, xx1, xx2, result);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* fitv1_ */

/* ******************************************************************** */
/*<       SUBROUTINE iterv1(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,res) >*/
/* Subroutine */ int iterv1_(doublereal *y, doublereal *y1, doublereal *y2, 
	doublereal *y12, doublereal *x1l, doublereal *x1u, doublereal *x2l, 
	doublereal *x2u, doublereal *x1, doublereal *x2, doublereal *res)
{
    static integer i__, j;
    static doublereal p[16]	/* was [4][4] */, t, u, d1, d2, t1, u1, ph[4],
	     ht[4], hu[4], d1d2;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       DIMENSION y(4),y1(4),y2(4),y12(4),p(4,4),c(4,4),ht(4),hu(4),ph(4) >*/
/*<       INTEGER I,J >*/
/*<       d1 = x1u-x1l >*/
    /* Parameter adjustments */
    --y12;
    --y2;
    --y1;
    --y;

    /* Function Body */
    d1 = *x1u - *x1l;
/*<       d2 = x2u-x2l >*/
    d2 = *x2u - *x2l;
/*<       d1d2 = d1*d2 >*/
    d1d2 = d1 * d2;
/* ***  local variables  **** */
/*<       t = (x1-x1l)/d1 >*/
    t = (*x1 - *x1l) / d1;
/*<       u = (x2-x2l)/d2 >*/
    u = (*x2 - *x2l) / d2;
/* ***  local derivatives  **** */
/*<       p(1,1) = y(1) >*/
    p[0] = y[1];
/*<       p(2,1) = y(2) >*/
    p[1] = y[2];
/*<       p(2,2) = y(3) >*/
    p[5] = y[3];
/*<       p(1,2) = y(4) >*/
    p[4] = y[4];
/*<       p(3,1) = y1(1)*d1 >*/
    p[2] = y1[1] * d1;
/*<       p(4,1) = y1(2)*d1 >*/
    p[3] = y1[2] * d1;
/*<       p(4,2) = y1(3)*d1 >*/
    p[7] = y1[3] * d1;
/*<       p(3,2) = y1(4)*d1 >*/
    p[6] = y1[4] * d1;
/*<       p(1,3) = y2(1)*d2 >*/
    p[8] = y2[1] * d2;
/*<       p(2,3) = y2(2)*d2 >*/
    p[9] = y2[2] * d2;
/*<       p(2,4) = y2(3)*d2 >*/
    p[13] = y2[3] * d2;
/*<       p(1,4) = y2(4)*d2 >*/
    p[12] = y2[4] * d2;
/*<       p(3,3) = y12(1)*d1d2 >*/
    p[10] = y12[1] * d1d2;
/*<       p(4,3) = y12(2)*d1d2 >*/
    p[11] = y12[2] * d1d2;
/*<       p(4,4) = y12(3)*d1d2 >*/
    p[15] = y12[3] * d1d2;
/*<       p(3,4) = y12(4)*d1d2 >*/
    p[14] = y12[4] * d1d2;
/* ***  Hermite polynomials  **** */
/*<       t1 = t-1. >*/
    t1 = t - (float)1.;
/*<       ht(1) = (2.*t+1.)*t1*t1 >*/
    ht[0] = (t * (float)2. + (float)1.) * t1 * t1;
/*<       ht(2) = t*t*(-2.*t+3.) >*/
    ht[1] = t * t * (t * (float)-2. + (float)3.);
/*<       ht(3) = t*t1*t1 >*/
    ht[2] = t * t1 * t1;
/*<       ht(4) = t*t*t1 >*/
    ht[3] = t * t * t1;
/*<       u1 = u-1 >*/
    u1 = u - 1;
/*<       hu(1) = (2.*u+1.)*u1*u1 >*/
    hu[0] = (u * (float)2. + (float)1.) * u1 * u1;
/*<       hu(2) = u*u*(-2.*u+3.) >*/
    hu[1] = u * u * (u * (float)-2. + (float)3.);
/*<       hu(3) = u*u1*u1 >*/
    hu[2] = u * u1 * u1;
/*<       hu(4) = u*u*u1 >*/
    hu[3] = u * u * u1;
/*<       do 10 I=1,4 >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<          ph(I) = 0.d0 >*/
	ph[i__ - 1] = 0.;
/*<          do 10 J=1,4 >*/
	for (j = 1; j <= 4; ++j) {
/*<  10         ph(I) = ph(I) + p(I,J)*hu(J) >*/
/* L10: */
	    ph[i__ - 1] += p[i__ + (j << 2) - 5] * hu[j - 1];
	}
    }
/*<       res = 0.d0 >*/
    *res = 0.;
/*<       do 20 I=1,4 >*/
    for (i__ = 1; i__ <= 4; ++i__) {
/*<  20      res = res + ht(I)*ph(I) >*/
/* L20: */
	*res += ht[i__ - 1] * ph[i__ - 1];
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* iterv1_ */

/* ******************************************************************** */
/* ***                   First derivative: df/dx                   **** */
/* ******************************************************************** */
/*<       DOUBLE PRECISION FUNCTION d1fv1(x1,x2,f1,f2) >*/
doublereal d1fv1_(doublereal *x1, doublereal *x2, doublereal *f1, doublereal *
	f2)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       d1fv1 = (f2-f1)/(x2-x1) >*/
    ret_val = (*f2 - *f1) / (*x2 - *x1);
/*<       END >*/
    return ret_val;
} /* d1fv1_ */

/* ******************************************************************** */
/* ***                Second derivative: d2f/dxdy                  **** */
/* ******************************************************************** */
/*<       DOUBLE PRECISION FUNCTION d2fv1(x1,x2,y1,y2,f11,f12,f21,f22) >*/
doublereal d2fv1_(doublereal *x1, doublereal *x2, doublereal *y1, doublereal *
	y2, doublereal *f11, doublereal *f12, doublereal *f21, doublereal *
	f22)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       d2fv1 = (f22 - f21 - f12 + f11)/((x2-x1)*(y2-y1)) >*/
    ret_val = (*f22 - *f21 - *f12 + *f11) / ((*x2 - *x1) * (*y2 - *y1));
/*<       END >*/
    return ret_val;
} /* d2fv1_ */

#ifdef __cplusplus
	}
#endif
