/* cjkhogrid.f -- translated by f2c (version 20200916).
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
    integer iread;
} ireadvh2_;

#define ireadvh2_1 ireadvh2_

Extern struct {
    doublereal gl[1836]	/* was [34][54] */, dn[1836]	/* was [34][54] */, 
	    up[1836]	/* was [34][54] */, st[1836]	/* was [34][54] */, 
	    ch[2700]	/* was [50][54] */, bt[2700]	/* was [50][54] */, 
	    chx[1200]	/* was [50][24] */, btx[1200]	/* was [50][24] */;
} partvh2_;

#define partvh2_1 partvh2_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__2 = 2;

/* ******************************************************************** */
/* ***   NLO parametrization of the parton densities in the real   **** */
/* ***            photon for the CJK NLO fit based on              **** */
/* ***    "A New 5 Flavour NLO Analysis and Parametrizations       **** */
/* ***       of Parton Distributions of the Real Photon"           **** */
/* ***           by F.Cornet, P.Jankowski & M.Krawczyk             **** */
/* ***                     hep-ph/0404063                          **** */
/* ******************************************************************** */
/* ***                                                             **** */
/* ***   valid for 10^(-5) < x < 1 and 1 < Q^2 < 2*10^5 GeV^2      **** */
/* ***                                                             **** */
/* ***       x      - Bjorken x variable                           **** */
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
/* ***         0.323  0.280  0.200 GeV                             **** */
/* ***                                                             **** */
/* ***  Grid parametrization utilizing the bicubic interpolation   **** */
/* ***  in the Hermite polynomials basis.                          **** */
/* ***                                                             **** */
/* ***  To use it one must add in ones main program:               **** */
/* ***        INTEGER IREAD                                        **** */
/* ***        common /IREADVH2/ IREAD                              **** */
/* ***        IREAD = 0                                            **** */
/* ***  This allows for fast multiple use of the parametrization.  **** */
/* ***                                                             **** */
/* ***                                                             **** */
/* ***   IOPT = 1 : light partons (gl,up,dn,str)                   **** */
/* ***        = 2 : all partons   (gl,up,dn,str,chm,bot)           **** */
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
/* ***  Last changes - 07 April 2004                               **** */
/* ******************************************************************** */
/*<       SUBROUTINE CJKHOGRID(IOPT,X,Q2,XPDF) >*/
/* Subroutine */ int cjkhogrid_(integer *iopt, doublereal *x, doublereal *q2, 
	doublereal *xpdf)
{
    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static doublereal i__, xx, qq2, xdn, xup, xchm, xbot, xglu, xstr;
    extern /* Subroutine */ int gridvh2_(integer *, doublereal *, doublereal *
	    , doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };


/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       INTEGER IOPT >*/
/*<        >*/
/*<       logical t >*/
/*<       dimension XPDF(-5:5) >*/
/*<       XX = X >*/
    /* Parameter adjustments */
    xpdf -= -5;

    /* Function Body */
    xx = *x;
/*<       QQ2 = Q2 >*/
    qq2 = *q2;
/*<       if ((XX.LE.1.d-5).OR.(XX.GE.1.d0)) then >*/
    if (xx <= 1e-5 || xx >= 1.) {
/*<          print *,'X out of range: ',XX >*/
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, "X out of range: ", (ftnlen)16);
	do_lio(&c__5, &c__1, (char *)&xx, (ftnlen)sizeof(doublereal));
	e_wsle();
/*<          stop >*/
	s_stop("", (ftnlen)0);
/*<       endif >*/
    }
/*<       if ((QQ2.LE.5.d-1).OR.(QQ2.GE.5.d5)) then >*/
    if (qq2 <= .5 || qq2 >= 5e5) {
/*<          print *,'Q2 out of range: ',QQ2 >*/
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, "Q2 out of range: ", (ftnlen)17);
	do_lio(&c__5, &c__1, (char *)&qq2, (ftnlen)sizeof(doublereal));
	e_wsle();
/*<          stop >*/
	s_stop("", (ftnlen)0);
/*<       endif >*/
    }
/*<       if (IOPT.EQ.1) then >*/
    if (*iopt == 1) {
/*<          call GRIDVH2(1,XX,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT) >*/
	gridvh2_(&c__1, &xx, &qq2, &xglu, &xdn, &xup, &xstr, &xchm, &xbot);
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
/*<          goto 1000 >*/
	goto L1000;
/*<       endif >*/
    }
/*<       if (IOPT.EQ.2) then >*/
    if (*iopt == 2) {
/*<          call GRIDVH2(2,XX,QQ2,XGLU,XDN,XUP,XSTR,XCHM,XBOT) >*/
	gridvh2_(&c__2, &xx, &qq2, &xglu, &xdn, &xup, &xstr, &xchm, &xbot);
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
/*<       endif >*/
    }
/*<  1000 continue >*/
L1000:
/*<       do 10 I=1,5 >*/
    for (i__ = 1.; i__ <= 5.; i__ += 1.) {
/*<  10      XPDF(-I) = XPDF(I) >*/
/* L10: */
	xpdf[(integer) (-i__)] = xpdf[(integer) i__];
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* cjkhogrid_ */

/* ******************************************************************** */
/* ******************************************************************** */
/*<       SUBROUTINE GRIDVH2(IOPT,XIN,Q2IN,XGLU,XDN,XUP,XSTR,XCHM,XBOT) >*/
/* Subroutine */ int gridvh2_(integer *iopt, doublereal *xin, doublereal *
	q2in, doublereal *xglu, doublereal *xdn, doublereal *xup, doublereal *
	xstr, doublereal *xchm, doublereal *xbot)
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
    extern /* Subroutine */ int findxvh2_(doublereal *, doublereal *, integer 
	    *), findq2vh2_(integer *, doublereal *, doublereal *, integer *);
    static integer i__, j, k, l, m, n;
    extern /* Subroutine */ int findxhvh2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, doublereal *, doublereal *);
    static doublereal x[4], q2[4];
    extern /* Subroutine */ int readtabvh2_(integer *);
    static doublereal mb, mc;
    static integer jh;
    static doublereal xb[4], xc[4], xx, mb2, mc2, q2h[4], qq2, chh[16]	/* 
	    was [4][4] */, dnh[16]	/* was [4][4] */, glh[16]	/* 
	    was [4][4] */, bth[16]	/* was [4][4] */, uph[16]	/* 
	    was [4][4] */, sth[16]	/* was [4][4] */;
    static integer nqq2, bord;
    static doublereal xmxb, xmxc;
    static integer nqq2h, bordb, bordc;
    static doublereal xmaxb, xmaxc;
    extern /* Subroutine */ int fitvh2_(integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, doublereal *, doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<        >*/
/*<       parameter (NX=52,NQ2=32,NQ2H=48,alfa=7.29735308D-3) >*/
/*<        >*/
/*<       common /IREADVH2/ IREAD >*/
/*<       common /PARTVH2/ gl,dn,up,st,ch,bt,chx,btx >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<        >*/
/*<       XX = XIN >*/
    xx = *xin;
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
    if (ireadvh2_1.iread == 0) {
/*<          call readtabvh2(IOPT) >*/
	readtabvh2_(iopt);
/*<       endif >*/
    }
/*<       IREAD = 100 >*/
    ireadvh2_1.iread = 100;
/* **************************************************************** */
/* ***   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  **** */
/* ***    searching for the I such that: x(I) < XX < x(I+1)    **** */
/* ***            for the light quarks and gluon               **** */
/* **************************************************************** */
/*<       NQQ2 = NQ2 >*/
    nqq2 = 32;
/*<       call findq2vh2(NQQ2,q2data,QQ2,J) >*/
    findq2vh2_(&nqq2, q2data, &qq2, &j);
/*<       call findxvh2(xdata,XX,I) >*/
    findxvh2_(xdata, &xx, &i__);
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
/*<       do 21 K=J-1,J+2 >*/
    i__1 = j + 2;
    for (k = j - 1; k <= i__1; ++k) {
/*<          do 21 L=I-1,I+2 >*/
	i__2 = i__ + 2;
	for (l = i__ - 1; l <= i__2; ++l) {
/*<             M = K+2-J >*/
	    m = k + 2 - j;
/*<             N = L+2-I >*/
	    n = l + 2 - i__;
/*<             glh(M,N) = gl(K,L) >*/
	    glh[m + (n << 2) - 5] = partvh2_1.gl[k + l * 34];
/*<             dnh(M,N) = dn(K,L) >*/
	    dnh[m + (n << 2) - 5] = partvh2_1.dn[k + l * 34];
/*<             uph(M,N) = up(K,L) >*/
	    uph[m + (n << 2) - 5] = partvh2_1.up[k + l * 34];
/*<  21         sth(M,N) = st(K,L) >*/
/* L21: */
	    sth[m + (n << 2) - 5] = partvh2_1.st[k + l * 34];
	}
    }
/*<       call fitvh2(bord,QQ2,XX,q2,x,glh,xglu) >*/
    fitvh2_(&bord, &qq2, &xx, q2, x, glh, xglu);
/*<       call fitvh2(bord,QQ2,XX,q2,x,dnh,xdn) >*/
    fitvh2_(&bord, &qq2, &xx, q2, x, dnh, xdn);
/*<       call fitvh2(bord,QQ2,XX,q2,x,uph,xup) >*/
    fitvh2_(&bord, &qq2, &xx, q2, x, uph, xup);
/*<       call fitvh2(bord,QQ2,XX,q2,x,sth,xstr) >*/
    fitvh2_(&bord, &qq2, &xx, q2, x, sth, xstr);
/*<       xglu = alfa*xglu >*/
    *xglu *= .00729735308;
/*<       xdn  = alfa*xdn >*/
    *xdn *= .00729735308;
/*<       xup  = alfa*xup >*/
    *xup *= .00729735308;
/*<       xstr = alfa*xstr >*/
    *xstr *= .00729735308;
/*<       if (IOPT.EQ.1) goto 1000 >*/
    if (*iopt == 1) {
	goto L1000;
    }
/* **************************************************************** */
/* ***   searching for the J such that: Q2(J) < QQ2 < Q2(J+1)  **** */
/* ***    searching for the I such that: x(I) < XX < x(I+1)    **** */
/* ***                  for the heavy quarks                   **** */
/* **************************************************************** */
/*<       NQQ2H = NQ2H >*/
    nqq2h = 48;
/*<       call findq2vh2(NQQ2H,q2hdata,QQ2,JH) >*/
    findq2vh2_(&nqq2h, q2hdata, &qq2, &jh);
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
/*<       if (XX.LT.xmaxc) then >*/
    if (xx < xmaxc) {
/*<          xmxc = 1.d0/(1.d0+4.d0*mc2/q2h(2)) >*/
	xmxc = 1. / (mc2 * 4. / q2h[1] + 1.);
/*<          call findxhvh2(xmxc,xdata,ch,chx,XX,Imaxc(JH),JH,bordc,xc,chh) >*/
	findxhvh2_(&xmxc, xdata, partvh2_1.ch, partvh2_1.chx, &xx, &imaxc[jh 
		- 1], &jh, &bordc, xc, chh);
/*<          call fitvh2(bordc,QQ2,XX,q2h,xc,chh,xchm) >*/
	fitvh2_(&bordc, &qq2, &xx, q2h, xc, chh, xchm);
/*<          xchm = alfa*xchm >*/
	*xchm *= .00729735308;
/*<       endif >*/
    }
/*<       if (XX.LT.xmaxb) then >*/
    if (xx < xmaxb) {
/*<          xmxb = 1.d0/(1.d0+4.d0*mb2/q2h(2)) >*/
	xmxb = 1. / (mb2 * 4. / q2h[1] + 1.);
/*<          call findxhvh2(xmxb,xdata,bt,btx,XX,Imaxb(JH),JH,bordb,xb,bth) >*/
	findxhvh2_(&xmxb, xdata, partvh2_1.bt, partvh2_1.btx, &xx, &imaxb[jh 
		- 1], &jh, &bordb, xb, bth);
/*<          call fitvh2(bordb,QQ2,XX,q2h,xb,bth,xbot) >*/
	fitvh2_(&bordb, &qq2, &xx, q2h, xb, bth, xbot);
/*<          xbot = alfa*xbot >*/
	*xbot *= .00729735308;
/*<       endif >*/
    }
/*<  1000 continue >*/
L1000:
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* gridvh2_ */

/* **************************************************************** */
/*<       SUBROUTINE READTABVH2(IOPT) >*/
/* Subroutine */ int readtabvh2_(integer *iopt)
{
    /* Format strings */
    static char fmt_100[] = "(f11.7,3(2x,f10.7))";
    static char fmt_200[] = "(f10.7,2x,f10.7)";
    static char fmt_300[] = "(23(f10.7,2x),f10.7)";

    /* System generated locals */
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(), f_clos(cllist *);

    /* Local variables */
    static integer i__, j;

    /* Fortran I/O blocks */
    static cilist io___52 = { 0, 10, 0, fmt_100, 0 };
    static cilist io___53 = { 0, 10, 0, fmt_200, 0 };
    static cilist io___54 = { 0, 10, 0, fmt_300, 0 };
    static cilist io___55 = { 0, 10, 0, fmt_300, 0 };


/*<       DOUBLE PRECISION gl,dn,up,st,ch,bt,chx,btx >*/
/*<       parameter (NX=52,NQ2=32,NQ2H=48) >*/
/*<        >*/
/*<       character*15 name >*/
/*<       common /PARTVH2/ gl,dn,up,st,ch,bt,chx,btx >*/
/*<       open(10,file='cjkhobest.dat',status='old') >*/
    o__1.oerr = 0;
    o__1.ounit = 10;
    o__1.ofnmlen = 13;
    o__1.ofnm = "cjkhobest.dat";
    o__1.orl = 0;
    o__1.osta = "old";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    f_open(&o__1);
/*<       do 1 I=0,NQ2+1 >*/
    for (i__ = 0; i__ <= 33; ++i__) {
/*<          do 1 J=0,NX+1 >*/
	for (j = 0; j <= 53; ++j) {
/*<             gl(I,J) = 0.d0 >*/
	    partvh2_1.gl[i__ + j * 34] = 0.;
/*<             dn(I,J) = 0.d0 >*/
	    partvh2_1.dn[i__ + j * 34] = 0.;
/*<             up(I,J) = 0.d0 >*/
	    partvh2_1.up[i__ + j * 34] = 0.;
/*<  1          st(I,J) = 0.d0 >*/
/* L1: */
	    partvh2_1.st[i__ + j * 34] = 0.;
	}
    }
/*<          do 10 I=1,NQ2 >*/
    for (i__ = 1; i__ <= 32; ++i__) {
/*<             do 10 J=1,NX >*/
	for (j = 1; j <= 52; ++j) {
/*<  10    >*/
/* L10: */
	    s_rsfe(&io___52);
	    do_fio(&c__1, (char *)&partvh2_1.gl[i__ + j * 34], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&partvh2_1.dn[i__ + j * 34], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&partvh2_1.up[i__ + j * 34], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&partvh2_1.st[i__ + j * 34], (ftnlen)sizeof(
		    doublereal));
	    e_rsfe();
	}
    }
/*<  100          format (F11.7,3(2X,F10.7)) >*/
/*<          if (IOPT.EQ.1) goto 1000 >*/
    if (*iopt == 1) {
	goto L1000;
    }
/*<          do 2 I=0,NQ2H+1 >*/
    for (i__ = 0; i__ <= 49; ++i__) {
/*<             do 2 J=0,NX+1 >*/
	for (j = 0; j <= 53; ++j) {
/*<                ch(I,J) = 0.d0 >*/
	    partvh2_1.ch[i__ + j * 50] = 0.;
/*<  2             bt(I,J) = 0.d0 >*/
/* L2: */
	    partvh2_1.bt[i__ + j * 50] = 0.;
	}
    }
/*<          do 20 I=1,NQ2H >*/
    for (i__ = 1; i__ <= 48; ++i__) {
/*<             do 20 J=1,NX >*/
	for (j = 1; j <= 52; ++j) {
/*<  20            read(10,200) ch(I,J),bt(I,J) >*/
/* L20: */
	    s_rsfe(&io___53);
	    do_fio(&c__1, (char *)&partvh2_1.ch[i__ + j * 50], (ftnlen)sizeof(
		    doublereal));
	    do_fio(&c__1, (char *)&partvh2_1.bt[i__ + j * 50], (ftnlen)sizeof(
		    doublereal));
	    e_rsfe();
	}
    }
/*<  200           format (F10.7,2X,F10.7) >*/
/*<          do 30 I=1,NQ2H >*/
    for (i__ = 1; i__ <= 48; ++i__) {
/*<             read(10,300) (chx(I,J),J=1,24) >*/
	s_rsfe(&io___54);
	for (j = 1; j <= 24; ++j) {
	    do_fio(&c__1, (char *)&partvh2_1.chx[i__ + j * 50 - 50], (ftnlen)
		    sizeof(doublereal));
	}
	e_rsfe();
/*<  30      continue >*/
/* L30: */
    }
/*<          do 40 I=1,NQ2H >*/
    for (i__ = 1; i__ <= 48; ++i__) {
/*<             read(10,300) (btx(I,J),J=1,24) >*/
	s_rsfe(&io___55);
	for (j = 1; j <= 24; ++j) {
	    do_fio(&c__1, (char *)&partvh2_1.btx[i__ + j * 50 - 50], (ftnlen)
		    sizeof(doublereal));
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
/*<                chx(I,J) = 0.d0 >*/
	    partvh2_1.chx[i__ + j * 50 - 50] = 0.;
/*<  31            btx(I,J) = 0.d0 >*/
/* L31: */
	    partvh2_1.btx[i__ + j * 50 - 50] = 0.;
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
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* readtabvh2_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDQ2VH2(NQ2,q2data,QQ2,I) >*/
/* Subroutine */ int findq2vh2_(integer *nq2, doublereal *q2data, doublereal *
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
} /* findq2vh2_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDXVH2(xdata,XX,I) >*/
/* Subroutine */ int findxvh2_(doublereal *xdata, doublereal *xx, integer *
	i__)
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
} /* findxvh2_ */

/* ******************************************************************** */
/* ***            Here I use the bisection method                  **** */
/* ******************************************************************** */
/*<       SUBROUTINE FINDXHVH2(xmax,xdata,hq,hqx,XX,Imax,J,bord,x,hqh) >*/
/* Subroutine */ int findxhvh2_(doublereal *xmax, doublereal *xdata, 
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
/*<       INTEGER I,J,K,L,M,N,il,iu,Imax,bord,NX,NQ2 >*/
/*<       parameter (NX=52,NQ2=48) >*/
/*<        >*/
/*<       data per/0.93d0,0.95d0,0.96d0,0.97d0,0.99d0,0.999d0,0d0/ >*/
    /* Parameter adjustments */
    hqh -= 5;
    --x;
    hqx -= 50;

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
/*<  30               hqh(M+2-J,L+2-I) = hq(M,L) >*/
/* L30: */
		    hqh[m + 2 - *j + (l + 2 - i__ << 2)] = hq[m + l * 50];
		}
	    }
/*<          else >*/
	} else {
/*<             do 40 M=J-1,J+2 >*/
	    i__2 = *j + 2;
	    for (m = *j - 1; m <= i__2; ++m) {
/*<                N = M+2-J >*/
		n = m + 2 - *j;
/*<                hqh(N,4) = hqx(M,6*N-5) >*/
		hqh[n + 16] = hqx[m + (n * 6 - 5) * 50];
/*<                do 40 L=I-1,I+1 >*/
		i__1 = i__ + 1;
		for (l = i__ - 1; l <= i__1; ++l) {
/*<  40               hqh(N,L+2-I) = hq(M,L) >*/
/* L40: */
		    hqh[n + (l + 2 - i__ << 2)] = hq[m + l * 50];
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
/*<             hqh(N,1) = hq(M,Imax-1) >*/
	    hqh[n + 4] = hq[m + (*imax - 1) * 50];
/*<             hqh(N,2) = hq(M,Imax) >*/
	    hqh[n + 8] = hq[m + *imax * 50];
/*<             hqh(N,3) = hqx(M,6*N-5) >*/
	    hqh[n + 12] = hqx[m + (n * 6 - 5) * 50];
/*<  50         hqh(N,4) = hqx(M,6*N-4) >*/
/* L50: */
	    hqh[n + 16] = hqx[m + (n * 6 - 4) * 50];
	}
/*<       elseif (K.EQ.2) then >*/
    } else if (k == 2) {
/*<          do 60 M=J-1,J+2 >*/
	i__1 = *j + 2;
	for (m = *j - 1; m <= i__1; ++m) {
/*<             N = M+2-J >*/
	    n = m + 2 - *j;
/*<             hqh(N,1) = hq(M,Imax) >*/
	    hqh[n + 4] = hq[m + *imax * 50];
/*<             do 60 I=2,4 >*/
	    for (i__ = 2; i__ <= 4; ++i__) {
/*<  60            hqh(N,I) = hqx(M,6*N+I-7) >*/
/* L60: */
		hqh[n + (i__ << 2)] = hqx[m + (n * 6 + i__ - 7) * 50];
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
/*<  70            hqh(N,I) = hqx(M,6*N+I-6) >*/
/* L70: */
		hqh[n + (i__ << 2)] = hqx[m + (n * 6 + i__ - 6) * 50];
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
/*<  80            hqh(N,I) = hqx(M,6*N+I-5) >*/
/* L80: */
		hqh[n + (i__ << 2)] = hqx[m + (n * 6 + i__ - 5) * 50];
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
/*<  90            hqh(N,I) = hqx(M,6*N+I-4) >*/
/* L90: */
		hqh[n + (i__ << 2)] = hqx[m + (n * 6 + i__ - 4) * 50];
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
/*<  100           hqh(N,I) = hqx(M,6*N+I-3) >*/
/* L100: */
		hqh[n + (i__ << 2)] = hqx[m + (n * 6 + i__ - 3) * 50];
	    }
	}
/*<       endif >*/
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* findxhvh2_ */

/* ******************************************************************** */
/* ***         Here I use the bicubic interpolation in             **** */
/* ***              the Hermit polynomials basis                   **** */
/* ******************************************************************** */
/*<       SUBROUTINE FITVH2(bord,xx1,xx2,x1,x2,yg,result) >*/
/* Subroutine */ int fitvh2_(integer *bord, doublereal *xx1, doublereal *xx2, 
	doublereal *x1, doublereal *x2, doublereal *yg, doublereal *result)
{
    static doublereal y[4], y1[4], y2[4], y12[4], x1l, x2l, x1u, x2u;
    extern doublereal d1fvh2_(doublereal *, doublereal *, doublereal *, 
	    doublereal *), d2fvh2_(doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    extern /* Subroutine */ int itervh2_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       dimension x1(4),x2(4),yg(4,4),y(4),y1(4),y2(4),y12(4) >*/
/*<       integer bord >*/
/*<       external d1fvh2,d2fvh2 >*/
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
/*<          y1(1) = d1fvh2(x1(1),x1(3),yg(1,2),yg(3,2)) >*/
	y1[0] = d1fvh2_(&x1[1], &x1[3], &yg[9], &yg[11]);
/*<          y1(2) = d1fvh2(x1(2),x1(4),yg(2,2),yg(4,2)) >*/
	y1[1] = d1fvh2_(&x1[2], &x1[4], &yg[10], &yg[12]);
/*<          y1(3) = d1fvh2(x1(2),x1(4),yg(2,3),yg(4,3)) >*/
	y1[2] = d1fvh2_(&x1[2], &x1[4], &yg[14], &yg[16]);
/*<          y1(4) = d1fvh2(x1(1),x1(3),yg(1,3),yg(3,3)) >*/
	y1[3] = d1fvh2_(&x1[1], &x1[3], &yg[13], &yg[15]);
/*<          y2(1) = d1fvh2(x2(1),x2(3),yg(2,1),yg(2,3)) >*/
	y2[0] = d1fvh2_(&x2[1], &x2[3], &yg[6], &yg[14]);
/*<          y2(2) = d1fvh2(x2(1),x2(3),yg(3,1),yg(3,3)) >*/
	y2[1] = d1fvh2_(&x2[1], &x2[3], &yg[7], &yg[15]);
/*<          y2(3) = d1fvh2(x2(2),x2(4),yg(3,2),yg(3,4)) >*/
	y2[2] = d1fvh2_(&x2[2], &x2[4], &yg[11], &yg[19]);
/*<          y2(4) = d1fvh2(x2(2),x2(4),yg(2,2),yg(2,4)) >*/
	y2[3] = d1fvh2_(&x2[2], &x2[4], &yg[10], &yg[18]);
/*<        >*/
	y12[0] = d2fvh2_(&x1[1], &x1[3], &x2[1], &x2[3], &yg[5], &yg[13], &yg[
		7], &yg[15]);
/*<        >*/
	y12[1] = d2fvh2_(&x1[2], &x1[4], &x2[1], &x2[3], &yg[6], &yg[14], &yg[
		8], &yg[16]);
/*<        >*/
	y12[2] = d2fvh2_(&x1[2], &x1[4], &x2[2], &x2[4], &yg[10], &yg[18], &
		yg[12], &yg[20]);
/*<        >*/
	y12[3] = d2fvh2_(&x1[1], &x1[3], &x2[2], &x2[4], &yg[9], &yg[17], &yg[
		11], &yg[19]);
/*<       else >*/
    } else {
/*<          y1(1) = d1fvh2(x1(2),x1(3),yg(2,2),yg(3,2)) >*/
	y1[0] = d1fvh2_(&x1[2], &x1[3], &yg[10], &yg[11]);
/*<          y1(2) = y1(1) >*/
	y1[1] = y1[0];
/*<          y1(3) = d1fvh2(x1(2),x1(3),yg(2,3),yg(3,3)) >*/
	y1[2] = d1fvh2_(&x1[2], &x1[3], &yg[14], &yg[15]);
/*<          y1(4) = y1(3) >*/
	y1[3] = y1[2];
/*<          y2(1) = d1fvh2(x2(2),x2(3),yg(2,2),yg(2,3)) >*/
	y2[0] = d1fvh2_(&x2[2], &x2[3], &yg[10], &yg[14]);
/*<          y2(2) = d1fvh2(x2(2),x2(3),yg(3,2),yg(3,3)) >*/
	y2[1] = d1fvh2_(&x2[2], &x2[3], &yg[11], &yg[15]);
/*<          y2(3) = y2(2) >*/
	y2[2] = y2[1];
/*<          y2(4) = y2(1) >*/
	y2[3] = y2[0];
/*<        >*/
	y12[0] = d2fvh2_(&x1[2], &x1[3], &x2[2], &x2[3], &yg[10], &yg[14], &
		yg[11], &yg[15]);
/*<          y12(2) = y12(1) >*/
	y12[1] = y12[0];
/*<          y12(3) = y12(1) >*/
	y12[2] = y12[0];
/*<          y12(4) = y12(1) >*/
	y12[3] = y12[0];
/*<       endif >*/
    }
/*<       call itervh2(y,y1,y2,y12,x1l,x1u,x2l,x2u,xx1,xx2,result) >*/
    itervh2_(y, y1, y2, y12, &x1l, &x1u, &x2l, &x2u, xx1, xx2, result);
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* fitvh2_ */

/* ******************************************************************** */
/*<       SUBROUTINE itervh2(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,res) >*/
/* Subroutine */ int itervh2_(doublereal *y, doublereal *y1, doublereal *y2, 
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
} /* itervh2_ */

/* ******************************************************************** */
/* ***                   First derivative: df/dx                   **** */
/* ******************************************************************** */
/*<       DOUBLE PRECISION FUNCTION d1fvh2(x1,x2,f1,f2) >*/
doublereal d1fvh2_(doublereal *x1, doublereal *x2, doublereal *f1, doublereal 
	*f2)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       d1fvh2 = (f2-f1)/(x2-x1) >*/
    ret_val = (*f2 - *f1) / (*x2 - *x1);
/*<       END >*/
    return ret_val;
} /* d1fvh2_ */

/* ******************************************************************** */
/* ***                Second derivative: d2f/dxdy                  **** */
/* ******************************************************************** */
/*<       DOUBLE PRECISION FUNCTION d2fvh2(x1,x2,y1,y2,f11,f12,f21,f22) >*/
doublereal d2fvh2_(doublereal *x1, doublereal *x2, doublereal *y1, doublereal 
	*y2, doublereal *f11, doublereal *f12, doublereal *f21, doublereal *
	f22)
{
    /* System generated locals */
    doublereal ret_val;

/*<       IMPLICIT DOUBLE PRECISION (a-z) >*/
/*<       d2fvh2 = (f22 - f21 - f12 + f11)/((x2-x1)*(y2-y1)) >*/
    ret_val = (*f22 - *f21 - *f12 + *f11) / ((*x2 - *x1) * (*y2 - *y1));
/*<       END >*/
    return ret_val;
} /* d2fvh2_ */

#ifdef __cplusplus
	}
#endif
