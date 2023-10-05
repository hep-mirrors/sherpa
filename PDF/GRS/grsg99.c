/* grsg99.f -- translated by f2c (version 20200916).
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

struct {
    integer iini;
} intini_;

#define intini_1 intini_

/* Table of constant values */

static integer c__1 = 1;
static doublereal c_b24 = .7;
static doublereal c_b25 = .3;
static integer c__9 = 9;
static integer c__2 = 2;

/* ******************************************************************** */
/*                                                                   * */
/*            LO AND NLO (REAL) PHOTON - PARAMETRIZATIONS            * */
/*                                                                   * */
/*                 FOR A DETAILED EXPLANATION SEE :                  * */
/*                 M. GLUECK, E.REYA, I. Schienbein                  * */
/*                    Phys. Rev. D60(1999)054019                     * */
/*                                                                   * */
/*        PROBLEMS/QUESTIONS TO schien@lpsc.in2p3.fr                 * */
/*                                                                   * */
/*   INPUT:   ISET = number of the parton set :                      * */
/*              ISET = 1  LEADING ORDER SET                          * */
/*                        (DATA FILE 'grsg99lo.grid')                * */
/*              ISET = 2  NEXT-TO-LEADING ORDER MSbar SET            * */
/*                        (DATA FILE 'grsg99m.grid')                 * */
/*              ISET = 3  NEXT-TO-LEADING ORDER DISgamma SET         * */
/*                        (DATA FILE 'grsg99d.grid')                 * */
/*                                                                   * */
/*            X  = Bjorken-x       (between  1.E-5 and 1   )         * */
/*            Q2 = scale in GeV**2 (between  0.4   and 1.E6)         * */
/*                                                                   * */
/*   OUTPUT:  UL = X*UP/ALPHA ; DL = x*DOWN/ALPHA                    * */
/*            SL = X*STRANGE SEA/ALPHA : GL = X*GLUON/ALPHA          * */
/*                                                                   * */
/*            Always x times the distribution is returned            * */
/*                                                                   * */
/*   COMMON:  The main program or the calling routine has to have    * */
/*            a common block  COMMON / INTINI / IINI , and  IINI     * */
/*            has always to be zero when GRSG99 is called for the    * */
/*            first time or when 'ISET' has been changed.            * */
/*                                                                   * */
/* 16.01.2001                                                        * */
/* ******************************************************************** */
/* ... */
/*<       SUBROUTINE GRSG99(ISET,X,Q2,UL,DL,SL,GL) >*/
/* Subroutine */ int grsg99_(integer *iset, doublereal *x, doublereal *q2, 
	doublereal *ul, doublereal *dl, doublereal *sl, doublereal *gl)
{
    /* Initialized data */

    static doublereal qs[34] = { .4,.45,.5,.55,.6,.75,1.,1.5,2.,2.5,4.,6.4,
	    10.,15.,25.,40.,64.,100.,180.,320.,580.,1e3,1800.,3200.,5800.,1e4,
	    1.8e4,3.2e4,5.8e4,1e5,1.8e5,3.2e5,5.8e5,1e6 };
    static doublereal xb[51] = { 1e-5,1.5e-5,2.2e-5,3.2e-5,4.8e-5,7e-5,1e-4,
	    1.5e-4,2.2e-4,3.2e-4,4.8e-4,7e-4,.001,.0015,.0022,.0032,.0048,
	    .007,.01,.015,.022,.032,.05,.075,.1,.125,.15,.175,.2,.225,.25,
	    .275,.3,.325,.35,.375,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.92,
	    .95,.98,.99 };

    /* Format strings */
    static char fmt_91[] = "(2x,\002PARTON INTERPOLATION: X OUT OF RANGE\002)"
	    ;
    static char fmt_92[] = "(2x,\002PARTON INTERPOLATION: Q2 OUT OF RANGE\
\002)";
    static char fmt_93[] = "(2x,\002PARTON INTERPOLATION: ISET OUT OF RANG\
E\002)";
    static char fmt_90[] = "(4(1pe12.5))";

    /* System generated locals */
    doublereal d__1, d__2;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer f_open(olist *), s_rsfe(cilist *), do_fio(integer *, char *, 
	    ftnlen), e_rsfe(), f_clos(cllist *);
    double pow_dd(doublereal *, doublereal *), log(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();

    /* Local variables */
    static integer m, n, na[2], iq, ix;
    static doublereal xt[2], xb0, xb1, xdf[1734]	/* was [51][34] */, 
	    xgf[1734]	/* was [51][34] */, xsf[1734]	/* was [51][34] */, 
	    xuf[1734]	/* was [51][34] */, arrf[85];
    extern doublereal fint_(integer *, doublereal *, integer *, doublereal *, 
	    doublereal *);
    static integer iiread;
    static doublereal parton[6800]	/* was [4][34][50] */;
    static integer isetsav;

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, fmt_91, 0 };
    static cilist io___4 = { 0, 6, 0, fmt_92, 0 };
    static cilist io___7 = { 0, 6, 0, fmt_93, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_90, 0 };
    static cilist io___22 = { 0, 6, 0, 0, 0 };
    static cilist io___23 = { 0, 6, 0, 0, 0 };


/* ... */
/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       INTEGER ISET,IINI,IIREAD, ISETSAV  >*/
/*<       PARAMETER (NPART=4, NX=51, NQ=34, NARG=2) >*/
/*<       DIMENSION XUF(NX,NQ), XDF(NX,NQ), XSF(NX,NQ), XGF(NX,NQ) >*/
/*<       DIMENSION PARTON (NPART,NQ,NX-1), QS(NQ), XB(NX) >*/
/*<       DIMENSION XT(NARG), NA(NARG), ARRF(NX+NQ) >*/
/*<       COMMON / INTINI / IINI >*/
/*<       SAVE XUF, XDF, XSF, XGF, NA, ARRF, ISETSAV >*/
/* ...BJORKEN-X AND Q**2 VALUES OF THE GRID : */
/*<        >*/
/*<        >*/
/* ...CHECK OF X AND Q2 VALUES : */
/*<        IF ( (X.LT.1.0D-5) .OR. (X.GT.1.0D0) ) THEN >*/
    if (*x < 1e-5 || *x > 1.) {
/*<            WRITE(6,91) >*/
	s_wsfe(&io___3);
	e_wsfe();
/*<   91       FORMAT (2X,'PARTON INTERPOLATION: X OUT OF RANGE') >*/
/*<            STOP >*/
	s_stop("", (ftnlen)0);
/*<        ENDIF >*/
    }
/*<        IF ( (Q2.LT.0.4D0) .OR. (Q2.GT.1.D6) ) THEN >*/
    if (*q2 < .4 || *q2 > 1e6) {
/*<            WRITE(6,92) >*/
	s_wsfe(&io___4);
	e_wsfe();
/*<   92       FORMAT (2X,'PARTON INTERPOLATION: Q2 OUT OF RANGE') >*/
/*<            STOP >*/
	s_stop("", (ftnlen)0);
/*<        ENDIF >*/
    }
/* ...INITIALIZATION : */
/*    SELECTION AND READING OF THE GRID : */
/*    FILE - NO. = 11 FOR LEADING ORDER ( FIRST NUMBER IN THE */
/*                                        GRID: 1.884E-03 ) */
/*    FILE - NO. = 22 FOR NEXT-TO-LEADING ORDER MSbar ( FIRST NUMBER IN THE */
/*                                                      GRID: 1.658E-03 ) */
/*    FILE - NO. = 33 FOR NEXT-TO-LEADING ORDER DISgamma ( FIRST NUMBER IN THE */
/*                                                         GRID: 1.658E-03 ) */
/*<       IF (IINI.NE.0) GOTO 16 >*/
    if (intini_1.iini != 0) {
	goto L16;
    }
/*<       IF (ISET.EQ.1) THEN >*/
    if (*iset == 1) {
/*<        ISETSAV=1   >*/
	isetsav = 1;
/*<        IIREAD=11 >*/
	iiread = 11;
/*<        OPEN(11,FILE='grsg99lo.grid') >*/
	o__1.oerr = 0;
	o__1.ounit = 11;
	o__1.ofnmlen = 13;
	o__1.ofnm = "grsg99lo.grid";
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
/*<       ELSE IF (ISET.EQ.2) THEN >*/
    } else if (*iset == 2) {
/*<        ISETSAV=2   >*/
	isetsav = 2;
/*<        IIREAD=22 >*/
	iiread = 22;
/*<        OPEN(22,FILE='grsg99m.grid') >*/
	o__1.oerr = 0;
	o__1.ounit = 22;
	o__1.ofnmlen = 12;
	o__1.ofnm = "grsg99m.grid";
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
/*<       ELSE IF (ISET.EQ.3) THEN >*/
    } else if (*iset == 3) {
/*<        ISETSAV=3   >*/
	isetsav = 3;
/*<        IIREAD=33 >*/
	iiread = 33;
/*<        OPEN(33,FILE='grsg99d.grid') >*/
	o__1.oerr = 0;
	o__1.ounit = 33;
	o__1.ofnmlen = 12;
	o__1.ofnm = "grsg99d.grid";
	o__1.orl = 0;
	o__1.osta = 0;
	o__1.oacc = 0;
	o__1.ofm = 0;
	o__1.oblnk = 0;
	f_open(&o__1);
/*<       ELSE >*/
    } else {
/*<         WRITE(6,93) >*/
	s_wsfe(&io___7);
	e_wsfe();
/*<   93    FORMAT (2X,'PARTON INTERPOLATION: ISET OUT OF RANGE') >*/
/*<         GOTO 60 >*/
	goto L60;
/*<       END IF >*/
    }

/*<        DO 15 M = 1, NX-1 >*/
    for (m = 1; m <= 50; ++m) {
/*<        DO 15 N = 1, NQ >*/
	for (n = 1; n <= 34; ++n) {
/*<        >*/
	    io___10.ciunit = iiread;
	    s_rsfe(&io___10);
	    do_fio(&c__1, (char *)&parton[(n + m * 34 << 2) - 140], (ftnlen)
		    sizeof(doublereal));
	    do_fio(&c__1, (char *)&parton[(n + m * 34 << 2) - 139], (ftnlen)
		    sizeof(doublereal));
	    do_fio(&c__1, (char *)&parton[(n + m * 34 << 2) - 138], (ftnlen)
		    sizeof(doublereal));
	    do_fio(&c__1, (char *)&parton[(n + m * 34 << 2) - 137], (ftnlen)
		    sizeof(doublereal));
	    e_rsfe();
/*<   90   FORMAT (4(1PE12.5)) >*/
/*<   15   CONTINUE >*/
/* L15: */
	}
    }

/*<       close(IIREAD) >*/
    cl__1.cerr = 0;
    cl__1.cunit = iiread;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*<       IINI = 1 >*/
    intini_1.iini = 1;
/* ....ARRAYS FOR THE INTERPOLATION SUBROUTINE : */
/*<       DO 10 IQ = 1, NQ >*/
    for (iq = 1; iq <= 34; ++iq) {
/*<       DO 20 IX = 1, NX-1 >*/
	for (ix = 1; ix <= 50; ++ix) {
/*<         XB0 = XB(IX) >*/
	    xb0 = xb[ix - 1];
/*<         XB1 = 1.D0-XB(IX) >*/
	    xb1 = 1. - xb[ix - 1];
/*<         XUF(IX,IQ) = PARTON(1,IQ,IX) / (XB1**3 * XB0**0.7) >*/
/* Computing 3rd power */
	    d__1 = xb1;
	    xuf[ix + iq * 51 - 52] = parton[(iq + ix * 34 << 2) - 140] / (
		    d__1 * (d__1 * d__1) * pow_dd(&xb0, &c_b24));
/*<         XDF(IX,IQ) = PARTON(2,IQ,IX) / (XB1**7 * XB0**0.3) >*/
/* Computing 7th power */
	    d__1 = xb1, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	    xdf[ix + iq * 51 - 52] = parton[(iq + ix * 34 << 2) - 139] / (
		    d__2 * (d__1 * d__1) * pow_dd(&xb0, &c_b25));
/*<         XSF(IX,IQ) = PARTON(3,IQ,IX) / (XB1**7 * XB0**0.3) >*/
/* Computing 7th power */
	    d__1 = xb1, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
	    xsf[ix + iq * 51 - 52] = parton[(iq + ix * 34 << 2) - 138] / (
		    d__2 * (d__1 * d__1) * pow_dd(&xb0, &c_b25));
/*<         XGF(IX,IQ) = PARTON(4,IQ,IX) / (XB1**5 * XB0**0.3) >*/
/* Computing 5th power */
	    d__1 = xb1, d__2 = d__1, d__1 *= d__1;
	    xgf[ix + iq * 51 - 52] = parton[(iq + ix * 34 << 2) - 137] / (
		    d__2 * (d__1 * d__1) * pow_dd(&xb0, &c_b25));
/*<   20  CONTINUE >*/
/* L20: */
	}
/*<         XUF(NX,IQ) = 0.D0 >*/
	xuf[iq * 51 - 1] = 0.;
/*<         XDF(NX,IQ) = 0.D0 >*/
	xdf[iq * 51 - 1] = 0.;
/*<         XSF(NX,IQ) = 0.D0 >*/
	xsf[iq * 51 - 1] = 0.;
/*<         XGF(NX,IQ) = 0.D0 >*/
	xgf[iq * 51 - 1] = 0.;
/*<   10  CONTINUE >*/
/* L10: */
    }
/*<       NA(1) = NX >*/
    na[0] = 51;
/*<       NA(2) = NQ >*/
    na[1] = 34;
/*<       DO 30 IX = 1, NX >*/
    for (ix = 1; ix <= 51; ++ix) {
/*<         ARRF(IX) = DLOG(XB(IX)) >*/
	arrf[ix - 1] = log(xb[ix - 1]);
/*<   30  CONTINUE >*/
/* L30: */
    }
/*<       DO 40 IQ = 1, NQ >*/
    for (iq = 1; iq <= 34; ++iq) {
/*<         ARRF(NX+IQ) = DLOG(QS(IQ)) >*/
	arrf[iq + 50] = log(qs[iq - 1]);
/*<   40  CONTINUE >*/
/* L40: */
    }
/*<   16  CONTINUE >*/
L16:
/*<       if (ISET .ne. ISETSAV) then  >*/
    if (*iset != isetsav) {
/*<          print*,'Warning : ISET has been changed' >*/
	s_wsle(&io___22);
	do_lio(&c__9, &c__1, "Warning : ISET has been changed", (ftnlen)31);
	e_wsle();
/*<          print*,'You should reinitialize the GRIDS ! by setting iini=0' >*/
	s_wsle(&io___23);
	do_lio(&c__9, &c__1, "You should reinitialize the GRIDS ! by setting\
 iini=0", (ftnlen)53);
	e_wsle();
/*<       end if    >*/
    }
/* ...INTERPOLATION : */
/*<       XT(1) = DLOG(X) >*/
    xt[0] = log(*x);
/*<       XT(2) = DLOG(Q2) >*/
    xt[1] = log(*q2);
/*<       UL = FINT(NARG,XT,NA,ARRF,XUF) * (1.D0-X)**3 * X**0.7 >*/
/* Computing 3rd power */
    d__1 = 1. - *x;
    *ul = fint_(&c__2, xt, na, arrf, xuf) * (d__1 * (d__1 * d__1)) * pow_dd(x,
	     &c_b24);
/*<       DL = FINT(NARG,XT,NA,ARRF,XDF) * (1.D0-X)**7 * X**0.3 >*/
/* Computing 7th power */
    d__1 = 1. - *x, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
    *dl = fint_(&c__2, xt, na, arrf, xdf) * (d__2 * (d__1 * d__1)) * pow_dd(x,
	     &c_b25);
/*<       SL = FINT(NARG,XT,NA,ARRF,XSF) * (1.D0-X)**7 * X**0.3 >*/
/* Computing 7th power */
    d__1 = 1. - *x, d__2 = d__1, d__1 *= d__1, d__2 *= d__1;
    *sl = fint_(&c__2, xt, na, arrf, xsf) * (d__2 * (d__1 * d__1)) * pow_dd(x,
	     &c_b25);
/*<       GL = FINT(NARG,XT,NA,ARRF,XGF) * (1.D0-X)**5 * X**0.3 >*/
/* Computing 5th power */
    d__1 = 1. - *x, d__2 = d__1, d__1 *= d__1;
    *gl = fint_(&c__2, xt, na, arrf, xgf) * (d__2 * (d__1 * d__1)) * pow_dd(x,
	     &c_b25);
/*<  60   RETURN >*/
L60:
    return 0;
/*<       END >*/
} /* grsg99_ */


/* ...CERN LIBRARY ROUTINE E104 (INTERPOLATION) : */

/*<       FUNCTION FINT(NARG,ARG,NENT,ENT,TABLE) >*/
doublereal fint_(integer *narg, doublereal *arg, integer *nent, doublereal *
	ent, doublereal *table)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal ret_val;

    /* Local variables */
    static doublereal d__[5];
    static integer i__, j, k, m, ja, jb, kd, il, jr;
    static doublereal fac;
    static integer iadr, ient[5], ifadr, ncomb[5];

/*<       IMPLICIT DOUBLE PRECISION (A-H,O-Z) >*/
/*<       DIMENSION ARG(5),NENT(5),ENT(63),TABLE(882) >*/
/*<       DIMENSION D(5),NCOMB(5),IENT(5) >*/
/*<       KD=1 >*/
    /* Parameter adjustments */
    --table;
    --ent;
    --nent;
    --arg;

    /* Function Body */
    kd = 1;
/*<       M=1 >*/
    m = 1;
/*<       JA=1 >*/
    ja = 1;
/*<          DO 5 I=1,NARG >*/
    i__1 = *narg;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       NCOMB(I)=1 >*/
	ncomb[i__ - 1] = 1;
/*<       JB=JA-1+NENT(I) >*/
	jb = ja - 1 + nent[i__];
/*<          DO 2 J=JA,JB >*/
	i__2 = jb;
	for (j = ja; j <= i__2; ++j) {
/*<       IF (ARG(I).LE.ENT(J)) GO TO 3 >*/
	    if (arg[i__] <= ent[j]) {
		goto L3;
	    }
/*<     2 CONTINUE >*/
/* L2: */
	}
/*<       J=JB >*/
	j = jb;
/*<     3 IF (J.NE.JA) GO TO 4 >*/
L3:
	if (j != ja) {
	    goto L4;
	}
/*<       J=J+1 >*/
	++j;
/*<     4 JR=J-1 >*/
L4:
	jr = j - 1;
/*<       D(I)=(ENT(J)-ARG(I))/(ENT(J)-ENT(JR)) >*/
	d__[i__ - 1] = (ent[j] - arg[i__]) / (ent[j] - ent[jr]);
/*<       IENT(I)=J-JA >*/
	ient[i__ - 1] = j - ja;
/*<       KD=KD+IENT(I)*M >*/
	kd += ient[i__ - 1] * m;
/*<       M=M*NENT(I) >*/
	m *= nent[i__];
/*<     5 JA=JB+1 >*/
/* L5: */
	ja = jb + 1;
    }
/*<       FINT=0.D0 >*/
    ret_val = 0.;
/*<    10 FAC=1.D0 >*/
L10:
    fac = 1.;
/*<       IADR=KD >*/
    iadr = kd;
/*<       IFADR=1 >*/
    ifadr = 1;
/*<          DO 15 I=1,NARG >*/
    i__1 = *narg;
    for (i__ = 1; i__ <= i__1; ++i__) {
/*<       IF (NCOMB(I).EQ.0) GO TO 12 >*/
	if (ncomb[i__ - 1] == 0) {
	    goto L12;
	}
/*<       FAC=FAC*(1.D0-D(I)) >*/
	fac *= 1. - d__[i__ - 1];
/*<       GO TO 15 >*/
	goto L15;
/*<    12 FAC=FAC*D(I) >*/
L12:
	fac *= d__[i__ - 1];
/*<       IADR=IADR-IFADR >*/
	iadr -= ifadr;
/*<    15 IFADR=IFADR*NENT(I) >*/
L15:
	ifadr *= nent[i__];
    }
/*<       FINT=FINT+FAC*TABLE(IADR) >*/
    ret_val += fac * table[iadr];
/*<       IL=NARG >*/
    il = *narg;
/*<    40 IF (NCOMB(IL).EQ.0) GO TO 80 >*/
L40:
    if (ncomb[il - 1] == 0) {
	goto L80;
    }
/*<       NCOMB(IL)=0 >*/
    ncomb[il - 1] = 0;
/*<       IF (IL.EQ.NARG) GO TO 10 >*/
    if (il == *narg) {
	goto L10;
    }
/*<       IL=IL+1 >*/
    ++il;
/*<          DO 50  K=IL,NARG >*/
    i__1 = *narg;
    for (k = il; k <= i__1; ++k) {
/*<    50 NCOMB(K)=1 >*/
/* L50: */
	ncomb[k - 1] = 1;
    }
/*<       GO TO 10 >*/
    goto L10;
/*<    80 IL=IL-1 >*/
L80:
    --il;
/*<       IF(IL.NE.0) GO TO 40 >*/
    if (il != 0) {
	goto L40;
    }
/*<       RETURN >*/
    return ret_val;
/*<       END >*/
} /* fint_ */

#ifdef __cplusplus
	}
#endif
