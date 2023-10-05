/* sasgam2.f -- translated by f2c (version 20200916).
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
    real xpvmd[13], xpanl[13], xpanh[13], xpbeh[13], xpdir[13];
} sascom_;

#define sascom_1 sascom_

struct {
    real vxpvmd[13], vxpanl[13], vxpanh[13], vxpdgm[13];
} sasval_;

#define sasval_1 sasval_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c_n3 = -3;
static integer c__5 = 5;
static integer c__0 = 0;
static doublereal c_b39 = .07407407407407407;
static doublereal c_b40 = .086956521739130432;
static doublereal c_b41 = .8;
static doublereal c_b42 = .76;
static doublereal c_b43 = .4;
static doublereal c_b44 = 1.76;
static doublereal c_b45 = 3.76;
static doublereal c_b46 = 1.3;
static doublereal c_b49 = .51;
static doublereal c_b50 = 1.37;
static doublereal c_b51 = .255;
static doublereal c_b52 = 2.37;
static doublereal c_b55 = .46;
static doublereal c_b56 = .64;
static doublereal c_b57 = .5;
static doublereal c_b58 = 2.6;

/* ...SaSgam version 2 - parton distributions of the photon */
/* ...by Gerhard A. Schuler and Torbjorn Sjostrand */
/* ...For further information see Z. Phys. C68 (1995) 607 */
/* ...and Phys. Lett. B376 (1996) 193. */
/* ...18 January 1996: original code. */
/* ...22 July 1996: calculation of BETA moved in SASBEH. */
/* !!!Note that one further call parameter - IP2 - has been added */
/* !!!to the SASGAM argument list compared with version 1. */
/* ...The user should only need to call the SASGAM routine, */
/* ...which in turn calls the auxiliary routines SASVMD, SASANO, */
/* ...SASBEH and SASDIR. The package is self-contained. */
/* ...One particular aspect of these parametrizations is that F2 for */
/* ...the photon is not obtained just as the charge-squared-weighted */
/* ...sum of quark distributions, but differ in the treatment of */
/* ...heavy flavours (in F2 the DIS relation W2 = Q2*(1-x)/x restricts */
/* ...the kinematics range of heavy-flavour production, but the same */
/* ...kinematics is not relevant e.g. for jet production) and, for the */
/* ...'MSbar' fits, in the addition of a Cgamma term related to the */
/* ...separation of direct processes. Schematically: */
/* ...PDF = VMD (rho, omega, phi) + anomalous (d, u, s, c, b). */
/* ...F2  = VMD (rho, omega, phi) + anomalous (d, u, s) + */
/* ...      Bethe-Heitler (c, b) (+ Cgamma (d, u, s)). */
/* ...The J/psi and Upsilon states have not been included in the VMD sum, */
/* ...but low c and b masses in the other components should compensate */
/* ...for this in a duality sense. */
/* ...The calling sequence is the following: */
/*     CALL SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM) */
/* ...with the following declaration statement: */
/*     DIMENSION XPDFGM(-6:6) */
/* ...and, optionally, further information in: */
/*     COMMON/SASCOM/XPVMD(-6:6),XPANL(-6:6),XPANH(-6:6),XPBEH(-6:6), */
/*    &XPDIR(-6:6) */
/*     COMMON/SASVAL/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6) */
/* ...Input:  ISET = 1 : SaS set 1D ('DIS',   Q0 = 0.6 GeV) */
/*                = 2 : SaS set 1M ('MSbar', Q0 = 0.6 GeV) */
/*                = 3 : SaS set 2D ('DIS',   Q0 =  2  GeV) */
/*                = 4 : SaS set 2M ('MSbar', Q0 =  2  GeV) */
/*           X : x value. */
/*           Q2 : Q2 value. */
/*           P2 : P2 value; should be = 0. for an on-shell photon. */
/*           IP2 : scheme used to evaluate off-shell anomalous component. */
/*               = 0 : recommended default, see = 7. */
/*               = 1 : dipole dampening by integration; very time-consuming. */
/*               = 2 : P_0^2 = max( Q_0^2, P^2 ) */
/*               = 3 : P'_0^2 = Q_0^2 + P^2. */
/*               = 4 : P_{eff} that preserves momentum sum. */
/*               = 5 : P_{int} that preserves momentum and average */
/*                     evolution range. */
/*               = 6 : P_{eff}, matched to P_0 in P2 -> Q2 limit. */
/*               = 7 : P_{int}, matched to P_0 in P2 -> Q2 limit. */
/* ...Output: F2GM : F2 value of the photon (including factors of alpha_em). */
/*           XPFDGM :  x times parton distribution functions of the photon, */
/*               with elements 0 = g, 1 = d, 2 = u, 3 = s, 4 = c, 5 = b, */
/*               6 = t (always empty!), - for antiquarks (result is same). */
/* ...The breakdown by component is stored in the commonblock SASCOM, */
/*               with elements as above. */
/*           XPVMD : rho, omega, phi VMD part only of output. */
/*           XPANL : d, u, s anomalous part only of output. */
/*           XPANH : c, b anomalous part only of output. */
/*           XPBEH : c, b Bethe-Heitler part only of output. */
/*           XPDIR : Cgamma (direct contribution) part only of output. */
/* ...The above arrays do not distinguish valence and sea contributions, */
/* ...although this information is available internally. The additional */
/* ...commonblock SASVAL provides the valence part only of the above */
/* ...distributions. Array names VXPVMD, VXPANL and VXPANH correspond */
/* ...to XPVMD, XPANL and XPANH, while XPBEH and XPDIR are valence only */
/* ...and therefore not given doubly. VXPDGM gives the sum of valence */
/* ...parts, and so matches XPDFGM. The difference, i.e. XPVMD-VXPVMD */
/* ...and so on, gives the sea part only. */
/*<       SUBROUTINE SASGAM(ISET,X,Q2,P2,IP2,F2GM,XPDFGM) >*/
/* Subroutine */ int sasgam_(integer *iset, real *x, real *q2, real *p2, 
	integer *ip2, real *f2gm, real *xpdfgm)
{
    /* Initialized data */

    static real pmc = (float)1.3;
    static real pmrho = (float).77;
    static real pmphi = (float)1.02;
    static integer nstep = 100;
    static real pmb = (float)4.6;
    static real aem = (float).007297;
    static real aem2pi = (float).0011614;
    static real alam = (float).2;
    static real fracu = (float).8;
    static real frho = (float)2.2;
    static real fomega = (float)23.6;
    static real fphi = (float)18.4;

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);
    double log(doublereal), exp(doublereal), sqrt(doublereal), pow_dd(
	    doublereal *, doublereal *);

    /* Local variables */
    static real q0;
    static integer kf;
    static real q02, q2a;
    static integer kfl;
    static real xpf2, p2mx, facq, facs, xpga[13], chsq, xpbh, p2mxa, p2mxb, 
	    facud, xfval;
    static integer istep;
    static real vxpga[13], q2step;
    extern /* Subroutine */ int sasbeh_(integer *, real *, real *, real *, 
	    real *, real *);
    static real facnor;
    extern /* Subroutine */ int sasano_(integer *, real *, real *, real *, 
	    real *, real *, real *), sasdir_(real *, real *, real *, real *, 
	    real *), sasvmd_(integer *, integer *, real *, real *, real *, 
	    real *, real *, real *);

    /* Fortran I/O blocks */
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };
    static cilist io___16 = { 0, 6, 0, 0, 0 };
    static cilist io___17 = { 0, 6, 0, 0, 0 };


/* ...Purpose: to construct the F2 and parton distributions of the photon */
/* ...by summing homogeneous (VMD) and inhomogeneous (anomalous) terms. */
/* ...For F2, c and b are included by the Bethe-Heitler formula; */
/* ...in the 'MSbar' scheme additionally a Cgamma term is added. */
/*<       DIMENSION XPDFGM(-6:6) >*/
/*<        >*/
/*<       COMMON/SASVAL/VXPVMD(-6:6),VXPANL(-6:6),VXPANH(-6:6),VXPDGM(-6:6) >*/
/*<       SAVE /SASCOM/,/SASVAL/ >*/
/* ...Temporary array. */
/*<       DIMENSION XPGA(-6:6), VXPGA(-6:6)    >*/
/* ...Charm and bottom masses (low to compensate for J/psi etc.). */
/*<       DATA PMC/1.3/, PMB/4.6/ >*/
    /* Parameter adjustments */
    xpdfgm -= -6;

    /* Function Body */
/* ...alpha_em and alpha_em/(2*pi). */
/*<       DATA AEM/0.007297/, AEM2PI/0.0011614/ >*/
/* ...Lambda value for 4 flavours. */
/*<       DATA ALAM/0.20/ >*/
/* ...Mixture u/(u+d), = 0.5 for incoherent and = 0.8 for coherent sum. */
/*<       DATA FRACU/0.8/ >*/
/* ...VMD couplings f_V**2/(4*pi). */
/*<       DATA FRHO/2.20/, FOMEGA/23.6/, FPHI/18.4/ >*/
/* ...Masses for rho (=omega) and phi. */
/*<       DATA PMRHO/0.770/, PMPHI/1.020/ >*/
/* ...Number of points in integration for IP2=1. */
/*<       DATA NSTEP/100/ >*/
/* ...Reset output. */
/*<       F2GM=0. >*/
    *f2gm = (float)0.;
/*<       DO 100 KFL=-6,6 >*/
    for (kfl = -6; kfl <= 6; ++kfl) {
/*<       XPDFGM(KFL)=0. >*/
	xpdfgm[kfl] = (float)0.;
/*<       XPVMD(KFL)=0. >*/
	sascom_1.xpvmd[kfl + 6] = (float)0.;
/*<       XPANL(KFL)=0. >*/
	sascom_1.xpanl[kfl + 6] = (float)0.;
/*<       XPANH(KFL)=0.  >*/
	sascom_1.xpanh[kfl + 6] = (float)0.;
/*<       XPBEH(KFL)=0. >*/
	sascom_1.xpbeh[kfl + 6] = (float)0.;
/*<       XPDIR(KFL)=0. >*/
	sascom_1.xpdir[kfl + 6] = (float)0.;
/*<       VXPVMD(KFL)=0. >*/
	sasval_1.vxpvmd[kfl + 6] = (float)0.;
/*<       VXPANL(KFL)=0. >*/
	sasval_1.vxpanl[kfl + 6] = (float)0.;
/*<       VXPANH(KFL)=0.  >*/
	sasval_1.vxpanh[kfl + 6] = (float)0.;
/*<       VXPDGM(KFL)=0. >*/
	sasval_1.vxpdgm[kfl + 6] = (float)0.;
/*<   100 CONTINUE >*/
/* L100: */
    }
/* ...Check that input sensible. */
/*<       IF(ISET.LE.0.OR.ISET.GE.5) THEN >*/
    if (*iset <= 0 || *iset >= 5) {
/*<         WRITE(*,*) ' FATAL ERROR: SaSgam called for unknown set' >*/
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, " FATAL ERROR: SaSgam called for unknown set", (
		ftnlen)43);
	e_wsle();
/*<         WRITE(*,*) ' ISET = ',ISET >*/
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, " ISET = ", (ftnlen)8);
	do_lio(&c__3, &c__1, (char *)&(*iset), (ftnlen)sizeof(integer));
	e_wsle();
/*<         STOP >*/
	s_stop("", (ftnlen)0);
/*<       ENDIF >*/
    }
/*<       IF(X.LE.0..OR.X.GT.1.) THEN  >*/
    if (*x <= (float)0. || *x > (float)1.) {
/*<         WRITE(*,*) ' FATAL ERROR: SaSgam called for unphysical x' >*/
	s_wsle(&io___16);
	do_lio(&c__9, &c__1, " FATAL ERROR: SaSgam called for unphysical x", (
		ftnlen)44);
	e_wsle();
/*<         WRITE(*,*) ' X = ',X >*/
	s_wsle(&io___17);
	do_lio(&c__9, &c__1, " X = ", (ftnlen)5);
	do_lio(&c__4, &c__1, (char *)&(*x), (ftnlen)sizeof(real));
	e_wsle();
/*<         STOP >*/
	s_stop("", (ftnlen)0);
/*<       ENDIF >*/
    }
/* ...Set Q0 cut-off parameter as function of set used. */
/*<       IF(ISET.LE.2) THEN >*/
    if (*iset <= 2) {
/*<         Q0=0.6 >*/
	q0 = (float).6;
/*<       ELSE >*/
    } else {
/*<         Q0=2. >*/
	q0 = (float)2.;
/*<       ENDIF  >*/
    }
/*<       Q02=Q0**2 >*/
/* Computing 2nd power */
    r__1 = q0;
    q02 = r__1 * r__1;
/* ...Scale choice for off-shell photon; common factors. */
/*<       Q2A=Q2 >*/
    q2a = *q2;
/*<       FACNOR=1. >*/
    facnor = (float)1.;
/*<       IF(IP2.EQ.1) THEN >*/
    if (*ip2 == 1) {
/*<         P2MX=P2+Q02 >*/
	p2mx = *p2 + q02;
/*<         Q2A=Q2+P2*Q02/MAX(Q02,Q2) >*/
	q2a = *q2 + *p2 * q02 / dmax(q02,*q2);
/*<         FACNOR=LOG(Q2/Q02)/NSTEP >*/
	facnor = log(*q2 / q02) / nstep;
/*<       ELSEIF(IP2.EQ.2) THEN >*/
    } else if (*ip2 == 2) {
/*<         P2MX=MAX(P2,Q02) >*/
	p2mx = dmax(*p2,q02);
/*<       ELSEIF(IP2.EQ.3) THEN >*/
    } else if (*ip2 == 3) {
/*<         P2MX=P2+Q02 >*/
	p2mx = *p2 + q02;
/*<         Q2A=Q2+P2*Q02/MAX(Q02,Q2) >*/
	q2a = *q2 + *p2 * q02 / dmax(q02,*q2);
/*<       ELSEIF(IP2.EQ.4) THEN >*/
    } else if (*ip2 == 4) {
/*<        >*/
	p2mx = *q2 * (q02 + *p2) / (*q2 + *p2) * exp(*p2 * (*q2 - q02) / ((*
		q2 + *p2) * (q02 + *p2)));
/*<       ELSEIF(IP2.EQ.5) THEN >*/
    } else if (*ip2 == 5) {
/*<        >*/
	p2mxa = *q2 * (q02 + *p2) / (*q2 + *p2) * exp(*p2 * (*q2 - q02) / ((*
		q2 + *p2) * (q02 + *p2)));
/*<         P2MX=Q0*SQRT(P2MXA) >*/
	p2mx = q0 * sqrt(p2mxa);
/*<         FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MX) >*/
	facnor = log(*q2 / p2mxa) / log(*q2 / p2mx);
/*<       ELSEIF(IP2.EQ.6) THEN >*/
    } else if (*ip2 == 6) {
/*<        >*/
	p2mx = *q2 * (q02 + *p2) / (*q2 + *p2) * exp(*p2 * (*q2 - q02) / ((*
		q2 + *p2) * (q02 + *p2)));
/*<         P2MX=MAX(0.,1.-P2/Q2)*P2MX+MIN(1.,P2/Q2)*MAX(P2,Q02) >*/
/* Computing MAX */
	r__1 = (float)0., r__2 = (float)1. - *p2 / *q2;
/* Computing MIN */
	r__3 = (float)1., r__4 = *p2 / *q2;
	p2mx = dmax(r__1,r__2) * p2mx + dmin(r__3,r__4) * dmax(*p2,q02);
/*<       ELSE >*/
    } else {
/*<        >*/
	p2mxa = *q2 * (q02 + *p2) / (*q2 + *p2) * exp(*p2 * (*q2 - q02) / ((*
		q2 + *p2) * (q02 + *p2)));
/*<         P2MX=Q0*SQRT(P2MXA) >*/
	p2mx = q0 * sqrt(p2mxa);
/*<         P2MXB=P2MX >*/
	p2mxb = p2mx;
/*<         P2MX=MAX(0.,1.-P2/Q2)*P2MX+MIN(1.,P2/Q2)*MAX(P2,Q02) >*/
/* Computing MAX */
	r__1 = (float)0., r__2 = (float)1. - *p2 / *q2;
/* Computing MIN */
	r__3 = (float)1., r__4 = *p2 / *q2;
	p2mx = dmax(r__1,r__2) * p2mx + dmin(r__3,r__4) * dmax(*p2,q02);
/*<         P2MXB=MAX(0.,1.-P2/Q2)*P2MXB+MIN(1.,P2/Q2)*P2MXA >*/
/* Computing MAX */
	r__1 = (float)0., r__2 = (float)1. - *p2 / *q2;
/* Computing MIN */
	r__3 = (float)1., r__4 = *p2 / *q2;
	p2mxb = dmax(r__1,r__2) * p2mxb + dmin(r__3,r__4) * p2mxa;
/*<         FACNOR=LOG(Q2/P2MXA)/LOG(Q2/P2MXB) >*/
	facnor = log(*q2 / p2mxa) / log(*q2 / p2mxb);
/*<       ENDIF >*/
    }
/* ...Call VMD parametrization for d quark and use to give rho, omega, */
/* ...phi. Note dipole dampening for off-shell photon. */
/*<       CALL SASVMD(ISET,1,X,Q2A,P2MX,ALAM,XPGA,VXPGA) >*/
    sasvmd_(iset, &c__1, x, &q2a, &p2mx, &alam, xpga, vxpga);
/*<       XFVAL=VXPGA(1) >*/
    xfval = vxpga[7];
/*<       XPGA(1)=XPGA(2) >*/
    xpga[7] = xpga[8];
/*<       XPGA(-1)=XPGA(-2) >*/
    xpga[5] = xpga[4];
/*<       FACUD=AEM*(1./FRHO+1./FOMEGA)*(PMRHO**2/(PMRHO**2+P2))**2 >*/
/* Computing 2nd power */
    r__2 = pmrho;
/* Computing 2nd power */
    r__3 = pmrho;
/* Computing 2nd power */
    r__1 = r__2 * r__2 / (r__3 * r__3 + *p2);
    facud = aem * ((float)1. / frho + (float)1. / fomega) * (r__1 * r__1);
/*<       FACS=AEM*(1./FPHI)*(PMPHI**2/(PMPHI**2+P2))**2 >*/
/* Computing 2nd power */
    r__2 = pmphi;
/* Computing 2nd power */
    r__3 = pmphi;
/* Computing 2nd power */
    r__1 = r__2 * r__2 / (r__3 * r__3 + *p2);
    facs = aem * ((float)1. / fphi) * (r__1 * r__1);
/*<       DO 110 KFL=-5,5 >*/
    for (kfl = -5; kfl <= 5; ++kfl) {
/*<       XPVMD(KFL)=(FACUD+FACS)*XPGA(KFL) >*/
	sascom_1.xpvmd[kfl + 6] = (facud + facs) * xpga[kfl + 6];
/*<   110 CONTINUE >*/
/* L110: */
    }
/*<       XPVMD(1)=XPVMD(1)+(1.-FRACU)*FACUD*XFVAL >*/
    sascom_1.xpvmd[7] += ((float)1. - fracu) * facud * xfval;
/*<       XPVMD(2)=XPVMD(2)+FRACU*FACUD*XFVAL >*/
    sascom_1.xpvmd[8] += fracu * facud * xfval;
/*<       XPVMD(3)=XPVMD(3)+FACS*XFVAL  >*/
    sascom_1.xpvmd[9] += facs * xfval;
/*<       XPVMD(-1)=XPVMD(-1)+(1.-FRACU)*FACUD*XFVAL >*/
    sascom_1.xpvmd[5] += ((float)1. - fracu) * facud * xfval;
/*<       XPVMD(-2)=XPVMD(-2)+FRACU*FACUD*XFVAL >*/
    sascom_1.xpvmd[4] += fracu * facud * xfval;
/*<       XPVMD(-3)=XPVMD(-3)+FACS*XFVAL  >*/
    sascom_1.xpvmd[3] += facs * xfval;
/*<       VXPVMD(1)=(1.-FRACU)*FACUD*XFVAL >*/
    sasval_1.vxpvmd[7] = ((float)1. - fracu) * facud * xfval;
/*<       VXPVMD(2)=FRACU*FACUD*XFVAL >*/
    sasval_1.vxpvmd[8] = fracu * facud * xfval;
/*<       VXPVMD(3)=FACS*XFVAL  >*/
    sasval_1.vxpvmd[9] = facs * xfval;
/*<       VXPVMD(-1)=(1.-FRACU)*FACUD*XFVAL >*/
    sasval_1.vxpvmd[5] = ((float)1. - fracu) * facud * xfval;
/*<       VXPVMD(-2)=FRACU*FACUD*XFVAL >*/
    sasval_1.vxpvmd[4] = fracu * facud * xfval;
/*<       VXPVMD(-3)=FACS*XFVAL  >*/
    sasval_1.vxpvmd[3] = facs * xfval;
/*<       IF(IP2.NE.1) THEN >*/
    if (*ip2 != 1) {
/* ...Anomalous parametrizations for different strategies */
/* ...for off-shell photons; except full integration. */
/* ...Call anomalous parametrization for d + u + s. */
/*<         CALL SASANO(-3,X,Q2A,P2MX,ALAM,XPGA,VXPGA) >*/
	sasano_(&c_n3, x, &q2a, &p2mx, &alam, xpga, vxpga);
/*<         DO 120 KFL=-5,5 >*/
	for (kfl = -5; kfl <= 5; ++kfl) {
/*<         XPANL(KFL)=FACNOR*XPGA(KFL) >*/
	    sascom_1.xpanl[kfl + 6] = facnor * xpga[kfl + 6];
/*<         VXPANL(KFL)=FACNOR*VXPGA(KFL) >*/
	    sasval_1.vxpanl[kfl + 6] = facnor * vxpga[kfl + 6];
/*<   120   CONTINUE >*/
/* L120: */
	}
/* ...Call anomalous parametrization for c and b. */
/*<         CALL SASANO(4,X,Q2A,P2MX,ALAM,XPGA,VXPGA) >*/
	sasano_(&c__4, x, &q2a, &p2mx, &alam, xpga, vxpga);
/*<         DO 130 KFL=-5,5 >*/
	for (kfl = -5; kfl <= 5; ++kfl) {
/*<         XPANH(KFL)=FACNOR*XPGA(KFL) >*/
	    sascom_1.xpanh[kfl + 6] = facnor * xpga[kfl + 6];
/*<         VXPANH(KFL)=FACNOR*VXPGA(KFL) >*/
	    sasval_1.vxpanh[kfl + 6] = facnor * vxpga[kfl + 6];
/*<   130   CONTINUE	 >*/
/* L130: */
	}
/*<         CALL SASANO(5,X,Q2A,P2MX,ALAM,XPGA,VXPGA) >*/
	sasano_(&c__5, x, &q2a, &p2mx, &alam, xpga, vxpga);
/*<         DO 140 KFL=-5,5 >*/
	for (kfl = -5; kfl <= 5; ++kfl) {
/*<         XPANH(KFL)=XPANH(KFL)+FACNOR*XPGA(KFL) >*/
	    sascom_1.xpanh[kfl + 6] += facnor * xpga[kfl + 6];
/*<         VXPANH(KFL)=VXPANH(KFL)+FACNOR*VXPGA(KFL) >*/
	    sasval_1.vxpanh[kfl + 6] += facnor * vxpga[kfl + 6];
/*<   140   CONTINUE	 >*/
/* L140: */
	}
/*<       ELSE >*/
    } else {
/* ...Special option: loop over flavours and integrate over k2. */
/*<         DO 170 KF=1,5 >*/
	for (kf = 1; kf <= 5; ++kf) {
/*<         DO 160 ISTEP=1,NSTEP >*/
	    i__1 = nstep;
	    for (istep = 1; istep <= i__1; ++istep) {
/*<         Q2STEP=Q02*(Q2/Q02)**((ISTEP-0.5)/NSTEP) >*/
		d__1 = (doublereal) (*q2 / q02);
		d__2 = (doublereal) ((istep - (float).5) / nstep);
		q2step = q02 * pow_dd(&d__1, &d__2);
/*<        >*/
/* Computing 2nd power */
		r__1 = pmc;
/* Computing 2nd power */
		r__2 = pmb;
		if (kf == 4 && q2step < r__1 * r__1 || kf == 5 && q2step < 
			r__2 * r__2) {
		    goto L160;
		}
/*<         CALL SASVMD(0,KF,X,Q2,Q2STEP,ALAM,XPGA,VXPGA) >*/
		sasvmd_(&c__0, &kf, x, q2, &q2step, &alam, xpga, vxpga);
/*<         FACQ=AEM2PI*(Q2STEP/(Q2STEP+P2))**2*FACNOR >*/
/* Computing 2nd power */
		r__1 = q2step / (q2step + *p2);
		facq = aem2pi * (r__1 * r__1) * facnor;
/*<         IF(MOD(KF,2).EQ.0) FACQ=FACQ*(8./9.) >*/
		if (kf % 2 == 0) {
		    facq *= (float).88888888888888884;
		}
/*<         IF(MOD(KF,2).EQ.1) FACQ=FACQ*(2./9.) >*/
		if (kf % 2 == 1) {
		    facq *= (float).22222222222222221;
		}
/*<         DO 150 KFL=-5,5 >*/
		for (kfl = -5; kfl <= 5; ++kfl) {
/*<         IF(KF.LE.3) XPANL(KFL)=XPANL(KFL)+FACQ*XPGA(KFL)  >*/
		    if (kf <= 3) {
			sascom_1.xpanl[kfl + 6] += facq * xpga[kfl + 6];
		    }
/*<         IF(KF.GE.4) XPANH(KFL)=XPANH(KFL)+FACQ*XPGA(KFL)  >*/
		    if (kf >= 4) {
			sascom_1.xpanh[kfl + 6] += facq * xpga[kfl + 6];
		    }
/*<         IF(KF.LE.3) VXPANL(KFL)=VXPANL(KFL)+FACQ*VXPGA(KFL)  >*/
		    if (kf <= 3) {
			sasval_1.vxpanl[kfl + 6] += facq * vxpga[kfl + 6];
		    }
/*<         IF(KF.GE.4) VXPANH(KFL)=VXPANH(KFL)+FACQ*VXPGA(KFL)  >*/
		    if (kf >= 4) {
			sasval_1.vxpanh[kfl + 6] += facq * vxpga[kfl + 6];
		    }
/*<   150   CONTINUE >*/
/* L150: */
		}
/*<   160   CONTINUE >*/
L160:
		;
	    }
/*<   170   CONTINUE >*/
/* L170: */
	}
/*<       ENDIF >*/
    }
/* ...Call Bethe-Heitler term expression for charm and bottom. */
/*<       CALL SASBEH(4,X,Q2,P2,PMC**2,XPBH) >*/
/* Computing 2nd power */
    r__2 = pmc;
    r__1 = r__2 * r__2;
    sasbeh_(&c__4, x, q2, p2, &r__1, &xpbh);
/*<       XPBEH(4)=XPBH >*/
    sascom_1.xpbeh[10] = xpbh;
/*<       XPBEH(-4)=XPBH >*/
    sascom_1.xpbeh[2] = xpbh;
/*<       CALL SASBEH(5,X,Q2,P2,PMB**2,XPBH) >*/
/* Computing 2nd power */
    r__2 = pmb;
    r__1 = r__2 * r__2;
    sasbeh_(&c__5, x, q2, p2, &r__1, &xpbh);
/*<       XPBEH(5)=XPBH >*/
    sascom_1.xpbeh[11] = xpbh;
/*<       XPBEH(-5)=XPBH >*/
    sascom_1.xpbeh[1] = xpbh;
/* ...For MSbar subtraction call C^gamma term expression for d, u, s. */
/*<       IF(ISET.EQ.2.OR.ISET.EQ.4) THEN >*/
    if (*iset == 2 || *iset == 4) {
/*<         CALL SASDIR(X,Q2,P2,Q02,XPGA) >*/
	sasdir_(x, q2, p2, &q02, xpga);
/*<         DO 180 KFL=-5,5 >*/
	for (kfl = -5; kfl <= 5; ++kfl) {
/*<         XPDIR(KFL)=XPGA(KFL) >*/
	    sascom_1.xpdir[kfl + 6] = xpga[kfl + 6];
/*<   180   CONTINUE >*/
/* L180: */
	}
/*<       ENDIF >*/
    }
/* ...Store result in output array. */
/*<       DO 190 KFL=-5,5 >*/
    for (kfl = -5; kfl <= 5; ++kfl) {
/*<       CHSQ=1./9. >*/
	chsq = (float).1111111111111111;
/*<       IF(IABS(KFL).EQ.2.OR.IABS(KFL).EQ.4) CHSQ=4./9. >*/
	if (abs(kfl) == 2 || abs(kfl) == 4) {
	    chsq = (float).44444444444444442;
	}
/*<       XPF2=XPVMD(KFL)+XPANL(KFL)+XPBEH(KFL)+XPDIR(KFL) >*/
	xpf2 = sascom_1.xpvmd[kfl + 6] + sascom_1.xpanl[kfl + 6] + 
		sascom_1.xpbeh[kfl + 6] + sascom_1.xpdir[kfl + 6];
/*<       IF(KFL.NE.0) F2GM=F2GM+CHSQ*XPF2     >*/
	if (kfl != 0) {
	    *f2gm += chsq * xpf2;
	}
/*<       XPDFGM(KFL)=XPVMD(KFL)+XPANL(KFL)+XPANH(KFL) >*/
	xpdfgm[kfl] = sascom_1.xpvmd[kfl + 6] + sascom_1.xpanl[kfl + 6] + 
		sascom_1.xpanh[kfl + 6];
/*<       VXPDGM(KFL)=VXPVMD(KFL)+VXPANL(KFL)+VXPANH(KFL) >*/
	sasval_1.vxpdgm[kfl + 6] = sasval_1.vxpvmd[kfl + 6] + sasval_1.vxpanl[
		kfl + 6] + sasval_1.vxpanh[kfl + 6];
/*<   190 CONTINUE      >*/
/* L190: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* sasgam_ */

/* ********************************************************************* */
/*<       SUBROUTINE SASVMD(ISET,KF,X,Q2,P2,ALAM,XPGA,VXPGA) >*/
/* Subroutine */ int sasvmd_(integer *iset, integer *kf, real *x, real *q2, 
	real *p2, real *alam, real *xpga, real *vxpga)
{
    /* Initialized data */

    static real pmc = (float)1.3;
    static real pmb = (float)4.6;

    /* System generated locals */
    real r__1, r__2, r__3, r__4, r__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal), exp(
	    doublereal);

    /* Local variables */
    static real s, s2, s3, s4, x1, xl;
    static integer kfa, kfl;
    static real sch;
    static integer nfp, nfq;
    static real sbt, sll, xsea, xchm, xval, xbot, xglu, alam3, alam5, p2eff, 
	    q2eff, xsea0, p2div, q2div;

/* ...Purpose: to evaluate the VMD parton distributions of a photon, */
/* ...evolved homogeneously from an initial scale P2 to Q2. */
/* ...Does not include dipole suppression factor. */
/* ...ISET is parton distribution set, see above; */
/* ...additionally ISET=0 is used for the evolution of an anomalous photon */
/* ...which branched at a scale P2 and then evolved homogeneously to Q2. */
/* ...ALAM is the 4-flavour Lambda, which is automatically converted */
/* ...to 3- and 5-flavour equivalents as needed. */
/*<       DIMENSION XPGA(-6:6), VXPGA(-6:6) >*/
/*<       DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/ >*/
    /* Parameter adjustments */
    vxpga -= -6;
    xpga -= -6;

    /* Function Body */
/* ...Reset output. */
/*<       DO 100 KFL=-6,6 >*/
    for (kfl = -6; kfl <= 6; ++kfl) {
/*<       XPGA(KFL)=0. >*/
	xpga[kfl] = (float)0.;
/*<       VXPGA(KFL)=0. >*/
	vxpga[kfl] = (float)0.;
/*<   100 CONTINUE >*/
/* L100: */
    }
/*<       KFA=IABS(KF) >*/
    kfa = abs(*kf);
/* ...Calculate Lambda; protect against unphysical Q2 and P2 input. */
/*<       ALAM3=ALAM*(PMC/ALAM)**(2./27.) >*/
    d__1 = (doublereal) (pmc / *alam);
    alam3 = *alam * pow_dd(&d__1, &c_b39);
/*<       ALAM5=ALAM*(ALAM/PMB)**(2./23.) >*/
    d__1 = (doublereal) (*alam / pmb);
    alam5 = *alam * pow_dd(&d__1, &c_b40);
/*<       P2EFF=MAX(P2,1.2*ALAM3**2) >*/
/* Computing MAX */
/* Computing 2nd power */
    r__3 = alam3;
    r__1 = *p2, r__2 = r__3 * r__3 * (float)1.2;
    p2eff = dmax(r__1,r__2);
/*<       IF(KFA.EQ.4) P2EFF=MAX(P2EFF,PMC**2) >*/
    if (kfa == 4) {
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmc;
	r__1 = p2eff, r__2 = r__3 * r__3;
	p2eff = dmax(r__1,r__2);
    }
/*<       IF(KFA.EQ.5) P2EFF=MAX(P2EFF,PMB**2) >*/
    if (kfa == 5) {
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmb;
	r__1 = p2eff, r__2 = r__3 * r__3;
	p2eff = dmax(r__1,r__2);
    }
/*<       Q2EFF=MAX(Q2,P2EFF) >*/
    q2eff = dmax(*q2,p2eff);
/* ...Find number of flavours at lower and upper scale. */
/*<       NFP=4 >*/
    nfp = 4;
/*<       IF(P2EFF.LT.PMC**2) NFP=3 >*/
/* Computing 2nd power */
    r__1 = pmc;
    if (p2eff < r__1 * r__1) {
	nfp = 3;
    }
/*<       IF(P2EFF.GT.PMB**2) NFP=5 >*/
/* Computing 2nd power */
    r__1 = pmb;
    if (p2eff > r__1 * r__1) {
	nfp = 5;
    }
/*<       NFQ=4 >*/
    nfq = 4;
/*<       IF(Q2EFF.LT.PMC**2) NFQ=3 >*/
/* Computing 2nd power */
    r__1 = pmc;
    if (q2eff < r__1 * r__1) {
	nfq = 3;
    }
/*<       IF(Q2EFF.GT.PMB**2) NFQ=5 >*/
/* Computing 2nd power */
    r__1 = pmb;
    if (q2eff > r__1 * r__1) {
	nfq = 5;
    }
/* ...Find s as sum of 3-, 4- and 5-flavour parts. */
/*<       S=0. >*/
    s = (float)0.;
/*<       IF(NFP.EQ.3) THEN >*/
    if (nfp == 3) {
/*<         Q2DIV=PMC**2 >*/
/* Computing 2nd power */
	r__1 = pmc;
	q2div = r__1 * r__1;
/*<         IF(NFQ.EQ.3) Q2DIV=Q2EFF >*/
	if (nfq == 3) {
	    q2div = q2eff;
	}
/*<         S=S+(6./27.)*LOG(LOG(Q2DIV/ALAM3**2)/LOG(P2EFF/ALAM3**2)) >*/
/* Computing 2nd power */
	r__1 = alam3;
/* Computing 2nd power */
	r__2 = alam3;
	s += log(log(q2div / (r__1 * r__1)) / log(p2eff / (r__2 * r__2))) * (
		float).22222222222222221;
/*<       ENDIF >*/
    }
/*<       IF(NFP.LE.4.AND.NFQ.GE.4) THEN >*/
    if (nfp <= 4 && nfq >= 4) {
/*<         P2DIV=P2EFF >*/
	p2div = p2eff;
/*<         IF(NFP.EQ.3) P2DIV=PMC**2 >*/
	if (nfp == 3) {
/* Computing 2nd power */
	    r__1 = pmc;
	    p2div = r__1 * r__1;
	}
/*<         Q2DIV=Q2EFF >*/
	q2div = q2eff;
/*<         IF(NFQ.EQ.5) Q2DIV=PMB**2  >*/
	if (nfq == 5) {
/* Computing 2nd power */
	    r__1 = pmb;
	    q2div = r__1 * r__1;
	}
/*<         S=S+(6./25.)*LOG(LOG(Q2DIV/ALAM**2)/LOG(P2DIV/ALAM**2)) >*/
/* Computing 2nd power */
	r__1 = *alam;
/* Computing 2nd power */
	r__2 = *alam;
	s += log(log(q2div / (r__1 * r__1)) / log(p2div / (r__2 * r__2))) * (
		float).23999999999999999;
/*<       ENDIF >*/
    }
/*<       IF(NFQ.EQ.5) THEN >*/
    if (nfq == 5) {
/*<         P2DIV=PMB**2 >*/
/* Computing 2nd power */
	r__1 = pmb;
	p2div = r__1 * r__1;
/*<         IF(NFP.EQ.5) P2DIV=P2EFF >*/
	if (nfp == 5) {
	    p2div = p2eff;
	}
/*<         S=S+(6./23.)*LOG(LOG(Q2EFF/ALAM5**2)/LOG(P2DIV/ALAM5**2)) >*/
/* Computing 2nd power */
	r__1 = alam5;
/* Computing 2nd power */
	r__2 = alam5;
	s += log(log(q2eff / (r__1 * r__1)) / log(p2div / (r__2 * r__2))) * (
		float).2608695652173913;
/*<       ENDIF >*/
    }
/* ...Calculate frequent combinations of x and s. */
/*<       X1=1.-X >*/
    x1 = (float)1. - *x;
/*<       XL=-LOG(X) >*/
    xl = -log(*x);
/*<       S2=S**2 >*/
/* Computing 2nd power */
    r__1 = s;
    s2 = r__1 * r__1;
/*<       S3=S**3 >*/
/* Computing 3rd power */
    r__1 = s;
    s3 = r__1 * (r__1 * r__1);
/*<       S4=S**4 >*/
/* Computing 4th power */
    r__1 = s, r__1 *= r__1;
    s4 = r__1 * r__1;
/* ...Evaluate homogeneous anomalous parton distributions below or */
/* ...above threshold. */
/*<       IF(ISET.EQ.0) THEN >*/
    if (*iset == 0) {
/*<        >*/
/* Computing 2nd power */
	r__1 = pmc;
/* Computing 2nd power */
	r__2 = pmb;
	if (*q2 <= *p2 || kfa == 4 && *q2 < r__1 * r__1 || kfa == 5 && *q2 < 
		r__2 * r__2) {
/*<         XVAL = X * 1.5 * (X**2+X1**2) >*/
/* Computing 2nd power */
	    r__1 = *x;
/* Computing 2nd power */
	    r__2 = x1;
	    xval = *x * (float)1.5 * (r__1 * r__1 + r__2 * r__2);
/*<         XGLU = 0. >*/
	    xglu = (float)0.;
/*<         XSEA = 0. >*/
	    xsea = (float)0.;
/*<       ELSE >*/
	} else {
/*<        >*/
/* Computing 2nd power */
	    r__1 = *x;
/* Computing 2nd power */
	    r__2 = x1;
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((float)1. / (s * (float)1.5 + (float)1.));
/* Computing 2nd power */
	    r__3 = *x;
	    d__3 = (doublereal) ((float)1. - r__3 * r__3);
	    d__4 = (doublereal) (s * (float)2.667);
	    xval = ((float)1.5 / ((float)1. - s * (float).197 + s2 * (float)
		    4.33) * (r__1 * r__1) + (s * (float)2.1 + (float)1.5) / (
		    s * (float)3.29 + (float)1.) * (r__2 * r__2) + s * (float)
		    5.23 / (s * (float)1.17 + (float)1. + s3 * (float)19.9) * 
		    *x * x1) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float)-2.03 / (s * (float)2.44 + (float)
		    1.));
	    d__3 = (doublereal) (x1 * xl);
	    d__4 = (doublereal) (s * (float)1.333);
/* Computing 2nd power */
	    r__1 = *x;
	    xglu = s * (float)4. / (s * (float)4.76 + (float)1. + s2 * (float)
		    15.2 + s4 * (float)29.3) * pow_dd(&d__1, &d__2) * pow_dd(&
		    d__3, &d__4) * ((r__1 * r__1 * (float)4. + *x * (float)7. 
		    + (float)4.) * x1 / (float)3. - *x * (float)2. * (*x + (
		    float)1.) * xl);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float)-1.54 / (s * (float)1.29 + (float)
		    1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (s * (float)2.667);
/* Computing 2nd power */
	    r__1 = *x;
/* Computing 2nd power */
	    r__2 = *x;
/* Computing 2nd power */
	    r__3 = xl;
	    xsea = s2 / (s * (float)4.54 + (float)1. + s2 * (float)8.19 + s3 *
		     (float)8.05) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &
		    d__4) * (((float)8. - *x * (float)73. + r__1 * r__1 * (
		    float)62.) * x1 / (float)9. + ((float)3. - r__2 * r__2 * (
		    float)8. / (float)3.) * *x * xl + (*x * (float)2. - (
		    float)1.) * *x * (r__3 * r__3));
/*<       ENDIF >*/
	}
/* ...Evaluate set 1D parton distributions below or above threshold. */
/*<       ELSEIF(ISET.EQ.1) THEN >*/
    } else if (*iset == 1) {
/*<        >*/
/* Computing 2nd power */
	r__1 = pmc;
/* Computing 2nd power */
	r__2 = pmb;
	if (*q2 <= *p2 || kfa == 4 && *q2 < r__1 * r__1 || kfa == 5 && *q2 < 
		r__2 * r__2) {
/*<         XVAL = 1.294 * X**0.80 * X1**0.76 >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xval = pow_dd(&d__1, &c_b41) * (float)1.294 * pow_dd(&d__2, &
		    c_b42);
/*<         XGLU = 1.273 * X**0.40 * X1**1.76 >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xglu = pow_dd(&d__1, &c_b43) * (float)1.273 * pow_dd(&d__2, &
		    c_b44);
/*<         XSEA = 0.100 * X1**3.76 >*/
	    d__1 = (doublereal) x1;
	    xsea = pow_dd(&d__1, &c_b45) * (float).1;
/*<       ELSE >*/
	} else {
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((float).8 - s * (float).13);
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (s * (float).667 + (float).76);
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)2.);
	    xval = (float)1.294 / (s * (float).252 + (float)1. + s2 * (float)
		    3.079) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4) * 
		    pow_dd(&d__5, &d__6);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float)-1.9 / (s * (float)3.6 + (float)
		    1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) xl;
	    d__5 = (doublereal) (s * (float)3. + (float).5);
	    d__6 = (doublereal) (*x);
	    d__7 = (doublereal) x1;
	    d__8 = (doublereal) (s * (float)3. + (float)1.76);
	    xglu = s * (float)7.9 / (s * (float)5.5 + (float)1.) * exp(s * (
		    float)-5.16) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &
		    c_b46) * pow_dd(&d__4, &d__5) + exp(s * (float)-10.) * (
		    float)1.273 * pow_dd(&d__6, &c_b43) * pow_dd(&d__7, &d__8)
		    ;
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s2 * (float)-7.32 / (s2 * (float)10.3 + (
		    float)1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) ((s * (float)15. + (float)3.76 + s2 * (float)
		    12.) / (s * (float)4. + (float)1.));
	    xsea = ((float).1 - s2 * (float).397 + s3 * (float)1.121) / (s2 * 
		    (float)5.61 + (float)1. + s3 * (float)5.26) * pow_dd(&
		    d__1, &d__2) * pow_dd(&d__3, &d__4);
/*<         XSEA0 = 0.100 * X1**3.76 >*/
	    d__1 = (doublereal) x1;
	    xsea0 = pow_dd(&d__1, &c_b45) * (float).1;
/*<       ENDIF >*/
	}
/* ...Evaluate set 1M parton distributions below or above threshold. */
/*<       ELSEIF(ISET.EQ.2) THEN >*/
    } else if (*iset == 2) {
/*<        >*/
/* Computing 2nd power */
	r__1 = pmc;
/* Computing 2nd power */
	r__2 = pmb;
	if (*q2 <= *p2 || kfa == 4 && *q2 < r__1 * r__1 || kfa == 5 && *q2 < 
		r__2 * r__2) {
/*<         XVAL = 0.8477 * X**0.51 * X1**1.37 >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xval = pow_dd(&d__1, &c_b49) * (float).8477 * pow_dd(&d__2, &
		    c_b50);
/*<         XGLU = 3.42 * X**0.255 * X1**2.37 >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xglu = pow_dd(&d__1, &c_b51) * (float)3.42 * pow_dd(&d__2, &c_b52)
		    ;
/*<         XSEA = 0. >*/
	    xsea = (float)0.;
/*<       ELSE >*/
	} else {
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float).21 + (float).51);
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) xl;
	    d__5 = (doublereal) (s * (float)2.667);
	    xval = (float).8477 / (s * (float)1.37 + (float)1. + s2 * (float)
		    2.18 + s3 * (float)3.73) * pow_dd(&d__1, &d__2) * pow_dd(&
		    d__3, &c_b50) * pow_dd(&d__4, &d__5);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (((float)-.013 - s * (float)1.8) / (s * (
		    float)3.14 + (float)1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (s * (float).4 + (float)2.37);
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)3.6 + (float).32);
	    d__7 = (doublereal) (*x);
	    d__8 = (doublereal) x1;
	    d__9 = (doublereal) (s * (float)3. + (float)2.37);
	    xglu = s * (float)24. / (s * (float)9.6 + (float)1. + s2 * (float)
		    .92 + s3 * (float)14.34) * exp(s * (float)-5.94) * pow_dd(
		    &d__1, &d__2) * pow_dd(&d__3, &d__4) * pow_dd(&d__5, &
		    d__6) + exp(s * (float)-12.) * (float)3.42 * pow_dd(&d__7,
		     &c_b51) * pow_dd(&d__8, &d__9);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (((float).13 - s * (float)2.9) / (s * (float)
		    5.44 + (float)1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (s * (float).5 + (float)3.45);
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)2.8);
	    xsea = s * (float).842 / (s * (float)21.3 + (float)1. - s2 * (
		    float)33.2 + s3 * (float)229.) * pow_dd(&d__1, &d__2) * 
		    pow_dd(&d__3, &d__4) * pow_dd(&d__5, &d__6);
/*<         XSEA0 = 0. >*/
	    xsea0 = (float)0.;
/*<       ENDIF >*/
	}
/* ...Evaluate set 2D parton distributions below or above threshold. */
/*<       ELSEIF(ISET.EQ.3) THEN >*/
    } else if (*iset == 3) {
/*<        >*/
/* Computing 2nd power */
	r__1 = pmc;
/* Computing 2nd power */
	r__2 = pmb;
	if (*q2 <= *p2 || kfa == 4 && *q2 < r__1 * r__1 || kfa == 5 && *q2 < 
		r__2 * r__2) {
/*<         XVAL = X**0.46 * X1**0.64 + 0.76 * X >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xval = pow_dd(&d__1, &c_b55) * pow_dd(&d__2, &c_b56) + *x * (
		    float).76;
/*<         XGLU = 1.925 * X1**2 >*/
/* Computing 2nd power */
	    r__1 = x1;
	    xglu = r__1 * r__1 * (float)1.925;
/*<         XSEA = 0.242 * X1**4 >*/
/* Computing 4th power */
	    r__1 = x1, r__1 *= r__1;
	    xsea = r__1 * r__1 * (float).242;
/*<       ELSE >*/
	} else {
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float).25 + (float).46);
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) ((s * (float).14 + (float).64 + s2 * (float)
		    5.) / (s + (float)1.));
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)1.9);
	    d__7 = (doublereal) x1;
	    d__8 = (doublereal) (s * (float)2.667);
	    xval = (s * (float).186 + (float)1.) / ((float)1. - s * (float)
		    .209 + s2 * (float)1.495) * pow_dd(&d__1, &d__2) * pow_dd(
		    &d__3, &d__4) * pow_dd(&d__5, &d__6) + (s * (float).4 + (
		    float).76) * *x * pow_dd(&d__7, &d__8);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((s * (float)-5.81 - s2 * (float)5.34) / (s * 
		    (float)29. + (float)1. - s2 * (float)4.26));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (((float)2. - s * (float)5.9) / (s * (float)
		    1.7 + (float)1.));
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)9.3 / (s * (float)1.7 + (float)1.)
		    );
	    xglu = (s * (float)5.55 + (float)1.925 + s2 * (float)147.) / ((
		    float)1. - s * (float)3.59 + s2 * (float)3.32) * exp(s * (
		    float)-18.67) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &
		    d__4) * pow_dd(&d__5, &d__6);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s2 * (float)-12.1 / (s * (float)2.62 + (
		    float)1. + s2 * (float)16.7));
/* Computing 4th power */
	    r__1 = x1, r__1 *= r__1;
	    d__3 = (doublereal) xl;
	    d__4 = (doublereal) s;
	    xsea = ((float).242 - s * (float).252 + s2 * (float)1.19) / ((
		    float)1. - s * (float).607 + s2 * (float)21.95) * pow_dd(&
		    d__1, &d__2) * (r__1 * r__1) * pow_dd(&d__3, &d__4);
/*<         XSEA0 = 0.242 * X1**4 >*/
/* Computing 4th power */
	    r__1 = x1, r__1 *= r__1;
	    xsea0 = r__1 * r__1 * (float).242;
/*<       ENDIF  >*/
	}
/* ...Evaluate set 2M parton distributions below or above threshold. */
/*<       ELSEIF(ISET.EQ.4) THEN >*/
    } else if (*iset == 4) {
/*<        >*/
/* Computing 2nd power */
	r__1 = pmc;
/* Computing 2nd power */
	r__2 = pmb;
	if (*q2 <= *p2 || kfa == 4 && *q2 < r__1 * r__1 || kfa == 5 && *q2 < 
		r__2 * r__2) {
/*<         XVAL = 1.168 * X**0.50 * X1**2.60 + 0.965 * X >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) x1;
	    xval = pow_dd(&d__1, &c_b57) * (float)1.168 * pow_dd(&d__2, &
		    c_b58) + *x * (float).965;
/*<         XGLU = 1.808 * X1**2 >*/
/* Computing 2nd power */
	    r__1 = x1;
	    xglu = r__1 * r__1 * (float)1.808;
/*<         XSEA = 0.209 * X1**4   >*/
/* Computing 4th power */
	    r__1 = x1, r__1 *= r__1;
	    xsea = r__1 * r__1 * (float).209;
/*<       ELSE >*/
	} else {
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((s * (float).208 + (float).5) / ((float)1. - 
		    s * (float).794 + s2 * (float)1.516));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) ((s * (float)7.6 + (float)2.6) / (s * (float)
		    5. + (float)1.));
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)5.15 / (s * (float)2. + (float)1.)
		    );
	    d__7 = (doublereal) x1;
	    d__8 = (doublereal) (s * (float)2.667);
	    xval = (s * (float)1.771 + (float)1.168 + s2 * (float)29.35) * 
		    exp(s * (float)-5.776) * pow_dd(&d__1, &d__2) * pow_dd(&
		    d__3, &d__4) * pow_dd(&d__5, &d__6) + (s * (float)22.35 + 
		    (float).965) / (s * (float)18.4 + (float)1.) * *x * 
		    pow_dd(&d__7, &d__8);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((s * (float)-5.35 - s2 * (float)10.11) / (s *
		     (float)31.71 + (float)1.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (((float)2. - s * (float)7.3 + s2 * (float)4.)
		     / (s * (float)2.5 + (float)1.));
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float)10.9 / (s * (float)2.5 + (float)
		    1.));
	    xglu = (s * (float)29.9 + (float)1.808) / (s * (float)26.4 + (
		    float)1.) * exp(s * (float)-5.28) * pow_dd(&d__1, &d__2) *
		     pow_dd(&d__3, &d__4) * pow_dd(&d__5, &d__6);
/*<        >*/
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((s * (float)-.373 - s2 * (float)7.71) / (s * 
		    (float).815 + (float)1. + s2 * (float)11.));
	    d__3 = (doublereal) x1;
	    d__4 = (doublereal) (s + (float)4.);
	    d__5 = (doublereal) xl;
	    d__6 = (doublereal) (s * (float).45);
	    xsea = (s2 * (float).644 + (float).209) / (s * (float).319 + (
		    float)1. + s2 * (float)17.6) * pow_dd(&d__1, &d__2) * 
		    pow_dd(&d__3, &d__4) * pow_dd(&d__5, &d__6);
/*<         XSEA0 = 0.209 * X1**4   >*/
/* Computing 4th power */
	    r__1 = x1, r__1 *= r__1;
	    xsea0 = r__1 * r__1 * (float).209;
/*<       ENDIF >*/
	}
/*<       ENDIF >*/
    }
/* ...Threshold factors for c and b sea. */
/*<       SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2)) >*/
/* Computing 2nd power */
    r__1 = *alam;
/* Computing 2nd power */
    r__2 = *alam;
    sll = log(log(q2eff / (r__1 * r__1)) / log(p2eff / (r__2 * r__2)));
/*<       XCHM=0.     >*/
    xchm = (float)0.;
/*<       IF(Q2.GT.PMC**2.AND.Q2.GT.1.001*P2EFF) THEN >*/
/* Computing 2nd power */
    r__1 = pmc;
    if (*q2 > r__1 * r__1 && *q2 > p2eff * (float)1.001) {
/*<         SCH=MAX(0.,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))   >*/
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmc;
/* Computing 2nd power */
	r__4 = *alam;
/* Computing 2nd power */
	r__5 = *alam;
	r__1 = (float)0., r__2 = log(log(r__3 * r__3 / (r__4 * r__4)) / log(
		p2eff / (r__5 * r__5)));
	sch = dmax(r__1,r__2);
/*<         IF(ISET.EQ.0) THEN >*/
	if (*iset == 0) {
/*<           XCHM=XSEA*(1.-(SCH/SLL)**2) >*/
/* Computing 2nd power */
	    r__1 = sch / sll;
	    xchm = xsea * ((float)1. - r__1 * r__1);
/*<         ELSE >*/
	} else {
/*<           XCHM=MAX(0.,XSEA-XSEA0*X1**(2.667*S))*(1.-SCH/SLL) >*/
/* Computing MAX */
	    d__1 = (doublereal) x1;
	    d__2 = (doublereal) (s * (float)2.667);
	    r__1 = (float)0., r__2 = xsea - xsea0 * pow_dd(&d__1, &d__2);
	    xchm = dmax(r__1,r__2) * ((float)1. - sch / sll);
/*<         ENDIF >*/
	}
/*<       ENDIF      >*/
    }
/*<       XBOT=0. >*/
    xbot = (float)0.;
/*<       IF(Q2.GT.PMB**2.AND.Q2.GT.1.001*P2EFF) THEN >*/
/* Computing 2nd power */
    r__1 = pmb;
    if (*q2 > r__1 * r__1 && *q2 > p2eff * (float)1.001) {
/*<         SBT=MAX(0.,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))   >*/
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmb;
/* Computing 2nd power */
	r__4 = *alam;
/* Computing 2nd power */
	r__5 = *alam;
	r__1 = (float)0., r__2 = log(log(r__3 * r__3 / (r__4 * r__4)) / log(
		p2eff / (r__5 * r__5)));
	sbt = dmax(r__1,r__2);
/*<         IF(ISET.EQ.0) THEN >*/
	if (*iset == 0) {
/*<           XBOT=XSEA*(1.-(SBT/SLL)**2) >*/
/* Computing 2nd power */
	    r__1 = sbt / sll;
	    xbot = xsea * ((float)1. - r__1 * r__1);
/*<         ELSE >*/
	} else {
/*<           XBOT=MAX(0.,XSEA-XSEA0*X1**(2.667*S))*(1.-SBT/SLL) >*/
/* Computing MAX */
	    d__1 = (doublereal) x1;
	    d__2 = (doublereal) (s * (float)2.667);
	    r__1 = (float)0., r__2 = xsea - xsea0 * pow_dd(&d__1, &d__2);
	    xbot = dmax(r__1,r__2) * ((float)1. - sbt / sll);
/*<         ENDIF   >*/
	}
/*<       ENDIF    >*/
    }
/* ...Fill parton distributions. */
/*<       XPGA(0)=XGLU >*/
    xpga[0] = xglu;
/*<       XPGA(1)=XSEA >*/
    xpga[1] = xsea;
/*<       XPGA(2)=XSEA >*/
    xpga[2] = xsea;
/*<       XPGA(3)=XSEA >*/
    xpga[3] = xsea;
/*<       XPGA(4)=XCHM >*/
    xpga[4] = xchm;
/*<       XPGA(5)=XBOT >*/
    xpga[5] = xbot;
/*<       XPGA(KFA)=XPGA(KFA)+XVAL >*/
    xpga[kfa] += xval;
/*<       DO 110 KFL=1,5 >*/
    for (kfl = 1; kfl <= 5; ++kfl) {
/*<       XPGA(-KFL)=XPGA(KFL) >*/
	xpga[-kfl] = xpga[kfl];
/*<   110 CONTINUE >*/
/* L110: */
    }
/*<       VXPGA(KFA)=XVAL >*/
    vxpga[kfa] = xval;
/*<       VXPGA(-KFA)=XVAL >*/
    vxpga[-kfa] = xval;
/*<       RETURN >*/
    return 0;
/*<       END  >*/
} /* sasvmd_ */

/* ********************************************************************* */
/*<       SUBROUTINE SASANO(KF,X,Q2,P2,ALAM,XPGA,VXPGA) >*/
/* Subroutine */ int sasano_(integer *kf, real *x, real *q2, real *p2, real *
	alam, real *xpga, real *vxpga)
{
    /* Initialized data */

    static real pmc = (float)1.3;
    static real pmb = (float)4.6;
    static real aem2pi = (float).0011614;

    /* System generated locals */
    integer i__1;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    static real s, xl, fac;
    static integer kfa, kfl;
    static real sch;
    static integer nfp, nfq;
    static real sbt, sll, snf3, snf4, chsq, xsea, xchm, snfp, snfq, xval, 
	    xbot, xglu, p2eff, q2eff, q2div, tdiff;
    static integer kflmn, kflmx;
    static real alamsq[3];

/* ...Purpose: to evaluate the parton distributions of the anomalous */
/* ...photon, inhomogeneously evolved from a scale P2 (where it vanishes) */
/* ...to Q2. */
/* ...KF=0 gives the sum over (up to) 5 flavours, */
/* ...KF<0 limits to flavours up to abs(KF), */
/* ...KF>0 is for flavour KF only. */
/* ...ALAM is the 4-flavour Lambda, which is automatically converted */
/* ...to 3- and 5-flavour equivalents as needed. */
/*<       DIMENSION XPGA(-6:6), VXPGA(-6:6), ALAMSQ(3:5) >*/
/*<       DATA PMC/1.3/, PMB/4.6/, AEM/0.007297/, AEM2PI/0.0011614/ >*/
    /* Parameter adjustments */
    vxpga -= -6;
    xpga -= -6;

    /* Function Body */
/* ...Reset output. */
/*<       DO 100 KFL=-6,6 >*/
    for (kfl = -6; kfl <= 6; ++kfl) {
/*<       XPGA(KFL)=0. >*/
	xpga[kfl] = (float)0.;
/*<       VXPGA(KFL)=0. >*/
	vxpga[kfl] = (float)0.;
/*<   100 CONTINUE >*/
/* L100: */
    }
/*<       IF(Q2.LE.P2) RETURN >*/
    if (*q2 <= *p2) {
	return 0;
    }
/*<       KFA=IABS(KF) >*/
    kfa = abs(*kf);
/* ...Calculate Lambda; protect against unphysical Q2 and P2 input. */
/*<       ALAMSQ(3)=(ALAM*(PMC/ALAM)**(2./27.))**2 >*/
    d__1 = (doublereal) (pmc / *alam);
/* Computing 2nd power */
    r__1 = *alam * pow_dd(&d__1, &c_b39);
    alamsq[0] = r__1 * r__1;
/*<       ALAMSQ(4)=ALAM**2 >*/
/* Computing 2nd power */
    r__1 = *alam;
    alamsq[1] = r__1 * r__1;
/*<       ALAMSQ(5)=(ALAM*(ALAM/PMB)**(2./23.))**2 >*/
    d__1 = (doublereal) (*alam / pmb);
/* Computing 2nd power */
    r__1 = *alam * pow_dd(&d__1, &c_b40);
    alamsq[2] = r__1 * r__1;
/*<       P2EFF=MAX(P2,1.2*ALAMSQ(3)) >*/
/* Computing MAX */
    r__1 = *p2, r__2 = alamsq[0] * (float)1.2;
    p2eff = dmax(r__1,r__2);
/*<       IF(KF.EQ.4) P2EFF=MAX(P2EFF,PMC**2) >*/
    if (*kf == 4) {
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmc;
	r__1 = p2eff, r__2 = r__3 * r__3;
	p2eff = dmax(r__1,r__2);
    }
/*<       IF(KF.EQ.5) P2EFF=MAX(P2EFF,PMB**2) >*/
    if (*kf == 5) {
/* Computing MAX */
/* Computing 2nd power */
	r__3 = pmb;
	r__1 = p2eff, r__2 = r__3 * r__3;
	p2eff = dmax(r__1,r__2);
    }
/*<       Q2EFF=MAX(Q2,P2EFF) >*/
    q2eff = dmax(*q2,p2eff);
/*<       XL=-LOG(X) >*/
    xl = -log(*x);
/* ...Find number of flavours at lower and upper scale. */
/*<       NFP=4 >*/
    nfp = 4;
/*<       IF(P2EFF.LT.PMC**2) NFP=3 >*/
/* Computing 2nd power */
    r__1 = pmc;
    if (p2eff < r__1 * r__1) {
	nfp = 3;
    }
/*<       IF(P2EFF.GT.PMB**2) NFP=5 >*/
/* Computing 2nd power */
    r__1 = pmb;
    if (p2eff > r__1 * r__1) {
	nfp = 5;
    }
/*<       NFQ=4 >*/
    nfq = 4;
/*<       IF(Q2EFF.LT.PMC**2) NFQ=3 >*/
/* Computing 2nd power */
    r__1 = pmc;
    if (q2eff < r__1 * r__1) {
	nfq = 3;
    }
/*<       IF(Q2EFF.GT.PMB**2) NFQ=5 >*/
/* Computing 2nd power */
    r__1 = pmb;
    if (q2eff > r__1 * r__1) {
	nfq = 5;
    }
/* ...Define range of flavour loop. */
/*<       IF(KF.EQ.0) THEN >*/
    if (*kf == 0) {
/*<         KFLMN=1 >*/
	kflmn = 1;
/*<         KFLMX=5 >*/
	kflmx = 5;
/*<       ELSEIF(KF.LT.0) THEN >*/
    } else if (*kf < 0) {
/*<         KFLMN=1 >*/
	kflmn = 1;
/*<         KFLMX=KFA >*/
	kflmx = kfa;
/*<       ELSE     >*/
    } else {
/*<         KFLMN=KFA >*/
	kflmn = kfa;
/*<         KFLMX=KFA >*/
	kflmx = kfa;
/*<       ENDIF >*/
    }
/* ...Loop over flavours the photon can branch into. */
/*<       DO 110 KFL=KFLMN,KFLMX >*/
    i__1 = kflmx;
    for (kfl = kflmn; kfl <= i__1; ++kfl) {
/* ...Light flavours: calculate t range and (approximate) s range. */
/*<       IF(KFL.LE.3.AND.(KFL.EQ.1.OR.KFL.EQ.KF)) THEN >*/
	if (kfl <= 3 && (kfl == 1 || kfl == *kf)) {
/*<         TDIFF=LOG(Q2EFF/P2EFF) >*/
	    tdiff = log(q2eff / p2eff);
/*<        >*/
	    s = (float)6. / ((float)33. - nfq * (float)2.) * log(log(q2eff / 
		    alamsq[nfq - 3]) / log(p2eff / alamsq[nfq - 3]));
/*<         IF(NFQ.GT.NFP) THEN >*/
	    if (nfq > nfp) {
/*<           Q2DIV=PMB**2 >*/
/* Computing 2nd power */
		r__1 = pmb;
		q2div = r__1 * r__1;
/*<           IF(NFQ.EQ.4) Q2DIV=PMC**2 >*/
		if (nfq == 4) {
/* Computing 2nd power */
		    r__1 = pmc;
		    q2div = r__1 * r__1;
		}
/*<        >*/
		snfq = (float)6. / ((float)33. - nfq * (float)2.) * log(log(
			q2div / alamsq[nfq - 3]) / log(p2eff / alamsq[nfq - 3]
			));
/*<        >*/
		snfp = (float)6. / ((float)33. - (nfq - 1) * (float)2.) * log(
			log(q2div / alamsq[nfq - 4]) / log(p2eff / alamsq[nfq 
			- 4]));
/*<           S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ) >*/
		s += log(q2div / p2eff) / log(q2eff / p2eff) * (snfp - snfq);
/*<         ENDIF >*/
	    }
/*<         IF(NFQ.EQ.5.AND.NFP.EQ.3) THEN >*/
	    if (nfq == 5 && nfp == 3) {
/*<           Q2DIV=PMC**2 >*/
/* Computing 2nd power */
		r__1 = pmc;
		q2div = r__1 * r__1;
/*<        >*/
		snf4 = log(log(q2div / alamsq[1]) / log(p2eff / alamsq[1])) * 
			(float).23999999999999999;
/*<        >*/
		snf3 = log(log(q2div / alamsq[0]) / log(p2eff / alamsq[0])) * 
			(float).22222222222222221;
/*<           S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNF3-SNF4) >*/
		s += log(q2div / p2eff) / log(q2eff / p2eff) * (snf3 - snf4);
/*<         ENDIF >*/
	    }
/* ...u and s quark do not need a separate treatment when d has been done. */
/*<       ELSEIF(KFL.EQ.2.OR.KFL.EQ.3) THEN >*/
	} else if (kfl == 2 || kfl == 3) {
/* ...Charm: as above, but only include range above c threshold. */
/*<       ELSEIF(KFL.EQ.4) THEN   >*/
	} else if (kfl == 4) {
/*<         IF(Q2.LE.PMC**2) GOTO 110 >*/
/* Computing 2nd power */
	    r__1 = pmc;
	    if (*q2 <= r__1 * r__1) {
		goto L110;
	    }
/*<         P2EFF=MAX(P2EFF,PMC**2) >*/
/* Computing MAX */
/* Computing 2nd power */
	    r__3 = pmc;
	    r__1 = p2eff, r__2 = r__3 * r__3;
	    p2eff = dmax(r__1,r__2);
/*<         Q2EFF=MAX(Q2EFF,P2EFF) >*/
	    q2eff = dmax(q2eff,p2eff);
/*<         TDIFF=LOG(Q2EFF/P2EFF) >*/
	    tdiff = log(q2eff / p2eff);
/*<        >*/
	    s = (float)6. / ((float)33. - nfq * (float)2.) * log(log(q2eff / 
		    alamsq[nfq - 3]) / log(p2eff / alamsq[nfq - 3]));
/*<         IF(NFQ.EQ.5.AND.NFP.EQ.4) THEN >*/
	    if (nfq == 5 && nfp == 4) {
/*<           Q2DIV=PMB**2 >*/
/* Computing 2nd power */
		r__1 = pmb;
		q2div = r__1 * r__1;
/*<        >*/
		snfq = (float)6. / ((float)33. - nfq * (float)2.) * log(log(
			q2div / alamsq[nfq - 3]) / log(p2eff / alamsq[nfq - 3]
			));
/*<        >*/
		snfp = (float)6. / ((float)33. - (nfq - 1) * (float)2.) * log(
			log(q2div / alamsq[nfq - 4]) / log(p2eff / alamsq[nfq 
			- 4]));
/*<           S=S+(LOG(Q2DIV/P2EFF)/LOG(Q2EFF/P2EFF))*(SNFP-SNFQ) >*/
		s += log(q2div / p2eff) / log(q2eff / p2eff) * (snfp - snfq);
/*<         ENDIF >*/
	    }
/* ...Bottom: as above, but only include range above b threshold. */
/*<       ELSEIF(KFL.EQ.5) THEN   >*/
	} else if (kfl == 5) {
/*<         IF(Q2.LE.PMB**2) GOTO 110 >*/
/* Computing 2nd power */
	    r__1 = pmb;
	    if (*q2 <= r__1 * r__1) {
		goto L110;
	    }
/*<         P2EFF=MAX(P2EFF,PMB**2) >*/
/* Computing MAX */
/* Computing 2nd power */
	    r__3 = pmb;
	    r__1 = p2eff, r__2 = r__3 * r__3;
	    p2eff = dmax(r__1,r__2);
/*<         Q2EFF=MAX(Q2,P2EFF) >*/
	    q2eff = dmax(*q2,p2eff);
/*<         TDIFF=LOG(Q2EFF/P2EFF) >*/
	    tdiff = log(q2eff / p2eff);
/*<        >*/
	    s = (float)6. / ((float)33. - nfq * (float)2.) * log(log(q2eff / 
		    alamsq[nfq - 3]) / log(p2eff / alamsq[nfq - 3]));
/*<       ENDIF >*/
	}
/* ...Evaluate flavour-dependent prefactor (charge^2 etc.). */
/*<       CHSQ=1./9. >*/
	chsq = (float).1111111111111111;
/*<       IF(KFL.EQ.2.OR.KFL.EQ.4) CHSQ=4./9. >*/
	if (kfl == 2 || kfl == 4) {
	    chsq = (float).44444444444444442;
	}
/*<       FAC=AEM2PI*2.*CHSQ*TDIFF >*/
	fac = aem2pi * (float)2. * chsq * tdiff;
/* ...Evaluate parton distributions (normalized to unit momentum sum). */
/*<       IF(KFL.EQ.1.OR.KFL.EQ.4.OR.KFL.EQ.5.OR.KFL.EQ.KF) THEN >*/
	if (kfl == 1 || kfl == 4 || kfl == 5 || kfl == *kf) {
/*<        >*/
/* Computing 2nd power */
	    r__1 = s;
/* Computing 2nd power */
	    r__2 = s;
/* Computing 2nd power */
	    r__3 = *x;
/* Computing 2nd power */
	    r__4 = s;
/* Computing 2nd power */
	    r__5 = s;
/* Computing 2nd power */
	    r__6 = (float)1. - *x;
/* Computing 2nd power */
	    r__7 = s;
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) ((float)1. / (s * (float).58 + (float)1.));
/* Computing 2nd power */
	    r__8 = *x;
	    d__3 = (doublereal) ((float)1. - r__8 * r__8);
	    d__4 = (doublereal) (s * (float)2.5 / (s * (float)10. + (float)1.)
		    );
	    xval = ((s * (float)2.49 + (float)1.5 + r__1 * r__1 * (float)26.9)
		     / (r__2 * r__2 * (float)32.3 + (float)1.) * (r__3 * r__3)
		     + ((float)1.5 - s * (float).49 + r__4 * r__4 * (float)
		    7.83) / (r__5 * r__5 * (float)7.68 + (float)1.) * (r__6 * 
		    r__6) + s * (float)1.5 / ((float)1. - s * (float)3.2 + 
		    r__7 * r__7 * (float)7.) * *x * ((float)1. - *x)) * 
		    pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4);
/*<        >*/
/* Computing 2nd power */
	    r__1 = s;
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float)-1.67 / (s * (float)2. + (float)
		    1.));
/* Computing 2nd power */
	    r__2 = *x;
	    d__3 = (doublereal) ((float)1. - r__2 * r__2);
	    d__4 = (doublereal) (s * (float)1.2);
/* Computing 2nd power */
	    r__3 = *x;
	    xglu = s * (float)2. / (s * (float)4. + (float)1. + r__1 * r__1 * 
		    (float)7.) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4) *
		     ((r__3 * r__3 * (float)4. + *x * (float)7. + (float)4.) *
		     ((float)1. - *x) / (float)3. - *x * (float)2. * (*x + (
		    float)1.) * xl);
/*<        >*/
/* Computing 2nd power */
	    r__1 = s;
/* Computing 2nd power */
	    r__2 = s;
/* Computing 3rd power */
	    r__3 = s;
	    d__1 = (doublereal) (*x);
	    d__2 = (doublereal) (s * (float)-1.18 / (s * (float)1.22 + (float)
		    1.));
	    d__3 = (doublereal) ((float)1. - *x);
	    d__4 = (doublereal) (s * (float)1.2);
/* Computing 2nd power */
	    r__4 = *x;
/* Computing 2nd power */
	    r__5 = *x;
/* Computing 2nd power */
	    r__6 = xl;
	    xsea = r__1 * r__1 * (float).333 / (s * (float)4.9 + (float)1. + 
		    r__2 * r__2 * (float)4.69 + r__3 * (r__3 * r__3) * (float)
		    21.4) * pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4) * (((
		    float)8. - *x * (float)73. + r__4 * r__4 * (float)62.) * (
		    (float)1. - *x) / (float)9. + ((float)3. - r__5 * r__5 * (
		    float)8. / (float)3.) * *x * xl + (*x * (float)2. - (
		    float)1.) * *x * (r__6 * r__6));
/* ...Threshold factors for c and b sea. */
/*<         SLL=LOG(LOG(Q2EFF/ALAM**2)/LOG(P2EFF/ALAM**2)) >*/
/* Computing 2nd power */
	    r__1 = *alam;
/* Computing 2nd power */
	    r__2 = *alam;
	    sll = log(log(q2eff / (r__1 * r__1)) / log(p2eff / (r__2 * r__2)))
		    ;
/*<         XCHM=0.     >*/
	    xchm = (float)0.;
/*<         IF(Q2.GT.PMC**2.AND.Q2.GT.1.001*P2EFF) THEN >*/
/* Computing 2nd power */
	    r__1 = pmc;
	    if (*q2 > r__1 * r__1 && *q2 > p2eff * (float)1.001) {
/*<           SCH=MAX(0.,LOG(LOG(PMC**2/ALAM**2)/LOG(P2EFF/ALAM**2)))   >*/
/* Computing MAX */
/* Computing 2nd power */
		r__3 = pmc;
/* Computing 2nd power */
		r__4 = *alam;
/* Computing 2nd power */
		r__5 = *alam;
		r__1 = (float)0., r__2 = log(log(r__3 * r__3 / (r__4 * r__4)) 
			/ log(p2eff / (r__5 * r__5)));
		sch = dmax(r__1,r__2);
/*<           XCHM=XSEA*(1.-(SCH/SLL)**3) >*/
/* Computing 3rd power */
		r__1 = sch / sll;
		xchm = xsea * ((float)1. - r__1 * (r__1 * r__1));
/*<         ENDIF      >*/
	    }
/*<         XBOT=0. >*/
	    xbot = (float)0.;
/*<         IF(Q2.GT.PMB**2.AND.Q2.GT.1.001*P2EFF) THEN >*/
/* Computing 2nd power */
	    r__1 = pmb;
	    if (*q2 > r__1 * r__1 && *q2 > p2eff * (float)1.001) {
/*<           SBT=MAX(0.,LOG(LOG(PMB**2/ALAM**2)/LOG(P2EFF/ALAM**2)))   >*/
/* Computing MAX */
/* Computing 2nd power */
		r__3 = pmb;
/* Computing 2nd power */
		r__4 = *alam;
/* Computing 2nd power */
		r__5 = *alam;
		r__1 = (float)0., r__2 = log(log(r__3 * r__3 / (r__4 * r__4)) 
			/ log(p2eff / (r__5 * r__5)));
		sbt = dmax(r__1,r__2);
/*<           XBOT=XSEA*(1.-(SBT/SLL)**3) >*/
/* Computing 3rd power */
		r__1 = sbt / sll;
		xbot = xsea * ((float)1. - r__1 * (r__1 * r__1));
/*<         ENDIF    >*/
	    }
/*<       ENDIF >*/
	}
/* ...Add contribution of each valence flavour. */
/*<       XPGA(0)=XPGA(0)+FAC*XGLU  >*/
	xpga[0] += fac * xglu;
/*<       XPGA(1)=XPGA(1)+FAC*XSEA >*/
	xpga[1] += fac * xsea;
/*<       XPGA(2)=XPGA(2)+FAC*XSEA >*/
	xpga[2] += fac * xsea;
/*<       XPGA(3)=XPGA(3)+FAC*XSEA >*/
	xpga[3] += fac * xsea;
/*<       XPGA(4)=XPGA(4)+FAC*XCHM >*/
	xpga[4] += fac * xchm;
/*<       XPGA(5)=XPGA(5)+FAC*XBOT >*/
	xpga[5] += fac * xbot;
/*<       XPGA(KFL)=XPGA(KFL)+FAC*XVAL >*/
	xpga[kfl] += fac * xval;
/*<       VXPGA(KFL)=VXPGA(KFL)+FAC*XVAL >*/
	vxpga[kfl] += fac * xval;
/*<   110 CONTINUE >*/
L110:
	;
    }
/*<       DO 120 KFL=1,5 >*/
    for (kfl = 1; kfl <= 5; ++kfl) {
/*<       XPGA(-KFL)=XPGA(KFL) >*/
	xpga[-kfl] = xpga[kfl];
/*<       VXPGA(-KFL)=VXPGA(KFL) >*/
	vxpga[-kfl] = vxpga[kfl];
/*<   120 CONTINUE >*/
/* L120: */
    }
/*<       RETURN >*/
    return 0;
/*<       END  >*/
} /* sasano_ */

/* ********************************************************************* */
/*<       SUBROUTINE SASBEH(KF,X,Q2,P2,PM2,XPBH) >*/
/* Subroutine */ int sasbeh_(integer *kf, real *x, real *q2, real *p2, real *
	pm2, real *xpbh)
{
    /* Initialized data */

    static real aem2pi = (float).0011614;

    /* System generated locals */
    real r__1, r__2, r__3, r__4;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal);

    /* Local variables */
    static real w2, xbi, xbl, rmq, rpq, beta, rpbe, chsq, beta2, sigbh, 
	    rpbesn;

/* ...Purpose: to evaluate the Bethe-Heitler cross section for */
/* ...heavy flavour production. */
/*<       DATA AEM2PI/0.0011614/ >*/
/* ...Reset output. */
/*<       XPBH=0. >*/
    *xpbh = (float)0.;
/*<       SIGBH=0. >*/
    sigbh = (float)0.;
/* ...Check kinematics limits. */
/*<       IF(X.GE.Q2/(4.*PM2+Q2+P2)) RETURN            >*/
    if (*x >= *q2 / (*pm2 * (float)4. + *q2 + *p2)) {
	return 0;
    }
/*<       W2=Q2*(1.-X)/X-P2 >*/
    w2 = *q2 * ((float)1. - *x) / *x - *p2;
/*<       BETA2=1.-4.*PM2/W2 >*/
    beta2 = (float)1. - *pm2 * (float)4. / w2;
/*<       IF(BETA2.LT.1E-10) RETURN >*/
    if (beta2 < (float)1e-10) {
	return 0;
    }
/*<       BETA=SQRT(BETA2) >*/
    beta = sqrt(beta2);
/*<       RMQ=4.*PM2/Q2 >*/
    rmq = *pm2 * (float)4. / *q2;
/* ...Simple case: P2 = 0. */
/*<       IF(P2.LT.1E-4) THEN >*/
    if (*p2 < (float)1e-4) {
/*<         IF(BETA.LT.0.99) THEN >*/
	if (beta < (float).99) {
/*<           XBL=LOG((1.+BETA)/(1.-BETA)) >*/
	    xbl = log((beta + (float)1.) / ((float)1. - beta));
/*<         ELSE >*/
	} else {
/*<           XBL=LOG((1.+BETA)**2*W2/(4.*PM2)) >*/
/* Computing 2nd power */
	    r__1 = beta + (float)1.;
	    xbl = log(r__1 * r__1 * w2 / (*pm2 * (float)4.));
/*<         ENDIF  >*/
	}
/*<        >*/
/* Computing 2nd power */
	r__1 = *x;
/* Computing 2nd power */
	r__2 = (float)1. - *x;
/* Computing 2nd power */
	r__3 = rmq;
/* Computing 2nd power */
	r__4 = *x;
	sigbh = beta * (*x * (float)8. * ((float)1. - *x) - (float)1. - rmq * 
		*x * ((float)1. - *x)) + xbl * (r__1 * r__1 + r__2 * r__2 + 
		rmq * *x * ((float)1. - *x * (float)3.) - r__3 * r__3 * (
		float).5 * (r__4 * r__4));
/* ...Complicated case: P2 > 0, based on approximation of */
/* ...C.T. Hill and G.G. Ross, Nucl. Phys. B148 (1979) 373 */
/*<       ELSE >*/
    } else {
/*<         RPQ=1.-4.*X**2*P2/Q2 >*/
/* Computing 2nd power */
	r__1 = *x;
	rpq = (float)1. - r__1 * r__1 * (float)4. * *p2 / *q2;
/*<         IF(RPQ.GT.1E-10) THEN >*/
	if (rpq > (float)1e-10) {
/*<           RPBE=SQRT(RPQ*BETA2) >*/
	    rpbe = sqrt(rpq * beta2);
/*<           IF(RPBE.LT.0.99) THEN >*/
	    if (rpbe < (float).99) {
/*<             XBL=LOG((1.+RPBE)/(1.-RPBE)) >*/
		xbl = log((rpbe + (float)1.) / ((float)1. - rpbe));
/*<             XBI=2.*RPBE/(1.-RPBE**2) >*/
/* Computing 2nd power */
		r__1 = rpbe;
		xbi = rpbe * (float)2. / ((float)1. - r__1 * r__1);
/*<           ELSE >*/
	    } else {
/*<             RPBESN=4.*PM2/W2+(4.*X**2*P2/Q2)*BETA2 >*/
/* Computing 2nd power */
		r__1 = *x;
		rpbesn = *pm2 * (float)4. / w2 + r__1 * r__1 * (float)4. * *
			p2 / *q2 * beta2;
/*<             XBL=LOG((1.+RPBE)**2/RPBESN) >*/
/* Computing 2nd power */
		r__1 = rpbe + (float)1.;
		xbl = log(r__1 * r__1 / rpbesn);
/*<             XBI=2.*RPBE/RPBESN >*/
		xbi = rpbe * (float)2. / rpbesn;
/*<           ENDIF >*/
	    }
/*<        >*/
/* Computing 2nd power */
	    r__1 = *x;
/* Computing 2nd power */
	    r__2 = (float)1. - *x;
/* Computing 2nd power */
	    r__3 = rmq;
/* Computing 2nd power */
	    r__4 = *x;
	    sigbh = beta * (*x * (float)6. * ((float)1. - *x) - (float)1.) + 
		    xbl * (r__1 * r__1 + r__2 * r__2 + rmq * *x * ((float)1. 
		    - *x * (float)3.) - r__3 * r__3 * (float).5 * (r__4 * 
		    r__4)) + xbi * (*x * (float)2. / *q2) * (*pm2 * *x * ((
		    float)2. - rmq) - *p2 * *x);
/*<         ENDIF                 >*/
	}
/*<       ENDIF >*/
    }
/* ...Multiply by charge-squared etc. to get parton distribution. */
/*<       CHSQ=1./9. >*/
    chsq = (float).1111111111111111;
/*<       IF(IABS(KF).EQ.2.OR.IABS(KF).EQ.4) CHSQ=4./9. >*/
    if (abs(*kf) == 2 || abs(*kf) == 4) {
	chsq = (float).44444444444444442;
    }
/*<       XPBH=3.*CHSQ*AEM2PI*X*SIGBH        >*/
    *xpbh = chsq * (float)3. * aem2pi * *x * sigbh;
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* sasbeh_ */

/* ********************************************************************* */
/*<        SUBROUTINE SASDIR(X,Q2,P2,Q02,XPGA) >*/
/* Subroutine */ int sasdir_(real *x, real *q2, real *p2, real *q02, real *
	xpga)
{
    /* Initialized data */

    static real aem2pi = (float).0011614;

    /* System generated locals */
    real r__1, r__2;

    /* Builtin functions */
    double log(doublereal);

    /* Local variables */
    static integer kf, kfl;
    static real cgam, xtmp;

/* ...Purpose: to evaluate the direct contribution, i.e. the C^gamma term, */
/* ...as needed in MSbar parametrizations. */
/*<       DIMENSION XPGA(-6:6) >*/
/*<       DATA PMC/1.3/, PMB/4.6/, AEM2PI/0.0011614/ >*/
    /* Parameter adjustments */
    xpga -= -6;

    /* Function Body */
/* ...Reset output. */
/*<       DO 100 KFL=-6,6 >*/
    for (kfl = -6; kfl <= 6; ++kfl) {
/*<       XPGA(KFL)=0. >*/
	xpga[kfl] = (float)0.;
/*<   100 CONTINUE >*/
/* L100: */
    }
/* ...Evaluate common x-dependent expression. */
/*<       XTMP = (X**2+(1.-X)**2) * (-LOG(X)) - 1. >*/
/* Computing 2nd power */
    r__1 = *x;
/* Computing 2nd power */
    r__2 = (float)1. - *x;
    xtmp = (r__1 * r__1 + r__2 * r__2) * (-log(*x)) - (float)1.;
/*<       CGAM = 3.*AEM2PI*X * (XTMP*(1.+P2/(P2+Q02)) + 6.*X*(1.-X)) >*/
    cgam = aem2pi * (float)3. * *x * (xtmp * (*p2 / (*p2 + *q02) + (float)1.) 
	    + *x * (float)6. * ((float)1. - *x));
/* ...d, u, s part by simple charge factor. */
/*<       XPGA(1)=(1./9.)*CGAM >*/
    xpga[1] = cgam * (float).1111111111111111;
/*<       XPGA(2)=(4./9.)*CGAM >*/
    xpga[2] = cgam * (float).44444444444444442;
/*<       XPGA(3)=(1./9.)*CGAM       >*/
    xpga[3] = cgam * (float).1111111111111111;
/* ...Also fill for antiquarks. */
/*<       DO 110 KF=1,5 >*/
    for (kf = 1; kf <= 5; ++kf) {
/*<       XPGA(-KF)=XPGA(KF) >*/
	xpga[-kf] = xpga[kf];
/*<   110 CONTINUE >*/
/* L110: */
    }
/*<       RETURN >*/
    return 0;
/*<       END >*/
} /* sasdir_ */

#ifdef __cplusplus
	}
#endif
