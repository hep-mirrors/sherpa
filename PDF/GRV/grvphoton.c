/* grvphoton.f -- translated by f2c (version 20200916).
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

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                 * */
/*      G R V - P H O T O N - P A R A M E T R I Z A T I O N S      * */
/*                                                                 * */
/*                 FOR A DETAILED EXPLANATION SEE :                * */
/*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/31             * */
/*             PUBLISHED IN PHYS. REV. D46 (1992) 1973             * */
/*                                                                 * */
/*    THE OUTPUT IS ALWAYS   1./ ALPHA(EM) * X * PARTON DENSITY    * */
/*                                                                 * */
/*   THE PARAMETRIZATIONS ARE FITTED TO THE PARTON DISTRIBUTIONS   * */
/*   FOR Q ** 2 BETWEEN MU ** 2 (=  0.25 / 0.30  GEV ** 2  IN LO   * */
/*   / HO) AND  1.E6 GEV ** 2  AND FOR X BETWEEN  1.E-5  AND  1.   * */
/*                                                                 * */
/*              HEAVY QUARK THRESHOLDS  Q(H) = M(H) :              * */
/*         M(C)  =  1.5,  M(B)  =  4.5,  M(T)  =  100  GEV         * */
/*                                                                 * */
/*      CORRESPONDING LAMBDA(F) VALUES FOR F ACTIVE FLAVOURS :     * */
/*      LO :   LAMBDA(3)  =  0.232,   LAMBDA(4)  =  0.200,         * */
/*             LAMBDA(5)  =  0.153,   LAMBDA(6)  =  0.082  GEV     * */
/*      HO :   LAMBDA(3)  =  0.248,   LAMBDA(4)  =  0.200,         * */
/*             LAMBDA(5)  =  0.131,   LAMBDA(6)  =  0.053  GEV     * */
/*                                                                 * */
/*      HO DISTRIBUTIONS REFER TO THE DIS(GAMMA) SCHEME, SEE :     * */
/*              M. GLUECK, E.REYA, A.VOGT: DO-TH 91/26             * */
/*              PUBLISHED IN PHYS. REV. D45 (1992) 3986            * */
/*                                                                 * */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*<        SUBROUTINE GRVGLO (X, Q2, UL, DL, SL, CL, BL, GL) >*/
/* Subroutine */ int grvglo_(real *x, real *q2, real *ul, real *dl, real *sl, 
	real *cl, real *bl, real *gl)
{
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static real c__, d__, e, s, s2, be, ag, bg, ak, al, bk, es, sf;
    extern doublereal fs_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static real ss, mu2, lam2;
    extern doublereal fone_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *);

/*<        IMPLICIT REAL (A - Z) >*/
/*<        MU2  = 0.25 >*/
    mu2 = (float).25;
/*<        LAM2 = 0.232 * 0.232 >*/
    lam2 = (float).053824000000000004;
/*<        S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2)) >*/
    s = log(log(*q2 / lam2) / log(mu2 / lam2));
/*<        SS = SQRT (S) >*/
    ss = sqrt(s);
/*<        S2 = S * S >*/
    s2 = s * s;
/* ...X * U = X * UBAR : */
/*<        AL =  1.717 >*/
    al = (float)1.717;
/*<        BE =  0.641 >*/
    be = (float).641;
/*<        AK =  0.500 - 0.176 * S >*/
    ak = (float).5 - s * (float).176;
/*<        BK = 15.00  - 5.687 * SS - 0.552 * S2 >*/
    bk = (float)15. - ss * (float)5.687 - s2 * (float).552;
/*<        AG =  0.235 + 0.046 * SS >*/
    ag = ss * (float).046 + (float).235;
/*<        BG =  0.082 - 0.051 * S  + 0.168 * S2 >*/
    bg = (float).082 - s * (float).051 + s2 * (float).168;
/*<        C  =   0.0  + 0.459 * S >*/
    c__ = s * (float).459 + (float)0.;
/*<        D  =  0.354 - 0.061 * S >*/
    d__ = (float).354 - s * (float).061;
/*<        E  =  4.899 + 1.678 * S >*/
    e = s * (float)1.678 + (float)4.899;
/*<        ES =  2.046 + 1.389 * S >*/
    es = s * (float)1.389 + (float)2.046;
/*<        UL = FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *ul = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * D = X * DBAR : */
/*<        AL =  1.549 >*/
    al = (float)1.549;
/*<        BE =  0.782 >*/
    be = (float).782;
/*<        AK =  0.496 + 0.026 * S >*/
    ak = s * (float).026 + (float).496;
/*<        BK =  0.685 - 0.580 * SS + 0.608 * S2 >*/
    bk = (float).685 - ss * (float).58 + s2 * (float).608;
/*<        AG =  0.233 + 0.302 * S >*/
    ag = s * (float).302 + (float).233;
/*<        BG =   0.0  - 0.818 * S  + 0.198 * S2 >*/
    bg = (float)0. - s * (float).818 + s2 * (float).198;
/*<        C  =  0.114 + 0.154 * S >*/
    c__ = s * (float).154 + (float).114;
/*<        D  =  0.405 - 0.195 * S  + 0.046 * S2 >*/
    d__ = (float).405 - s * (float).195 + s2 * (float).046;
/*<        E  =  4.807 + 1.226 * S >*/
    e = s * (float)1.226 + (float)4.807;
/*<        ES =  2.166 + 0.664 * S >*/
    es = s * (float).664 + (float)2.166;
/*<        DL  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *dl = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * G : */
/*<        AL =  0.676 >*/
    al = (float).676;
/*<        BE =  1.089 >*/
    be = (float)1.089;
/*<        AK =  0.462 - 0.524 * SS >*/
    ak = (float).462 - ss * (float).524;
/*<        BK =  5.451              - 0.804 * S2 >*/
    bk = (float)5.451 - s2 * (float).804;
/*<        AG =  0.535 - 0.504 * SS + 0.288 * S2 >*/
    ag = (float).535 - ss * (float).504 + s2 * (float).288;
/*<        BG =  0.364 - 0.520 * S >*/
    bg = (float).364 - s * (float).52;
/*<        C  = -0.323              + 0.115 * S2 >*/
    c__ = s2 * (float).115 - (float).323;
/*<        D  =  0.233 + 0.790 * S  - 0.139 * S2 >*/
    d__ = s * (float).79 + (float).233 - s2 * (float).139;
/*<        E  =  0.893 + 1.968 * S >*/
    e = s * (float)1.968 + (float).893;
/*<        ES =  3.432 + 0.392 * S >*/
    es = s * (float).392 + (float)3.432;
/*<        GL =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *gl = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * S = X * SBAR : */
/*<        SF =   0.0 >*/
    sf = (float)0.;
/*<        AL =  1.609 >*/
    al = (float)1.609;
/*<        BE =  0.962 >*/
    be = (float).962;
/*<        AK =  0.470              - 0.099 * S2 >*/
    ak = (float).47 - s2 * (float).099;
/*<        BK =  3.246 >*/
    bk = (float)3.246;
/*<        AG =  0.121 - 0.068 * SS >*/
    ag = (float).121 - ss * (float).068;
/*<        BG = -0.090 + 0.074 * S >*/
    bg = s * (float).074 - (float).09;
/*<        C  =  0.062 + 0.034 * S >*/
    c__ = s * (float).034 + (float).062;
/*<        D  =   0.0  + 0.226 * S  - 0.060 * S2 >*/
    d__ = s * (float).226 + (float)0. - s2 * (float).06;
/*<        E  =  4.288 + 1.707 * S >*/
    e = s * (float)1.707 + (float)4.288;
/*<        ES =  2.122 + 0.656 * S >*/
    es = s * (float).656 + (float)2.122;
/*<        SL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *sl = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * C = X * CBAR : */
/*<        SF =  0.888 >*/
    sf = (float).888;
/*<        AL =  0.970 >*/
    al = (float).97;
/*<        BE =  0.545 >*/
    be = (float).545;
/*<        AK =  1.254 - 0.251 * S >*/
    ak = (float)1.254 - s * (float).251;
/*<        BK =  3.932              - 0.327 * S2 >*/
    bk = (float)3.932 - s2 * (float).327;
/*<        AG =  0.658 + 0.202 * S >*/
    ag = s * (float).202 + (float).658;
/*<        BG = -0.699 >*/
    bg = (float)-.699;
/*<        C  =  0.965 >*/
    c__ = (float).965;
/*<        D  =   0.0  + 0.141 * S  - 0.027 * S2 >*/
    d__ = s * (float).141 + (float)0. - s2 * (float).027;
/*<        E  =  4.911 + 0.969 * S >*/
    e = s * (float).969 + (float)4.911;
/*<        ES =  2.796 + 0.952 * S >*/
    es = s * (float).952 + (float)2.796;
/*<        CL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *cl = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * B = X * BBAR : */
/*<        SF =  1.351 >*/
    sf = (float)1.351;
/*<        AL =  1.016 >*/
    al = (float)1.016;
/*<        BE =  0.338 >*/
    be = (float).338;
/*<        AK =  1.961 - 0.370 * S >*/
    ak = (float)1.961 - s * (float).37;
/*<        BK =  0.923 + 0.119 * S >*/
    bk = s * (float).119 + (float).923;
/*<        AG =  0.815 + 0.207 * S >*/
    ag = s * (float).207 + (float).815;
/*<        BG = -2.275 >*/
    bg = (float)-2.275;
/*<        C  =  1.480 >*/
    c__ = (float)1.48;
/*<        D  = -0.223 + 0.173 * S >*/
    d__ = s * (float).173 - (float).223;
/*<        E  =  5.426 + 0.623 * S >*/
    e = s * (float).623 + (float)5.426;
/*<        ES =  3.819 + 0.901 * S >*/
    es = s * (float).901 + (float)3.819;
/*<        BL =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *bl = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/*<        RETURN >*/
    return 0;
/*<        END >*/
} /* grvglo_ */



/*<        SUBROUTINE GRVGHO (X, Q2, UH, DH, SH, CH, BH, GH) >*/
/* Subroutine */ int grvgho_(real *x, real *q2, real *uh, real *dh, real *sh, 
	real *ch, real *bh, real *gh)
{
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static real c__, d__, e, s, s2, be, ag, bg, al, ak, bk, es, sf;
    extern doublereal fs_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static real ss, mu2, lam2;
    extern doublereal fone_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *);

/*<        IMPLICIT REAL (A - Z) >*/
/*<        MU2  = 0.3 >*/
    mu2 = (float).3;
/*<        LAM2 = 0.248 * 0.248 >*/
    lam2 = (float).061503999999999996;
/*<        S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2)) >*/
    s = log(log(*q2 / lam2) / log(mu2 / lam2));
/*<        SS = SQRT (S) >*/
    ss = sqrt(s);
/*<        S2 = S * S >*/
    s2 = s * s;
/* ...X * U = X * UBAR : */
/*<        AL =  0.583 >*/
    al = (float).583;
/*<        BE =  0.688 >*/
    be = (float).688;
/*<        AK =  0.449 - 0.025 * S  - 0.071 * S2 >*/
    ak = (float).449 - s * (float).025 - s2 * (float).071;
/*<        BK =  5.060 - 1.116 * SS >*/
    bk = (float)5.06 - ss * (float)1.116;
/*<        AG =  0.103 >*/
    ag = (float).103;
/*<        BG =  0.319 + 0.422 * S >*/
    bg = s * (float).422 + (float).319;
/*<        C  =  1.508 + 4.792 * S  - 1.963 * S2 >*/
    c__ = s * (float)4.792 + (float)1.508 - s2 * (float)1.963;
/*<        D  =  1.075 + 0.222 * SS - 0.193 * S2 >*/
    d__ = ss * (float).222 + (float)1.075 - s2 * (float).193;
/*<        E  =  4.147 + 1.131 * S >*/
    e = s * (float)1.131 + (float)4.147;
/*<        ES =  1.661 + 0.874 * S >*/
    es = s * (float).874 + (float)1.661;
/*<        UH =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *uh = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * D = X * DBAR : */
/*<        AL =  0.591 >*/
    al = (float).591;
/*<        BE =  0.698 >*/
    be = (float).698;
/*<        AK =  0.442 - 0.132 * S  - 0.058 * S2 >*/
    ak = (float).442 - s * (float).132 - s2 * (float).058;
/*<        BK =  5.437 - 1.916 * SS >*/
    bk = (float)5.437 - ss * (float)1.916;
/*<        AG =  0.099 >*/
    ag = (float).099;
/*<        BG =  0.311 - 0.059 * S >*/
    bg = (float).311 - s * (float).059;
/*<        C  =  0.800 + 0.078 * S  - 0.100 * S2 >*/
    c__ = s * (float).078 + (float).8 - s2 * (float).1;
/*<        D  =  0.862 + 0.294 * SS - 0.184 * S2 >*/
    d__ = ss * (float).294 + (float).862 - s2 * (float).184;
/*<        E  =  4.202 + 1.352 * S >*/
    e = s * (float)1.352 + (float)4.202;
/*<        ES =  1.841 + 0.990 * S >*/
    es = s * (float).99 + (float)1.841;
/*<        DH  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *dh = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * G : */
/*<        AL =  1.161 >*/
    al = (float)1.161;
/*<        BE =  1.591 >*/
    be = (float)1.591;
/*<        AK =  0.530 - 0.742 * SS + 0.025 * S2 >*/
    ak = (float).53 - ss * (float).742 + s2 * (float).025;
/*<        BK =  5.662 >*/
    bk = (float)5.662;
/*<        AG =  0.533 - 0.281 * SS + 0.218 * S2 >*/
    ag = (float).533 - ss * (float).281 + s2 * (float).218;
/*<        BG =  0.025 - 0.518 * S  + 0.156 * S2 >*/
    bg = (float).025 - s * (float).518 + s2 * (float).156;
/*<        C  = -0.282              + 0.209 * S2 >*/
    c__ = s2 * (float).209 - (float).282;
/*<        D  =  0.107 + 1.058 * S  - 0.218 * S2 >*/
    d__ = s * (float)1.058 + (float).107 - s2 * (float).218;
/*<        E  =   0.0  + 2.704 * S >*/
    e = s * (float)2.704 + (float)0.;
/*<        ES =  3.071 - 0.378 * S >*/
    es = (float)3.071 - s * (float).378;
/*<        GH =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *gh = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * S = X * SBAR : */
/*<        SF =   0.0 >*/
    sf = (float)0.;
/*<        AL =  0.635 >*/
    al = (float).635;
/*<        BE =  0.456 >*/
    be = (float).456;
/*<        AK =  1.770 - 0.735 * SS - 0.079 * S2 >*/
    ak = (float)1.77 - ss * (float).735 - s2 * (float).079;
/*<        BK =  3.832 >*/
    bk = (float)3.832;
/*<        AG =  0.084 - 0.023 * S >*/
    ag = (float).084 - s * (float).023;
/*<        BG =  0.136 >*/
    bg = (float).136;
/*<        C  =  2.119 - 0.942 * S  + 0.063 * S2 >*/
    c__ = (float)2.119 - s * (float).942 + s2 * (float).063;
/*<        D  =  1.271 + 0.076 * S  - 0.190 * S2 >*/
    d__ = s * (float).076 + (float)1.271 - s2 * (float).19;
/*<        E  =  4.604 + 0.737 * S >*/
    e = s * (float).737 + (float)4.604;
/*<        ES =  1.641 + 0.976 * S >*/
    es = s * (float).976 + (float)1.641;
/*<        SH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *sh = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * C = X * CBAR : */
/*<        SF =  0.820 >*/
    sf = (float).82;
/*<        AL =  0.926 >*/
    al = (float).926;
/*<        BE =  0.152 >*/
    be = (float).152;
/*<        AK =  1.142 - 0.175 * S >*/
    ak = (float)1.142 - s * (float).175;
/*<        BK =  3.276 >*/
    bk = (float)3.276;
/*<        AG =  0.504 + 0.317 * S >*/
    ag = s * (float).317 + (float).504;
/*<        BG = -0.433 >*/
    bg = (float)-.433;
/*<        C  =  3.334 >*/
    c__ = (float)3.334;
/*<        D  =  0.398 + 0.326 * S  - 0.107 * S2 >*/
    d__ = s * (float).326 + (float).398 - s2 * (float).107;
/*<        E  =  5.493 + 0.408 * S >*/
    e = s * (float).408 + (float)5.493;
/*<        ES =  2.426 + 1.277 * S >*/
    es = s * (float)1.277 + (float)2.426;
/*<        CH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *ch = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * B = X * BBAR : */
/*<        SF =  1.297 >*/
    sf = (float)1.297;
/*<        AL =  0.969 >*/
    al = (float).969;
/*<        BE =  0.266 >*/
    be = (float).266;
/*<        AK =  1.953 - 0.391 * S >*/
    ak = (float)1.953 - s * (float).391;
/*<        BK =  1.657 - 0.161 * S >*/
    bk = (float)1.657 - s * (float).161;
/*<        AG =  1.076 + 0.034 * S >*/
    ag = s * (float).034 + (float)1.076;
/*<        BG = -2.015 >*/
    bg = (float)-2.015;
/*<        C  =  1.662 >*/
    c__ = (float)1.662;
/*<        D  =  0.353 + 0.016 * S >*/
    d__ = s * (float).016 + (float).353;
/*<        E  =  5.713 + 0.249 * S >*/
    e = s * (float).249 + (float)5.713;
/*<        ES =  3.456 + 0.673 * S >*/
    es = s * (float).673 + (float)3.456;
/*<        BH =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *bh = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/*<        RETURN >*/
    return 0;
/*<        END >*/
} /* grvgho_ */



/*<        SUBROUTINE GRVGH0 (X, Q2, U0, D0, S0, C0, B0, G0) >*/
/* Subroutine */ int grvgh0_(real *x, real *q2, real *u0, real *d0, real *s0, 
	real *c0, real *b0, real *g0)
{
    /* Builtin functions */
    double log(doublereal), sqrt(doublereal);

    /* Local variables */
    static real c__, d__, e, s, s2, be, ag, bg, ak, al, bk, es, sf;
    extern doublereal fs_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *, real *);
    static real ss, mu2, lam2;
    extern doublereal fone_(real *, real *, real *, real *, real *, real *, 
	    real *, real *, real *, real *, real *, real *);

/*<        IMPLICIT REAL (A - Z) >*/
/*<        MU2  = 0.3 >*/
    mu2 = (float).3;
/*<        LAM2 = 0.248 * 0.248 >*/
    lam2 = (float).061503999999999996;
/*<        S  = ALOG (ALOG(Q2/LAM2) / ALOG(MU2/LAM2)) >*/
    s = log(log(*q2 / lam2) / log(mu2 / lam2));
/*<        SS = SQRT (S) >*/
    ss = sqrt(s);
/*<        S2 = S * S >*/
    s2 = s * s;
/* ...X * U = X * UBAR : */
/*<        AL =  1.447 >*/
    al = (float)1.447;
/*<        BE =  0.848 >*/
    be = (float).848;
/*<        AK =  0.527 + 0.200 * S  - 0.107 * S2 >*/
    ak = s * (float).2 + (float).527 - s2 * (float).107;
/*<        BK =  7.106 - 0.310 * SS - 0.786 * S2 >*/
    bk = (float)7.106 - ss * (float).31 - s2 * (float).786;
/*<        AG =  0.197 + 0.533 * S >*/
    ag = s * (float).533 + (float).197;
/*<        BG =  0.062 - 0.398 * S  + 0.109 * S2 >*/
    bg = (float).062 - s * (float).398 + s2 * (float).109;
/*<        C  =          0.755 * S  - 0.112 * S2 >*/
    c__ = s * (float).755 - s2 * (float).112;
/*<        D  =  0.318 - 0.059 * S >*/
    d__ = (float).318 - s * (float).059;
/*<        E  =  4.225 + 1.708 * S >*/
    e = s * (float)1.708 + (float)4.225;
/*<        ES =  1.752 + 0.866 * S >*/
    es = s * (float).866 + (float)1.752;
/*<        U0 =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *u0 = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * D = X * DBAR : */
/*<        AL =  1.424 >*/
    al = (float)1.424;
/*<        BE =  0.770 >*/
    be = (float).77;
/*<        AK =  0.500 + 0.067 * SS - 0.055 * S2 >*/
    ak = ss * (float).067 + (float).5 - s2 * (float).055;
/*<        BK =  0.376 - 0.453 * SS + 0.405 * S2 >*/
    bk = (float).376 - ss * (float).453 + s2 * (float).405;
/*<        AG =  0.156 + 0.184 * S >*/
    ag = s * (float).184 + (float).156;
/*<        BG =   0.0  - 0.528 * S  + 0.146 * S2 >*/
    bg = (float)0. - s * (float).528 + s2 * (float).146;
/*<        C  =  0.121 + 0.092 * S >*/
    c__ = s * (float).092 + (float).121;
/*<        D  =  0.379 - 0.301 * S  + 0.081 * S2 >*/
    d__ = (float).379 - s * (float).301 + s2 * (float).081;
/*<        E  =  4.346 + 1.638 * S >*/
    e = s * (float)1.638 + (float)4.346;
/*<        ES =  1.645 + 1.016 * S >*/
    es = s * (float)1.016 + (float)1.645;
/*<        D0  =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *d0 = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * G : */
/*<        AL =  0.661 >*/
    al = (float).661;
/*<        BE =  0.793 >*/
    be = (float).793;
/*<        AK =  0.537 - 0.600 * SS >*/
    ak = (float).537 - ss * (float).6;
/*<        BK =  6.389              - 0.953 * S2 >*/
    bk = (float)6.389 - s2 * (float).953;
/*<        AG =  0.558 - 0.383 * SS + 0.261 * S2 >*/
    ag = (float).558 - ss * (float).383 + s2 * (float).261;
/*<        BG =   0.0  - 0.305 * S >*/
    bg = (float)0. - s * (float).305;
/*<        C  = -0.222              + 0.078 * S2 >*/
    c__ = s2 * (float).078 - (float).222;
/*<        D  =  0.153 + 0.978 * S  - 0.209 * S2 >*/
    d__ = s * (float).978 + (float).153 - s2 * (float).209;
/*<        E  =  1.429 + 1.772 * S >*/
    e = s * (float)1.772 + (float)1.429;
/*<        ES =  3.331 + 0.806 * S >*/
    es = s * (float).806 + (float)3.331;
/*<        G0 =  FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *g0 = fone_(x, &s, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * S = X * SBAR : */
/*<        SF =   0.0 >*/
    sf = (float)0.;
/*<        AL =  1.578 >*/
    al = (float)1.578;
/*<        BE =  0.863 >*/
    be = (float).863;
/*<        AK =  0.622 + 0.332 * S  - 0.300 * S2 >*/
    ak = s * (float).332 + (float).622 - s2 * (float).3;
/*<        BK =  2.469 >*/
    bk = (float)2.469;
/*<        AG =  0.211 - 0.064 * SS - 0.018 * S2 >*/
    ag = (float).211 - ss * (float).064 - s2 * (float).018;
/*<        BG = -0.215 + 0.122 * S >*/
    bg = s * (float).122 - (float).215;
/*<        C  =  0.153 >*/
    c__ = (float).153;
/*<        D  =   0.0  + 0.253 * S  - 0.081 * S2 >*/
    d__ = s * (float).253 + (float)0. - s2 * (float).081;
/*<        E  =  3.990 + 2.014 * S >*/
    e = s * (float)2.014 + (float)3.99;
/*<        ES =  1.720 + 0.986 * S >*/
    es = s * (float).986 + (float)1.72;
/*<        S0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *s0 = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * C = X * CBAR : */
/*<        SF =  0.820 >*/
    sf = (float).82;
/*<        AL =  0.929 >*/
    al = (float).929;
/*<        BE =  0.381 >*/
    be = (float).381;
/*<        AK =  1.228 - 0.231 * S >*/
    ak = (float)1.228 - s * (float).231;
/*<        BK =  3.806             - 0.337 * S2 >*/
    bk = (float)3.806 - s2 * (float).337;
/*<        AG =  0.932 + 0.150 * S >*/
    ag = s * (float).15 + (float).932;
/*<        BG = -0.906 >*/
    bg = (float)-.906;
/*<        C  =  1.133 >*/
    c__ = (float)1.133;
/*<        D  =   0.0  + 0.138 * S  - 0.028 * S2 >*/
    d__ = s * (float).138 + (float)0. - s2 * (float).028;
/*<        E  =  5.588 + 0.628 * S >*/
    e = s * (float).628 + (float)5.588;
/*<        ES =  2.665 + 1.054 * S >*/
    es = s * (float)1.054 + (float)2.665;
/*<        C0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *c0 = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/* ...X * B = X * BBAR : */
/*<        SF =  1.297 >*/
    sf = (float)1.297;
/*<        AL =  0.970 >*/
    al = (float).97;
/*<        BE =  0.207 >*/
    be = (float).207;
/*<        AK =  1.719 - 0.292 * S >*/
    ak = (float)1.719 - s * (float).292;
/*<        BK =  0.928 + 0.096 * S >*/
    bk = s * (float).096 + (float).928;
/*<        AG =  0.845 + 0.178 * S >*/
    ag = s * (float).178 + (float).845;
/*<        BG = -2.310 >*/
    bg = (float)-2.31;
/*<        C  =  1.558 >*/
    c__ = (float)1.558;
/*<        D  = -0.191 + 0.151 * S >*/
    d__ = s * (float).151 - (float).191;
/*<        E  =  6.089 + 0.282 * S >*/
    e = s * (float).282 + (float)6.089;
/*<        ES =  3.379 + 1.062 * S >*/
    es = s * (float)1.062 + (float)3.379;
/*<        B0 =  FS (X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
    *b0 = fs_(x, &s, &sf, &al, &be, &ak, &bk, &ag, &bg, &c__, &d__, &e, &es);
/*<        RETURN >*/
    return 0;
/*<        END >*/
} /* grvgh0_ */



/*       FUNCTION F(X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) */
/*       IMPLICIT REAL (A - Z) */
/*       SX = SQRT (X) */
/*       LX = ALOG (1./X) */
/*       F  = (X**AK * (AG + BG * SX + C * X**BK)  +  S**AL */
/*     1       * EXP (-E + SQRT (ES * S**BE * LX))) * (1.- X)**D */
/*       RETURN */
/*       END */

/*      renamed function F to FONE */
/*<        FUNCTION FONE (X, S, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
doublereal fone_(real *x, real *s, real *al, real *be, real *ak, real *bk, 
	real *ag, real *bg, real *c__, real *d__, real *e, real *es)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static real lx, sx;

/*<        IMPLICIT REAL (A - Z) >*/
/*<        SX = SQRT (X) >*/
    sx = sqrt(*x);
/*<        LX = ALOG (1./X) >*/
    lx = log((float)1. / *x);
/*<        >*/
    d__1 = (doublereal) (*x);
    d__2 = (doublereal) (*ak);
    d__3 = (doublereal) (*x);
    d__4 = (doublereal) (*bk);
    d__5 = (doublereal) (*s);
    d__6 = (doublereal) (*al);
    d__7 = (doublereal) (*s);
    d__8 = (doublereal) (*be);
    d__9 = (doublereal) ((float)1. - *x);
    d__10 = (doublereal) (*d__);
    ret_val = (pow_dd(&d__1, &d__2) * (*ag + *bg * sx + *c__ * pow_dd(&d__3, &
	    d__4)) + pow_dd(&d__5, &d__6) * exp(-(*e) + sqrt(*es * pow_dd(&
	    d__7, &d__8) * lx))) * pow_dd(&d__9, &d__10);
/*<        RETURN >*/
    return ret_val;
/*<        END >*/
} /* fone_ */


/*<        FUNCTION FS(X, S, SF, AL, BE, AK, BK, AG, BG, C, D, E, ES) >*/
doublereal fs_(real *x, real *s, real *sf, real *al, real *be, real *ak, real 
	*bk, real *ag, real *bg, real *c__, real *d__, real *e, real *es)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8, d__9, d__10;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), pow_dd(doublereal *, doublereal 
	    *), exp(doublereal);

    /* Local variables */
    static real ds, lx, sx;

/*<        IMPLICIT REAL (A - Z) >*/
/*<        IF (S .LE. SF) THEN >*/
    if (*s <= *sf) {
/*<           FS = 0.0 >*/
	ret_val = (float)0.;
/*<        ELSE >*/
    } else {
/*<           SX = SQRT (X) >*/
	sx = sqrt(*x);
/*<           LX = ALOG (1./X) >*/
	lx = log((float)1. / *x);
/*<           DS = S - SF >*/
	ds = *s - *sf;
/*<        >*/
	d__1 = (doublereal) (*x);
	d__2 = (doublereal) (*ak);
	d__3 = (doublereal) (*x);
	d__4 = (doublereal) (*bk);
	d__5 = (doublereal) ds;
	d__6 = (doublereal) (*al);
	d__7 = (doublereal) (*s);
	d__8 = (doublereal) (*be);
	d__9 = (doublereal) ((float)1. - *x);
	d__10 = (doublereal) (*d__);
	ret_val = (ds * pow_dd(&d__1, &d__2) * (*ag + *bg * sx + *c__ * 
		pow_dd(&d__3, &d__4)) + pow_dd(&d__5, &d__6) * exp(-(*e) + 
		sqrt(*es * pow_dd(&d__7, &d__8) * lx))) * pow_dd(&d__9, &
		d__10);
/*<        END IF >*/
    }
/*<        RETURN >*/
    return ret_val;
/*<        END >*/
} /* fs_ */

#ifdef __cplusplus
	}
#endif
