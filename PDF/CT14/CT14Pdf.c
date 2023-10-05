/* CT14Pdf.f -- translated by f2c (version 20200916).
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
    integer nx, nt, nfmx, mxval;
} ctqpar2_;

#define ctqpar2_1 ctqpar2_

union {
    struct {
	doublereal alfaq, qalfa;
	integer ipk, iorder, nfl;
    } _1;
    struct {
	doublereal alfaq, qalfa;
	integer ipk, iorder, nfl0;
    } _2;
} qcdtbl_;

#define qcdtbl_1 (qcdtbl_._1)
#define qcdtbl_2 (qcdtbl_._2)

struct setchange_1_ {
    integer isetch, ipdsset, ipdsformat;
};

#define setchange_1 (*(struct setchange_1_ *) &setchange_)

union {
    struct {
	doublereal qini0, qmax0, xmin0;
    } _1;
    struct {
	doublereal qini, qmax, xmin;
    } _2;
} xqrange_;

#define xqrange_1 (xqrange_._1)
#define xqrange_2 (xqrange_._2)

struct {
    doublereal qbase, xv[202], tv[41], upd[88440], alscteq[41];
} ctqpar1_;

#define ctqpar1_1 ctqpar1_

struct {
    doublereal amass[6];
} masstbl_;

#define masstbl_1 masstbl_

/* Initialized data */

struct {
    integer fill_1[1];
    integer e_2[2];
    } setchange_ = { {0}, 0, 0 };


/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static integer c__5 = 5;
static integer c__3 = 3;

/* ============================================================================ */
/*                CTEQ-TEA Parton Distribution Functions: version 2014 */
/*                       December 31, 2014 */

/*   When using these PDFs, please cite the references below */


/*   This package provides a standard interface for CT10, CT12 */
/*   (unpublished), and CT14 parton distribution functions. */

/*   The following sets of CTEQ PDF table files can be computed */
/*    with this program: */
/*               PDF                             References */
/*   (1) 1+50 sets of CT10 NNLO PDF's;             [1] */
/*   (2' ) 1+52 sets of CT10 NLO PDF's;            [2] */
/*   (2'') 1+52 sets of CT10W NLO PDF's;           [2] */
/*   (3) 4 sets of CT10W NNLO and NLO PDF's        [2] */
/*       with alternative alpha_s values; */
/*   (4) 1+56 sets of CT14 NNLO PDF's;             [3] */
/*   (5) 1+56 sets of CT14 NLO PDF's;              [3] */
/*   (6) 2 sets of CT14 LO PDF's;                  [3] */
/*   (7) 11 CT14 NNLO sets and 11 CT14 NLO sets    [3] */
/*       with alternative alpha_s values */
/*   (8) 3 CT14 NNLO and 3 CT14 NLO sets           [3] */
/*       with up to 3, 4, and 6 active quark flavors */
/*   (9) 4 CT14 NNLO sets with intrinsic charm     [X] */
/*   References */
/*   [1] J. Gao, M. Guzzi, J. Huston, H.-L. Lai, Z. Li, P. M. Nadolsky, */
/*       J. Pumplin, D. Stump, C.-P. Yuan,  arXiv:1302.6246 [hep-ph] */
/*   [2] H.-L. Lai, M. Guzzi, J. Huston, Z. Li, P. M. Nadolsky, */
/*       J. Pumplin, and C.-P. Yuan, arXiv: 1007.2241 [hep-ph] */
/*   [3] S. Dulat, T.-J. Hou, J. Gao, M. Guzzi, J. Huston, */
/*       P. M. Nadolsky, J. Pumplin, C. Schmidt, D. Stump, and */
/*       C.-P. Yuan, arXiv:1506.XXXX */

/* =========================================================================== */
/*   The table grids are generated for */
/*    *  10^-9 < x < 1 and 1.3 < Q < 10^5 (GeV). */

/*   PDF values outside of the above range are returned using extrapolation. */

/*   The Table_Files are assumed to be in the working directory. */

/*   Before using the PDF, it is necessary to do the initialization by */
/*       Call SetCT14(Iset) */
/*   where Tablefile is a 40-character text string with the name of the */
/*   the desired PDF specified in the above table.table (.pds) file */

/*   Other provided functions include: */
/*   The function CT14Pdf (Iparton, X, Q) */
/*     returns the parton distribution inside the proton for parton [Iparton] */
/*     at [X] Bjorken_X and scale [Q] (GeV) in PDF set [Iset]. */
/*     Iparton  is the parton label (5, 4, 3, 2, 1, 0, -1, ......, -5) */
/*                              for (b, c, s, d, u, g, u_bar, ..., b_bar), */

/*   The function CT14Alphas (Q) */
/*     returns the value of the QCD coupling strength alpha_s(Q) at */
/*     an energy scale Q. The alpha_s value is obtained from the interpolation */
/*     of a numerical table of alpha_s included in each .pds file. */
/*     In the PDF fit, the alpha_s values are obtained at each momentum */
/*     scale by evolution in the HOPPET program at the respective QCD order */
/*     (NLO or NNLO). The table of alpha_s values at discrete Q values */
/*     is included in the input .pds file. The function CT14Alphas */
/*     estimates alpha_s at an arbitrary Q value, which agrees */
/*     with the direct evolution by HOPPET within a fraction of percent */
/*     point at typical Q. */

/*   The function CT14Mass(i) */
/*     returns the value of the quark mass for the i-th flavor. */
/*     The flavors are: */
/*     1  2  3  4  5  6 */
/*     u  d  s  c  b  t */

/*   Values of various PDF parameters assumed in the computation of the */
/*    PDFs can be obtained by */
/*     Call CT14GetPars( xmin,Qini,Qmax,Nloops,Nfl), */
/*   which returns */
/*     xmin, the minimal value of x; */
/*     Qmin,  the initial Q scale for the PDF evolution; */
/*     Qmax,  the maximal Q scale included in the PDF table; */
/*     Nloop, the number of QCD loops (order of the PDF in the QCD coupling); */
/*     Nfl,   the maximal number of quark flavors assumed in the PDF and */
/*            alpha_s evolution. */
/*   These programs, as provided, are in double precision.  By removing the */
/*   "Implicit Double Precision" lines, they can also be run in single */
/*   precision. */
/*   If you have detailed questions concerning these CT14 distributions, */
/*   or if you find problems/bugs using this package, direct inquires to */
/*   nadolsky@smu.edu. */

/* =========================================================================== */
/*<       Function CT14Pdf (Iparton, X, Q) >*/
doublereal ct14pdf_(integer *iparton, doublereal *x, doublereal *q)
{
    /* Initialized data */

    static logical warn = TRUE_;
    static doublereal qsml = .3;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();

    /* Local variables */
    extern doublereal partonx12_(integer *, doublereal *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___3 = { 0, 6, 0, 0, 0 };
    static cilist io___4 = { 0, 6, 0, 0, 0 };
    static cilist io___5 = { 0, 6, 0, 0, 0 };
    static cilist io___6 = { 0, 6, 0, 0, 0 };


/*<       Implicit Double Precision (A-H,O-Z) >*/
/*<       Logical Warn >*/
/*<       integer isetch, ipdsformat >*/
/*<        >*/
/*<       Data Warn /.true./ >*/
/*<       Data Qsml /.3d0/ >*/
/*<       save Warn >*/
/*<        >*/
    if (setchange_1.ipdsset != 1) {
	s_stop("CT14Pdf: the PDF table was not initialized", (ftnlen)42);
    }
/*<       If (X .lt. 0d0 .or. X .gt. 1D0) Then >*/
    if (*x < 0. || *x > 1.) {
/*<         Print *, 'X out of range in CT14Pdf: ', X >*/
	s_wsle(&io___3);
	do_lio(&c__9, &c__1, "X out of range in CT14Pdf: ", (ftnlen)27);
	do_lio(&c__5, &c__1, (char *)&(*x), (ftnlen)sizeof(doublereal));
	e_wsle();
/*<         CT14Pdf = 0D0 >*/
	ret_val = 0.;
/*<         Return >*/
	return ret_val;
/*<       Endif >*/
    }
/*<       If (Q .lt. Qsml) Then >*/
    if (*q < qsml) {
/*<         Print *, 'Q out of range in CT14Pdf: ', Q >*/
	s_wsle(&io___4);
	do_lio(&c__9, &c__1, "Q out of range in CT14Pdf: ", (ftnlen)27);
	do_lio(&c__5, &c__1, (char *)&(*q), (ftnlen)sizeof(doublereal));
	e_wsle();
/*<         Stop >*/
	s_stop("", (ftnlen)0);
/*<       Endif >*/
    }
/*<       If (abs(Iparton).gt. NfMx) Then >*/
    if (abs(*iparton) > ctqpar2_1.nfmx) {
/*<         If (Warn) Then >*/
	if (warn) {
/*        print a warning for calling extra flavor */
/*<           Warn = .false. >*/
	    warn = FALSE_;
/*<           Print *, 'Warning: Iparton out of range in CT14Pdf! ' >*/
	    s_wsle(&io___5);
	    do_lio(&c__9, &c__1, "Warning: Iparton out of range in CT14Pdf! ",
		     (ftnlen)42);
	    e_wsle();
/*<           Print *, 'Iparton, MxFlvN0: ', Iparton, NfMx >*/
	    s_wsle(&io___6);
	    do_lio(&c__9, &c__1, "Iparton, MxFlvN0: ", (ftnlen)18);
	    do_lio(&c__3, &c__1, (char *)&(*iparton), (ftnlen)sizeof(integer))
		    ;
	    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nfmx, (ftnlen)sizeof(
		    integer));
	    e_wsle();
/*<         Endif >*/
	}
/*<         CT14Pdf = 0D0 >*/
	ret_val = 0.;
/*<       else >*/
    } else {
/*<         CT14Pdf = PartonX12 (Iparton, X, Q) >*/
	ret_val = partonx12_(iparton, x, q);
/*<         if (CT14Pdf.lt.0D0) CT14Pdf = 0D0 >*/
	if (ret_val < 0.) {
	    ret_val = 0.;
	}
/*<       endif                     !if (abs(Iparton... >*/
    }
/*<       Return >*/
    return ret_val;
/*                             ******************** */
/*<       End >*/
} /* ct14pdf_ */

/*<       Subroutine SetCT14(Tablefile)     >*/
/* Subroutine */ int setct14_(char *tablefile, ftnlen tablefile_len)
{
    /* System generated locals */
    integer i__1;
    olist o__1;
    cllist cl__1;

    /* Builtin functions */
    integer f_open(olist *), f_clos(cllist *), s_wsle(cilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer iu;
    extern integer nextun_();
    extern /* Subroutine */ int readpds0_(integer *);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };


/*<       Implicit Double Precision (A-H,O-Z) >*/
/*<       Character Tablefile*40 >*/
/*<       Common /Setchange/ Isetch, ipdsset, ipdsformat >*/
/*<       data ipdsset, ipdsformat/0,0/ >*/
/*<       save >*/
/*<       IU= NextUn() >*/
    iu = nextun_();
/*<       Open(IU, File=Tablefile, Status='OLD', Err=100) >*/
    o__1.oerr = 1;
    o__1.ounit = iu;
    o__1.ofnmlen = 40;
    o__1.ofnm = tablefile;
    o__1.orl = 0;
    o__1.osta = "OLD";
    o__1.oacc = 0;
    o__1.ofm = 0;
    o__1.oblnk = 0;
    i__1 = f_open(&o__1);
    if (i__1 != 0) {
	goto L100;
    }
/*<       Call Readpds0 (IU) >*/
    readpds0_(&iu);
/*<       Close (IU) >*/
    cl__1.cerr = 0;
    cl__1.cunit = iu;
    cl__1.csta = 0;
    f_clos(&cl__1);
/*<       Isetch=1  >*/
    setchange_1.isetch = 1;
/*<       ipdsset=1 >*/
    setchange_1.ipdsset = 1;
/*<       Return >*/
    return 0;
/*<  100   >*/
L100:
    s_wsle(&io___8);
    do_lio(&c__9, &c__1, " Data file ", (ftnlen)11);
    do_lio(&c__9, &c__1, tablefile, (ftnlen)40);
    do_lio(&c__9, &c__1, " cannot be opened in SetCT14!!", (ftnlen)30);
    e_wsle();
/*<       Stop >*/
    s_stop("", (ftnlen)0);
/*                             ******************** */
/*<       End >*/
    return 0;
} /* setct14_ */

/*<       subroutine CT14GetPars(xmin,Qini,Qmax,Nloops,Nfl) >*/
/* Subroutine */ int ct14getpars_(doublereal *xmin, doublereal *qini, 
	doublereal *qmax, integer *nloops, integer *nfl)
{
/* Get various parameters associated with the PDF grid */
/* Output: xmin  is the minimal value of x */
/*         Qmin  is the initial Q scale */
/*         Qmax  is the maximal Q scale */
/*         Nloop is the number of QCD loops */
/*         Nfl   is the maximal number of quark flavors */
/*<       implicit none >*/
/*<       double precision Qini0, Qmax0, Xmin0, xmin, Qini, Qmax >*/
/*<       integer Nloops, Ipk, Iorder, Nfl,Nfl0 >*/
/*<       double precision AlfaQ, Qalfa >*/
/*<       common / XQrange / Qini0, Qmax0, Xmin0 >*/
/*<       common / QCDtbl /  AlfaQ, Qalfa, Ipk, Iorder, Nfl0 >*/
/*<       Qini=Qini0 >*/
    *qini = xqrange_1.qini0;
/*<       Qmax=Qmax0 >*/
    *qmax = xqrange_1.qmax0;
/*<       Xmin=Xmin0 >*/
    *xmin = xqrange_1.xmin0;
/*<       Nloops=Iorder-1 >*/
    *nloops = qcdtbl_2.iorder - 1;
/*<       Nfl=Nfl0 >*/
    *nfl = qcdtbl_2.nfl0;
/*<       return  >*/
    return 0;
/*<       end >*/
} /* ct14getpars_ */

/*<       Function CT14Alphas (QQ) >*/
doublereal ct14alphas_(doublereal *qq)
{
    /* Initialized data */

    static integer jq = 0;
    static doublereal q = -1.;

    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);
    integer s_wsle(cilist *), e_wsle(), do_lio(integer *, integer *, char *, 
	    ftnlen);
    double log(doublereal);

    /* Local variables */
    extern /* Subroutine */ int polint4f_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer jm, ju;
    static doublereal tt;
    static integer jlq;
    static doublereal alsout;

    /* Fortran I/O blocks */
    static cilist io___11 = { 0, 6, 0, 0, 0 };
    static cilist io___12 = { 0, 6, 0, 0, 0 };
    static cilist io___13 = { 0, 6, 0, 0, 0 };
    static cilist io___14 = { 0, 6, 0, 0, 0 };
    static cilist io___15 = { 0, 6, 0, 0, 0 };


/*<       Implicit Double Precision (A-H,O-Z) >*/
/*<       PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4) >*/
/*<       PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX) >*/
/*<       double precision Alsout >*/
/*<        >*/
/*<       Data Q, JQ /-1D0, 0/ >*/
/*<       save >*/
/*<        >*/
    if (setchange_1.ipdsset != 1) {
	s_stop("CT14Alphas: the PDF table was not initialized", (ftnlen)45);
    }
/*<       if (ipdsformat.lt.11) then >*/
    if (setchange_1.ipdsformat < 11) {
/*<         print * >*/
	s_wsle(&io___11);
	e_wsle();
/*<        >*/
	s_wsle(&io___12);
	do_lio(&c__9, &c__1, "STOP in CT14alphas: the PDF table file has an \
older format", (ftnlen)58);
	e_wsle();
/*<        >*/
	s_wsle(&io___13);
	do_lio(&c__9, &c__1, "and does not include the table of QCD coupling\
 values.", (ftnlen)54);
	e_wsle();
/*<        >*/
	s_wsle(&io___14);
	do_lio(&c__9, &c__1, "You can still compute the PDFs, but do not call"
		, (ftnlen)47);
	e_wsle();
/*<        >*/
	s_wsle(&io___15);
	do_lio(&c__9, &c__1, "the CT14alphas function for the interpolation \
of alpha_s.", (ftnlen)57);
	e_wsle();
/*<         stop >*/
	s_stop("", (ftnlen)0);
/*<       endif >*/
    }
/*<       Q = QQ >*/
    q = *qq;
/*<       tt = log(log(Q/qBase)) >*/
    tt = log(log(q / ctqpar1_1.qbase));
/*         --------------   Find lower end of interval containing Q, i.e., */
/*                          get jq such that qv(jq) .le. q .le. qv(jq+1)... */
/*<       JLq = -1 >*/
    jlq = -1;
/*<       JU = NT+1 >*/
    ju = ctqpar2_1.nt + 1;
/*<  13   If (JU-JLq .GT. 1) Then >*/
L13:
    if (ju - jlq > 1) {
/*<         JM = (JU+JLq) / 2 >*/
	jm = (ju + jlq) / 2;
/*<         If (tt .GE. TV(JM)) Then >*/
	if (tt >= ctqpar1_1.tv[jm]) {
/*<             JLq = JM >*/
	    jlq = jm;
/*<           Else >*/
	} else {
/*<             JU = JM >*/
	    ju = jm;
/*<           Endif >*/
	}
/*<           Goto 13 >*/
	goto L13;
/*<        Endif >*/
    }
/*<       If     (JLq .LE. 0) Then >*/
    if (jlq <= 0) {
/*<          Jq = 0 >*/
	jq = 0;
/*<       Elseif (JLq .LE. Nt-2) Then >*/
    } else if (jlq <= ctqpar2_1.nt - 2) {
/*                                  keep q in the middle, as shown above */
/*<          Jq = JLq - 1 >*/
	jq = jlq - 1;
/*<       Else >*/
    } else {
/*                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq. */
/*<         Jq = Nt - 3 >*/
	jq = ctqpar2_1.nt - 3;
/*<       Endif >*/
    }
/*                                 This is the interpolation variable in Q */
/*<       Call Polint4F (TV(jq), AlsCTEQ(jq), tt, Alsout) >*/
    polint4f_(&ctqpar1_1.tv[jq], &ctqpar1_1.alscteq[jq], &tt, &alsout);
/*<       CT14Alphas = Alsout >*/
    ret_val = alsout;
/*<       Return >*/
    return ret_val;
/*                                       ******************** */
/*<       End >*/
} /* ct14alphas_ */

/*<       function CT14Mass(i) >*/
doublereal ct14mass_(integer *i__)
{
    /* System generated locals */
    doublereal ret_val;

    /* Builtin functions */
    /* Subroutine */ int s_stop(char *, ftnlen);

/*     Returns the value of the quark mass for the i-th flavor */
/*     The flavors are: */
/*     1  2  3  4  5  6 */
/*     u  d  s  c  b  t */
/*<       implicit none >*/
/*<       double precision CT14Mass, Amass >*/
/*<       integer  Isetch, ipdsset, i, ipdsformat >*/
/*<        >*/
/*<        >*/
    if (setchange_1.ipdsset != 1) {
	s_stop("CT14Mass: the PDF table was not initialized", (ftnlen)43);
    }
/*<       CT14Mass = Amass(i) >*/
    ret_val = masstbl_1.amass[*i__ - 1];
/*<       return  >*/
    return ret_val;
/*<       end >*/
} /* ct14mass_ */

/*<       Subroutine Readpds0 (Nu) >*/
/* Subroutine */ int readpds0_(integer *nu)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_rsfe(cilist *), do_fio(integer *, char *, ftnlen), e_rsfe(), 
	    s_cmp(char *, char *, ftnlen, ftnlen), s_rsle(cilist *), do_lio(
	    integer *, integer *, char *, ftnlen), e_rsle(), i_dnnt(
	    doublereal *);
    double exp(doublereal);
    integer s_wsle(cilist *), e_wsle();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, n0;
    static doublereal aa, fl;
    static integer ng;
    static doublereal dr, qv[41];
    static integer nblk;
    static char line[80];
    static integer iret, npts;
    static doublereal dummy, qbase1, qbase2, aimass, alambda, fswitch;

    /* Fortran I/O blocks */
    static cilist io___21 = { 0, 0, 0, "(A)", 0 };
    static cilist io___23 = { 0, 0, 0, "(A)", 0 };
    static cilist io___24 = { 0, 0, 0, 0, 0 };
    static cilist io___27 = { 0, 0, 0, "(A)", 0 };
    static cilist io___28 = { 0, 0, 0, 0, 0 };
    static cilist io___32 = { 0, 0, 0, 0, 0 };
    static cilist io___33 = { 0, 0, 0, 0, 0 };
    static cilist io___36 = { 0, 0, 0, "(A)", 0 };
    static cilist io___37 = { 0, 0, 0, 0, 0 };
    static cilist io___39 = { 0, 0, 0, "(A)", 0 };
    static cilist io___40 = { 0, 0, 0, 0, 0 };
    static cilist io___42 = { 0, 0, 0, "(A)", 0 };
    static cilist io___43 = { 0, 0, 0, "(A)", 0 };
    static cilist io___44 = { 0, 0, 0, 0, 0 };
    static cilist io___46 = { 0, 0, 0, 0, 0 };
    static cilist io___49 = { 0, 6, 0, 0, 0 };
    static cilist io___50 = { 0, 6, 0, 0, 0 };
    static cilist io___51 = { 0, 0, 0, "(A)", 0 };
    static cilist io___52 = { 0, 0, 0, 0, 0 };
    static cilist io___56 = { 0, 0, 0, "(A)", 0 };
    static cilist io___58 = { 1, 0, 1, 0, 0 };


/*<       Implicit Double Precision (A-H,O-Z) >*/
/*<       Character Line*80 >*/
/*<       integer ipdsformat >*/
/*<       PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4) >*/
/*<       PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX) >*/
/*<       double precision qv(0:mxq) >*/
/*<        >*/
/*<       Read  (Nu, '(A)') Line >*/
    io___21.ciunit = *nu;
    s_rsfe(&io___21);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       Read  (Nu, '(A)') Line >*/
    io___23.ciunit = *nu;
    s_rsfe(&io___23);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       if (Line(1:11) .eq. '  ipk, Ordr') then !post-CT10 .pds format; >*/
    if (s_cmp(line, "  ipk, Ordr", (ftnlen)11, (ftnlen)11) == 0) {
/* Set alphas(MZ) at scale Zm, quark masses, and evolution type */
/*<         ipdsformat = 10           !Post-CT10 .pds format >*/
	setchange_1.ipdsformat = 10;
/*<         Read (Nu, *) ipk, Dr, Qalfa, AlfaQ, (amass(i),i=1,6)  >*/
	io___24.ciunit = *nu;
	s_rsle(&io___24);
	do_lio(&c__3, &c__1, (char *)&qcdtbl_1.ipk, (ftnlen)sizeof(integer));
	do_lio(&c__5, &c__1, (char *)&dr, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&qcdtbl_1.qalfa, (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&qcdtbl_1.alfaq, (ftnlen)sizeof(
		doublereal));
	for (i__ = 1; i__ <= 6; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&masstbl_1.amass[i__ - 1], (ftnlen)
		    sizeof(doublereal));
	}
	e_rsle();
/*<         Iorder = Nint(Dr)         >*/
	qcdtbl_1.iorder = i_dnnt(&dr);
/*<         read (Nu, '(A)') Line >*/
	io___27.ciunit = *nu;
	s_rsfe(&io___27);
	do_fio(&c__1, line, (ftnlen)80);
	e_rsfe();
/*<         if (Line(1:7) .eq. '  IMASS' ) then >*/
	if (s_cmp(line, "  IMASS", (ftnlen)7, (ftnlen)7) == 0) {
/*<           ipdsformat = 11         !CT12 .pds format >*/
	    setchange_1.ipdsformat = 11;
/*<           read (Nu, *) aimass, fswitch, N0, N0, N0, Nfmx, MxVal >*/
	    io___28.ciunit = *nu;
	    s_rsle(&io___28);
	    do_lio(&c__5, &c__1, (char *)&aimass, (ftnlen)sizeof(doublereal));
	    do_lio(&c__5, &c__1, (char *)&fswitch, (ftnlen)sizeof(doublereal))
		    ;
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nfmx, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.mxval, (ftnlen)sizeof(
		    integer));
	    e_rsle();
/*<           Nfl=Nfmx >*/
	    qcdtbl_1.nfl = ctqpar2_1.nfmx;
/*<         else                      !Pre-CT12 format >*/
	} else {
/*<           Read  (Nu, *) N0, N0, N0, NfMx, MxVal >*/
	    io___32.ciunit = *nu;
	    s_rsle(&io___32);
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nfmx, (ftnlen)sizeof(
		    integer));
	    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.mxval, (ftnlen)sizeof(
		    integer));
	    e_rsle();
/*<         endif                     !Line(1:7) >*/
	}
/*<       else                        !old .pds format;       >*/
    } else {
/*<         ipdsformat = 6            !CTEQ6.6 .pds format; alpha_s  is not  >*/
	setchange_1.ipdsformat = 6;
/*<         Read (Nu, *) Dr, fl, Alambda, (amass(i),i=1,6)  !set Lambda_QCD >*/
	io___33.ciunit = *nu;
	s_rsle(&io___33);
	do_lio(&c__5, &c__1, (char *)&dr, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&fl, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&alambda, (ftnlen)sizeof(doublereal));
	for (i__ = 1; i__ <= 6; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&masstbl_1.amass[i__ - 1], (ftnlen)
		    sizeof(doublereal));
	}
	e_rsle();
/*<         Iorder = Nint(Dr) >*/
	qcdtbl_1.iorder = i_dnnt(&dr);
/*<         Nfl = Nint(fl) >*/
	qcdtbl_1.nfl = i_dnnt(&fl);
/*<         Read  (Nu, '(A)') Line >*/
	io___36.ciunit = *nu;
	s_rsfe(&io___36);
	do_fio(&c__1, line, (ftnlen)80);
	e_rsfe();
/*<         Read  (Nu, *) dummy,dummy,dummy, NfMx, MxVal, N0 >*/
	io___37.ciunit = *nu;
	s_rsle(&io___37);
	do_lio(&c__5, &c__1, (char *)&dummy, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&dummy, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&dummy, (ftnlen)sizeof(doublereal));
	do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nfmx, (ftnlen)sizeof(integer))
		;
	do_lio(&c__3, &c__1, (char *)&ctqpar2_1.mxval, (ftnlen)sizeof(integer)
		);
	do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
	e_rsle();
/*<       endif                       !Line(1:11... >*/
    }
/*<       Read  (Nu, '(A)') Line >*/
    io___39.ciunit = *nu;
    s_rsfe(&io___39);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       Read  (Nu, *) NX,  NT, N0, NG, N0 >*/
    io___40.ciunit = *nu;
    s_rsle(&io___40);
    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nx, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ctqpar2_1.nt, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&ng, (ftnlen)sizeof(integer));
    do_lio(&c__3, &c__1, (char *)&n0, (ftnlen)sizeof(integer));
    e_rsle();
/*<       if (ng.gt.0) Read  (Nu, '(A)') (Line, i=1,ng+1) >*/
    if (ng > 0) {
	io___42.ciunit = *nu;
	s_rsfe(&io___42);
	i__1 = ng + 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    do_fio(&c__1, line, (ftnlen)80);
	}
	e_rsfe();
    }
/*<       Read  (Nu, '(A)') Line >*/
    io___43.ciunit = *nu;
    s_rsfe(&io___43);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       if (ipdsformat.ge.11) then  !CT12 format with alpha_s values >*/
    if (setchange_1.ipdsformat >= 11) {
/*<         Read  (Nu, *) QINI, QMAX, (qv(i),TV(I), AlsCTEQ(I), I =0, NT) >*/
	io___44.ciunit = *nu;
	s_rsle(&io___44);
	do_lio(&c__5, &c__1, (char *)&xqrange_2.qini, (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&xqrange_2.qmax, (ftnlen)sizeof(
		doublereal));
	i__1 = ctqpar2_1.nt;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&qv[i__], (ftnlen)sizeof(doublereal))
		    ;
	    do_lio(&c__5, &c__1, (char *)&ctqpar1_1.tv[i__], (ftnlen)sizeof(
		    doublereal));
	    do_lio(&c__5, &c__1, (char *)&ctqpar1_1.alscteq[i__], (ftnlen)
		    sizeof(doublereal));
	}
	e_rsle();
/*<       else                        !pre-CT12 format >*/
    } else {
/*<         Read  (Nu, *) QINI, QMAX, (qv(i),TV(I), I =0, NT) >*/
	io___46.ciunit = *nu;
	s_rsle(&io___46);
	do_lio(&c__5, &c__1, (char *)&xqrange_2.qini, (ftnlen)sizeof(
		doublereal));
	do_lio(&c__5, &c__1, (char *)&xqrange_2.qmax, (ftnlen)sizeof(
		doublereal));
	i__1 = ctqpar2_1.nt;
	for (i__ = 0; i__ <= i__1; ++i__) {
	    do_lio(&c__5, &c__1, (char *)&qv[i__], (ftnlen)sizeof(doublereal))
		    ;
	    do_lio(&c__5, &c__1, (char *)&ctqpar1_1.tv[i__], (ftnlen)sizeof(
		    doublereal));
	}
	e_rsle();
/*<       endif                       !ipdsformat.ge.11 >*/
    }
/* check that qBase is consistent with the definition of Tv(0:nQ) for 2 values of Qv */
/*<       qbase1 = Qv(1)/Exp(Exp(Tv(1))) >*/
    qbase1 = qv[1] / exp(exp(ctqpar1_1.tv[1]));
/*<       qbase2 = Qv(nT)/Exp(Exp(Tv(NT))) >*/
    qbase2 = qv[ctqpar2_1.nt] / exp(exp(ctqpar1_1.tv[ctqpar2_1.nt]));
/*<       if (abs(qbase1-qbase2).gt.1e-5) then >*/
    if ((d__1 = qbase1 - qbase2, abs(d__1)) > (float)1e-5) {
/*<         print *, 'Readpds0: something wrong with qbase' >*/
	s_wsle(&io___49);
	do_lio(&c__9, &c__1, "Readpds0: something wrong with qbase", (ftnlen)
		36);
	e_wsle();
/*<         print *,'qbase1, qbase2=',qbase1,qbase2 >*/
	s_wsle(&io___50);
	do_lio(&c__9, &c__1, "qbase1, qbase2=", (ftnlen)15);
	do_lio(&c__5, &c__1, (char *)&qbase1, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&qbase2, (ftnlen)sizeof(doublereal));
	e_wsle();
/*<         stop >*/
	s_stop("", (ftnlen)0);
/*<       else >*/
    } else {
/*<         qbase=(qbase1+qbase2)/2.0d0 >*/
	ctqpar1_1.qbase = (qbase1 + qbase2) / 2.;
/*<       endif                     !abs(qbase1... >*/
    }
/*<       Read  (Nu, '(A)') Line >*/
    io___51.ciunit = *nu;
    s_rsfe(&io___51);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       Read  (Nu, *) XMIN, aa, (XV(I), I =1, NX) >*/
    io___52.ciunit = *nu;
    s_rsle(&io___52);
    do_lio(&c__5, &c__1, (char *)&xqrange_2.xmin, (ftnlen)sizeof(doublereal));
    do_lio(&c__5, &c__1, (char *)&aa, (ftnlen)sizeof(doublereal));
    i__1 = ctqpar2_1.nx;
    for (i__ = 1; i__ <= i__1; ++i__) {
	do_lio(&c__5, &c__1, (char *)&ctqpar1_1.xv[i__], (ftnlen)sizeof(
		doublereal));
    }
    e_rsle();
/*<       XV(0)=0D0 >*/
    ctqpar1_1.xv[0] = 0.;
/*<       Nblk = (NX+1) * (NT+1) >*/
    nblk = (ctqpar2_1.nx + 1) * (ctqpar2_1.nt + 1);
/*<       Npts =  Nblk  * (NfMx+1+MxVal) >*/
    npts = nblk * (ctqpar2_1.nfmx + 1 + ctqpar2_1.mxval);
/*<       Read  (Nu, '(A)') Line >*/
    io___56.ciunit = *nu;
    s_rsfe(&io___56);
    do_fio(&c__1, line, (ftnlen)80);
    e_rsfe();
/*<       Read  (Nu, *, IOSTAT=IRET) (UPD(I), I=1,Npts) >*/
    io___58.ciunit = *nu;
    iret = s_rsle(&io___58);
    if (iret != 0) {
	goto L100001;
    }
    i__1 = npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	iret = do_lio(&c__5, &c__1, (char *)&ctqpar1_1.upd[i__ - 1], (ftnlen)
		sizeof(doublereal));
	if (iret != 0) {
	    goto L100001;
	}
    }
    iret = e_rsle();
L100001:
/*<       Return >*/
    return 0;
/*                        **************************** */
/*<       End >*/
} /* readpds0_ */

/*<       Function PartonX12 (IPRTN, XX, QQ) >*/
doublereal partonx12_(integer *iprtn, doublereal *xx, doublereal *qq)
{
    /* Initialized data */

    static doublereal onep = 1.00001;
    static doublereal xpow = .3;
    static integer nqvec = 4;
    static doublereal x = -1.;
    static doublereal q = -1.;
    static integer jx = 0;
    static integer jq = 0;

    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe();
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern /* Subroutine */ int polint4f_(doublereal *, doublereal *, 
	    doublereal *, doublereal *);
    static integer i__;
    static doublereal g1;
    static integer j1;
    static doublereal g4, h00, ff, s12, t12, s23, t13, t23, t24, t34;
    static integer jm, ju;
    static doublereal s13, s24, s34, ss, tt;
    static integer ip, it;
    static doublereal fx, sf2, sf3, tf2, tf3, sy2, sy3, ty2, ty3, fij[4], 
	    s1213, s2434;
    static integer jlq, jlx;
    static doublereal tmp, tmp1, tmp2, fvec[4], sdet, tdet;
    static integer jtmp;
    static doublereal svec1, svec2, svec3, svec4, tvec1, tvec2, tvec3, tvec4, 
	    xvpow[202], const1, const2, const3, const4, const5, const6;

    /* Fortran I/O blocks */
    static cilist io___73 = { 0, 6, 0, "(A,1pE12.4)", 0 };
    static cilist io___74 = { 0, 6, 0, "(A,1pE12.4)", 0 };


/*  Given the parton distribution function in the array U in */
/*  COMMON / PEVLDT / , this routine interpolates to find */
/*  the parton distribution at an arbitray point in x and q. */

/*<       Implicit Double Precision (A-H,O-Z) >*/
/*<       PARAMETER (MXX = 201, MXQ = 40, MXF = 6, MaxVal=4) >*/
/*<       PARAMETER (MXPQX = (MXF+1+MaxVal) * MXQ * MXX) >*/
/*<        >*/
/*<       Dimension fvec(4), fij(4) >*/
/*<       Dimension xvpow(0:mxx) >*/
/*<       Data OneP / 1.00001 / >*/
/*<       Data xpow / 0.3d0 /       !**** choice of interpolation variable >*/
/*<       Data nqvec / 4 / >*/
/*<       Data ientry / 0 / >*/
/*<       Data X, Q, JX, JQ /-1D0, -1D0, 0, 0/ >*/
/*<       Save xvpow >*/
/*<       Save X, Q, JX, JQ, JLX, JLQ >*/
/*<       Save ss, const1, const2, const3, const4, const5, const6 >*/
/*<       Save sy2, sy3, s23, tt, t12, t13, t23, t24, t34, ty2, ty3 >*/
/*<       Save tmp1, tmp2, tdet >*/
/* store the powers used for interpolation on first call... */
/*<       if(Isetch .eq. 1) then >*/
    if (setchange_1.isetch == 1) {
/*<          Isetch = 0 >*/
	setchange_1.isetch = 0;
/*<          xvpow(0) = 0D0 >*/
	xvpow[0] = 0.;
/*<          do i = 1, nx >*/
	i__1 = ctqpar2_1.nx;
	for (i__ = 1; i__ <= i__1; ++i__) {
/*<             xvpow(i) = xv(i)**xpow >*/
	    xvpow[i__] = pow_dd(&ctqpar1_1.xv[i__], &xpow);
/*<          enddo >*/
	}
/*<       elseIf((XX.eq.X).and.(QQ.eq.Q)) then >*/
    } else if (*xx == x && *qq == q) {
/*<       	goto 99 >*/
	goto L99;
/*<       endif >*/
    }
/*<       X = XX >*/
    x = *xx;
/*<       Q = QQ >*/
    q = *qq;
/*<       tt = log(log(Q/qBase)) >*/
    tt = log(log(q / ctqpar1_1.qbase));
/*      -------------    find lower end of interval containing x, i.e., */
/*                       get jx such that xv(jx) .le. x .le. xv(jx+1)... */
/*<       JLx = -1 >*/
    jlx = -1;
/*<       JU = Nx+1 >*/
    ju = ctqpar2_1.nx + 1;
/*<  11   If (JU-JLx .GT. 1) Then >*/
L11:
    if (ju - jlx > 1) {
/*<          JM = (JU+JLx) / 2 >*/
	jm = (ju + jlx) / 2;
/*<          If (X .Ge. XV(JM)) Then >*/
	if (x >= ctqpar1_1.xv[jm]) {
/*<             JLx = JM >*/
	    jlx = jm;
/*<          Else >*/
	} else {
/*<             JU = JM >*/
	    ju = jm;
/*<          Endif >*/
	}
/*<          Goto 11 >*/
	goto L11;
/*<       Endif >*/
    }
/*                     Ix    0   1   2      Jx  JLx         Nx-2     Nx */
/*                           |---|---|---|...|---|-x-|---|...|---|---| */
/*                     x     0  Xmin               x                 1 */

/*<       If     (JLx .LE. -1) Then >*/
    if (jlx <= -1) {
/*<         Print '(A,1pE12.4)','Severe error: x <= 0 in PartonX12! x = ',x >*/
	s_wsfe(&io___73);
	do_fio(&c__1, "Severe error: x <= 0 in PartonX12! x = ", (ftnlen)39);
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(doublereal));
	e_wsfe();
/*<         Stop >*/
	s_stop("", (ftnlen)0);
/*<       ElseIf (JLx .Eq. 0) Then >*/
    } else if (jlx == 0) {
/*<          Jx = 0 >*/
	jx = 0;
/*<       Elseif (JLx .LE. Nx-2) Then >*/
    } else if (jlx <= ctqpar2_1.nx - 2) {
/*                For interrior points, keep x in the middle, as shown above */
/*<          Jx = JLx - 1 >*/
	jx = jlx - 1;
/*<       Elseif (JLx.Eq.Nx-1 .or. x.LT.OneP) Then >*/
    } else if (jlx == ctqpar2_1.nx - 1 || x < onep) {
/*                  We tolerate a slight over-shoot of one (OneP=1.00001), */
/*              perhaps due to roundoff or whatever, but not more than that. */
/*                                      Keep at least 4 points >= Jx */
/*<          Jx = JLx - 2 >*/
	jx = jlx - 2;
/*<       Else >*/
    } else {
/*<         Print '(A,1pE12.4)','Severe error: x > 1 in PartonX12! x = ',x >*/
	s_wsfe(&io___74);
	do_fio(&c__1, "Severe error: x > 1 in PartonX12! x = ", (ftnlen)38);
	do_fio(&c__1, (char *)&x, (ftnlen)sizeof(doublereal));
	e_wsfe();
/*<         Stop >*/
	s_stop("", (ftnlen)0);
/*<       Endif >*/
    }
/*          ---------- Note: JLx uniquely identifies the x-bin; Jx does not. */
/*                       This is the variable to be interpolated in */
/*<       ss = x**xpow >*/
    ss = pow_dd(&x, &xpow);
/*<       If (JLx.Ge.2 .and. JLx.Le.Nx-2) Then >*/
    if (jlx >= 2 && jlx <= ctqpar2_1.nx - 2) {
/*     initiation work for "interior bins": store the lattice points in s... */
/*<       svec1 = xvpow(jx) >*/
	svec1 = xvpow[jx];
/*<       svec2 = xvpow(jx+1) >*/
	svec2 = xvpow[jx + 1];
/*<       svec3 = xvpow(jx+2) >*/
	svec3 = xvpow[jx + 2];
/*<       svec4 = xvpow(jx+3) >*/
	svec4 = xvpow[jx + 3];
/*<       s12 = svec1 - svec2 >*/
	s12 = svec1 - svec2;
/*<       s13 = svec1 - svec3 >*/
	s13 = svec1 - svec3;
/*<       s23 = svec2 - svec3 >*/
	s23 = svec2 - svec3;
/*<       s24 = svec2 - svec4 >*/
	s24 = svec2 - svec4;
/*<       s34 = svec3 - svec4 >*/
	s34 = svec3 - svec4;
/*<       sy2 = ss - svec2 >*/
	sy2 = ss - svec2;
/*<       sy3 = ss - svec3 >*/
	sy3 = ss - svec3;
/* constants needed for interpolating in s at fixed t lattice points... */
/*<       const1 = s13/s23 >*/
	const1 = s13 / s23;
/*<       const2 = s12/s23 >*/
	const2 = s12 / s23;
/*<       const3 = s34/s23 >*/
	const3 = s34 / s23;
/*<       const4 = s24/s23 >*/
	const4 = s24 / s23;
/*<       s1213 = s12 + s13 >*/
	s1213 = s12 + s13;
/*<       s2434 = s24 + s34 >*/
	s2434 = s24 + s34;
/*<       sdet = s12*s34 - s1213*s2434 >*/
	sdet = s12 * s34 - s1213 * s2434;
/*<       tmp = sy2*sy3/sdet >*/
	tmp = sy2 * sy3 / sdet;
/*<       const5 = (s34*sy2-s2434*sy3)*tmp/s12 >*/
	const5 = (s34 * sy2 - s2434 * sy3) * tmp / s12;
/*<       const6 = (s1213*sy2-s12*sy3)*tmp/s34 >*/
	const6 = (s1213 * sy2 - s12 * sy3) * tmp / s34;
/*<       EndIf >*/
    }
/*         --------------Now find lower end of interval containing Q, i.e., */
/*                          get jq such that qv(jq) .le. q .le. qv(jq+1)... */
/*<       JLq = -1 >*/
    jlq = -1;
/*<       JU = NT+1 >*/
    ju = ctqpar2_1.nt + 1;
/*<  12   If (JU-JLq .GT. 1) Then >*/
L12:
    if (ju - jlq > 1) {
/*<          JM = (JU+JLq) / 2 >*/
	jm = (ju + jlq) / 2;
/*<          If (tt .GE. TV(JM)) Then >*/
	if (tt >= ctqpar1_1.tv[jm]) {
/*<             JLq = JM >*/
	    jlq = jm;
/*<          Else >*/
	} else {
/*<             JU = JM >*/
	    ju = jm;
/*<          Endif >*/
	}
/*<          Goto 12 >*/
	goto L12;
/*<        Endif >*/
    }
/*<       If     (JLq .LE. 0) Then >*/
    if (jlq <= 0) {
/*<          Jq = 0 >*/
	jq = 0;
/*<       Elseif (JLq .LE. Nt-2) Then >*/
    } else if (jlq <= ctqpar2_1.nt - 2) {
/*                                  keep q in the middle, as shown above */
/*<          Jq = JLq - 1 >*/
	jq = jlq - 1;
/*<       Else >*/
    } else {
/*                         JLq .GE. Nt-1 case:  Keep at least 4 points >= Jq. */
/*<         Jq = Nt - 3 >*/
	jq = ctqpar2_1.nt - 3;
/*<       Endif >*/
    }
/*                                   This is the interpolation variable in Q */
/*<       If (JLq.GE.1 .and. JLq.LE.Nt-2) Then >*/
    if (jlq >= 1 && jlq <= ctqpar2_1.nt - 2) {
/*                                        store the lattice points in t... */
/*<       tvec1 = Tv(jq) >*/
	tvec1 = ctqpar1_1.tv[jq];
/*<       tvec2 = Tv(jq+1) >*/
	tvec2 = ctqpar1_1.tv[jq + 1];
/*<       tvec3 = Tv(jq+2) >*/
	tvec3 = ctqpar1_1.tv[jq + 2];
/*<       tvec4 = Tv(jq+3) >*/
	tvec4 = ctqpar1_1.tv[jq + 3];
/*<       t12 = tvec1 - tvec2 >*/
	t12 = tvec1 - tvec2;
/*<       t13 = tvec1 - tvec3 >*/
	t13 = tvec1 - tvec3;
/*<       t23 = tvec2 - tvec3 >*/
	t23 = tvec2 - tvec3;
/*<       t24 = tvec2 - tvec4 >*/
	t24 = tvec2 - tvec4;
/*<       t34 = tvec3 - tvec4 >*/
	t34 = tvec3 - tvec4;
/*<       ty2 = tt - tvec2 >*/
	ty2 = tt - tvec2;
/*<       ty3 = tt - tvec3 >*/
	ty3 = tt - tvec3;
/*<       tmp1 = t12 + t13 >*/
	tmp1 = t12 + t13;
/*<       tmp2 = t24 + t34 >*/
	tmp2 = t24 + t34;
/*<       tdet = t12*t34 - tmp1*tmp2 >*/
	tdet = t12 * t34 - tmp1 * tmp2;
/*<       EndIf >*/
    }
/* get the pdf function values at the lattice points... */
/*<  99   If (Iprtn .Gt. MxVal) Then >*/
L99:
    if (*iprtn > ctqpar2_1.mxval) {
/*<          Ip = - Iprtn >*/
	ip = -(*iprtn);
/*<       Else >*/
    } else {
/*<          Ip = Iprtn >*/
	ip = *iprtn;
/*<       EndIf >*/
    }
/*<       jtmp = ((Ip + NfMx)*(NT+1)+(jq-1))*(NX+1)+jx+1 >*/
    jtmp = ((ip + ctqpar2_1.nfmx) * (ctqpar2_1.nt + 1) + (jq - 1)) * (
	    ctqpar2_1.nx + 1) + jx + 1;
/*<       Do it = 1, nqvec >*/
    i__1 = nqvec;
    for (it = 1; it <= i__1; ++it) {
/*<         J1  = jtmp + it*(NX+1) >*/
	j1 = jtmp + it * (ctqpar2_1.nx + 1);
/*<        If (Jx .Eq. 0) Then >*/
	if (jx == 0) {
/*                      For the first 4 x points, interpolate x^2*f(x,Q) */
/*                      This applies to the two lowest bins JLx = 0, 1 */
/*            We can not put the JLx.eq.1 bin into the "interrior" section */
/*                           (as we do for q), since Upd(J1) is undefined. */
/*<          fij(1) = 0 >*/
	    fij[0] = 0.;
/*<          fij(2) = Upd(J1+1) * XV(1)**2 >*/
/* Computing 2nd power */
	    d__1 = ctqpar1_1.xv[1];
	    fij[1] = ctqpar1_1.upd[j1] * (d__1 * d__1);
/*<          fij(3) = Upd(J1+2) * XV(2)**2 >*/
/* Computing 2nd power */
	    d__1 = ctqpar1_1.xv[2];
	    fij[2] = ctqpar1_1.upd[j1 + 1] * (d__1 * d__1);
/*<          fij(4) = Upd(J1+3) * XV(3)**2 >*/
/* Computing 2nd power */
	    d__1 = ctqpar1_1.xv[3];
	    fij[3] = ctqpar1_1.upd[j1 + 2] * (d__1 * d__1);

/*                 Use Polint which allows x to be anywhere w.r.t. the grid */
/*<          Call Polint4F (XVpow(0), Fij(1), ss, Fx) >*/
	    polint4f_(xvpow, fij, &ss, &fx);
/*<          If (x .GT. 0D0)  Fvec(it) =  Fx / x**2 >*/
	    if (x > 0.) {
/* Computing 2nd power */
		d__1 = x;
		fvec[it - 1] = fx / (d__1 * d__1);
	    }
/*                                              Pdf is undefined for x.eq.0 */
/*<        ElseIf  (JLx .Eq. Nx-1) Then >*/
	} else if (jlx == ctqpar2_1.nx - 1) {
/*                                                This is the highest x bin: */
/*<         Call Polint4F (XVpow(Nx-3), Upd(J1), ss, Fx) >*/
	    polint4f_(&xvpow[ctqpar2_1.nx - 3], &ctqpar1_1.upd[j1 - 1], &ss, &
		    fx);
/*<         Fvec(it) = Fx >*/
	    fvec[it - 1] = fx;
/*<        Else >*/
	} else {
/*                       for all interior points, use Jon's in-line function */
/*                              This applied to (JLx.Ge.2 .and. JLx.Le.Nx-2) */
/*<          sf2 = Upd(J1+1) >*/
	    sf2 = ctqpar1_1.upd[j1];
/*<          sf3 = Upd(J1+2) >*/
	    sf3 = ctqpar1_1.upd[j1 + 1];
/*<          g1 =  sf2*const1 - sf3*const2 >*/
	    g1 = sf2 * const1 - sf3 * const2;
/*<          g4 = -sf2*const3 + sf3*const4 >*/
	    g4 = -sf2 * const3 + sf3 * const4;
/*<        >*/
	    fvec[it - 1] = (const5 * (ctqpar1_1.upd[j1 - 1] - g1) + const6 * (
		    ctqpar1_1.upd[j1 + 2] - g4) + sf2 * sy3 - sf3 * sy2) / 
		    s23;
/*<        Endif >*/
	}
/*<       enddo >*/
    }
/*                                   We now have the four values Fvec(1:4) */
/*     interpolate in t... */
/*<       If (JLq .LE. 0) Then >*/
    if (jlq <= 0) {
/*                         1st Q-bin, as well as extrapolation to lower Q */
/*<         Call Polint4F (TV(0), Fvec(1), tt, ff) >*/
	polint4f_(ctqpar1_1.tv, fvec, &tt, &ff);
/*<       ElseIf (JLq .GE. Nt-1) Then >*/
    } else if (jlq >= ctqpar2_1.nt - 1) {
/*                         Last Q-bin, as well as extrapolation to higher Q */
/*<         Call Polint4F (TV(Nt-3), Fvec(1), tt, ff) >*/
	polint4f_(&ctqpar1_1.tv[ctqpar2_1.nt - 3], fvec, &tt, &ff);
/*<       Else >*/
    } else {
/*                         Interrior bins : (JLq.GE.1 .and. JLq.LE.Nt-2) */
/*       which include JLq.Eq.1 and JLq.Eq.Nt-2, since Upd is defined for */
/*                         the full range QV(0:Nt)  (in contrast to XV) */
/*<         tf2 = fvec(2) >*/
	tf2 = fvec[1];
/*<         tf3 = fvec(3) >*/
	tf3 = fvec[2];
/*<         g1 = ( tf2*t13 - tf3*t12) / t23 >*/
	g1 = (tf2 * t13 - tf3 * t12) / t23;
/*<         g4 = (-tf2*t34 + tf3*t24) / t23 >*/
	g4 = (-tf2 * t34 + tf3 * t24) / t23;
/*<        >*/
	h00 = (t34 * ty2 - tmp2 * ty3) * (fvec[0] - g1) / t12 + (tmp1 * ty2 - 
		t12 * ty3) * (fvec[3] - g4) / t34;
/*<         ff = (h00*ty2*ty3/tdet + tf2*ty3 - tf3*ty2) / t23 >*/
	ff = (h00 * ty2 * ty3 / tdet + tf2 * ty3 - tf3 * ty2) / t23;
/*<       EndIf >*/
    }
/*<       PartonX12 = ff >*/
    ret_val = ff;
/*<       Return >*/
    return ret_val;
/*                                       ******************** */
/*<       End >*/
} /* partonx12_ */

/*<       SUBROUTINE POLINT4F (XA,YA,X,Y) >*/
/* Subroutine */ int polint4f_(doublereal *xa, doublereal *ya, doublereal *x, 
	doublereal *y)
{
    static doublereal w, c1, d1, d2, c2, d3, h1, h2, h3, h4, c3, cc1, cd1, 
	    cd2, cc2, dd1, dc1, den;

/*<       IMPLICIT DOUBLE PRECISION (A-H, O-Z) >*/
/*  The POLINT4 routine is based on the POLINT routine from "Numerical Recipes", */
/*  but assuming N=4, and ignoring the error estimation. */
/*  suggested by Z. Sullivan. */
/*<       DIMENSION XA(*),YA(*) >*/
/*<       H1=XA(1)-X >*/
    /* Parameter adjustments */
    --ya;
    --xa;

    /* Function Body */
    h1 = xa[1] - *x;
/*<       H2=XA(2)-X >*/
    h2 = xa[2] - *x;
/*<       H3=XA(3)-X >*/
    h3 = xa[3] - *x;
/*<       H4=XA(4)-X >*/
    h4 = xa[4] - *x;
/*<       W=YA(2)-YA(1) >*/
    w = ya[2] - ya[1];
/*<       DEN=W/(H1-H2) >*/
    den = w / (h1 - h2);
/*<       D1=H2*DEN >*/
    d1 = h2 * den;
/*<       C1=H1*DEN >*/
    c1 = h1 * den;
/*<       W=YA(3)-YA(2) >*/
    w = ya[3] - ya[2];
/*<       DEN=W/(H2-H3) >*/
    den = w / (h2 - h3);
/*<       D2=H3*DEN >*/
    d2 = h3 * den;
/*<       C2=H2*DEN >*/
    c2 = h2 * den;
/*<       W=YA(4)-YA(3) >*/
    w = ya[4] - ya[3];
/*<       DEN=W/(H3-H4) >*/
    den = w / (h3 - h4);
/*<       D3=H4*DEN >*/
    d3 = h4 * den;
/*<       C3=H3*DEN >*/
    c3 = h3 * den;
/*<       W=C2-D1 >*/
    w = c2 - d1;
/*<       DEN=W/(H1-H3) >*/
    den = w / (h1 - h3);
/*<       CD1=H3*DEN >*/
    cd1 = h3 * den;
/*<       CC1=H1*DEN >*/
    cc1 = h1 * den;
/*<       W=C3-D2 >*/
    w = c3 - d2;
/*<       DEN=W/(H2-H4) >*/
    den = w / (h2 - h4);
/*<       CD2=H4*DEN >*/
    cd2 = h4 * den;
/*<       CC2=H2*DEN >*/
    cc2 = h2 * den;
/*<       W=CC2-CD1 >*/
    w = cc2 - cd1;
/*<       DEN=W/(H1-H4) >*/
    den = w / (h1 - h4);
/*<       DD1=H4*DEN >*/
    dd1 = h4 * den;
/*<       DC1=H1*DEN >*/
    dc1 = h1 * den;
/*<       If((H3+H4).lt.0D0) Then >*/
    if (h3 + h4 < 0.) {
/*<          Y=YA(4)+D3+CD2+DD1 >*/
	*y = ya[4] + d3 + cd2 + dd1;
/*<       Elseif((H2+H3).lt.0D0) Then >*/
    } else if (h2 + h3 < 0.) {
/*<          Y=YA(3)+D2+CD1+DC1 >*/
	*y = ya[3] + d2 + cd1 + dc1;
/*<       Elseif((H1+H2).lt.0D0) Then >*/
    } else if (h1 + h2 < 0.) {
/*<          Y=YA(2)+C2+CD1+DC1 >*/
	*y = ya[2] + c2 + cd1 + dc1;
/*<       ELSE >*/
    } else {
/*<          Y=YA(1)+C1+CC1+DC1 >*/
	*y = ya[1] + c1 + cc1 + dc1;
/*<       ENDIF >*/
    }
/*<       RETURN >*/
    return 0;
/*               ************************* */
/*<       END >*/
} /* polint4f_ */

/*<       Function NextUn() >*/
integer nextun_()
{
    /* System generated locals */
    integer ret_val;
    inlist ioin__1;

    /* Builtin functions */
    integer f_inqu(inlist *);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer n;
    static logical ex;

/*                                 Returns an unallocated FORTRAN i/o unit. */
/*<       Logical EX >*/

/*<       Do 10 N = 10, 300 >*/
    for (n = 10; n <= 300; ++n) {
/*<          INQUIRE (UNIT=N, OPENED=EX) >*/
	ioin__1.inerr = 0;
	ioin__1.inunit = n;
	ioin__1.infile = 0;
	ioin__1.inex = 0;
	ioin__1.inopen = &ex;
	ioin__1.innum = 0;
	ioin__1.innamed = 0;
	ioin__1.inname = 0;
	ioin__1.inacc = 0;
	ioin__1.inseq = 0;
	ioin__1.indir = 0;
	ioin__1.infmt = 0;
	ioin__1.inform = 0;
	ioin__1.inunf = 0;
	ioin__1.inrecl = 0;
	ioin__1.innrec = 0;
	ioin__1.inblank = 0;
	f_inqu(&ioin__1);
/*<          If (.NOT. EX) then >*/
	if (! ex) {
/*<             NextUn = N >*/
	    ret_val = n;
/*<             Return >*/
	    return ret_val;
/*<          Endif >*/
	}
/*<  10   Continue >*/
/* L10: */
    }
/*<       Stop ' There is no available I/O unit. ' >*/
    s_stop(" There is no available I/O unit. ", (ftnlen)33);
/*               ************************* */
/*<       End >*/
    return ret_val;
} /* nextun_ */

#ifdef __cplusplus
	}
#endif
