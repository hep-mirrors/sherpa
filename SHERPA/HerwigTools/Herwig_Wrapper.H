#ifndef SHERPA_HerwigTools_Herwig_Wrapper_H
#define SHERPA_HerwigTools_Herwig_Wrapper_H

#include <ctype.h>
#include <cstring>

const int nmxhep    = 4000;

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevt_;
#define hepevt hepevt_

extern struct {
  double asfixd, clq[7][6], coss, costh, ctmax, disf[13][2], emlst;
  double emmax, emmin, empow, emsca, epoln[3], gcoef[7], gpoln, omega0, phomas;
  double ppoln[3],ptmax, ptmin, ptpow, q2max, q2min, q2pow, q2wwmn, q2wwmx, qlim;
  double sins, thmax, y4jt, tmnisr, tqwt, xx[2], xlmin, xxmin, ybmax, ybmin, yjmax;
  double yjmin, ywwmax, ywwmin, whmin, zjmax, zmxisr;
  int iaphig, ibrn[2], ibsh;
  int ico[10], idcmf, idn[10], iflmax, iflmin, ihpro, ipro, mapq[6], maxfl;
  bool bgshat, colisr, fstevt, fstwgt, genev, hvfcen, tpol, durham;
} hwhard_;
#define hwhard hwhard_

extern struct {
  double avwgt, evwgt, gamwt, tlout, wbigst;
  double wgtmax, wgtsum, wsqsum;
  int idhw[nmxhep], ierror, istat, lwevt;
  int maxer, maxpr;
  bool nowgt;
  int nrn[2], numer;
  int numeru, nwgts;
  bool gensof;
} hwevnt_;
#define hwevnt hwevnt_

extern struct {
  char part1[8], part2[8];
} hwbmch_;
#define hwbmch hwbmch_

extern struct {
  int ipart1, ipart2;
} hwbeam_;
#define hwbeam hwbeam_

extern struct {
  double ebeam1, ebeam2, pbeam1, pbeam2;
  int iproc, maxev;
} hwproc_;
#define hwproc hwproc_

const int maxhrp = 100;

extern struct {
  bool lhsoft, lhglsf;
  double lhwgt[maxhrp], lwgts[maxhrp], lhxsct[maxhrp];
  double lhxerr[maxhrp], lhxmax[maxhrp], lhmxsm;
  int lhiwgt[maxhrp], lhnevt[maxhrp], ityplh;
} hwgupr_;
#define hwgupr hwgupr_

extern struct {
  char autpdf[2][20], bdecay[4];
} hwprch_;
#define hwprch hwprch_

const int modmax = 5;

extern struct {
  double alpfac, brhig[12], enhanc[12], gammax, rhohep[3][nmxhep];
  int iophig, modbos[modmax];
} hwbosc_;
#define hwbosc hwbosc_ 

extern struct {
  double afch[16][2], alphem, b1lim, betaf, btclm;
  double cafac, cffac, clmax, clpow, clsmr[2], cspeed;
  double ensof, etamix, f0mix, f1mix, f2mix;
  double gamh, gamw, gamz, gamzp;
  double gev2nb, h1mix, pdiqk, pgsmx, pgspl[4], phimix, pifac;
  double prsof, psplt[2], ptrms, pxrms, qcdl3, qcdl5, qcdlam, qdiqk, qfch[16], qg;
  double qspac, qv, scabi, swein, tmtop, vfch[16][2], vckm[3][3], vgcut, vqcut;
  double vpcut, zbinm, effmin, omhmix, et2mix, ph3mix, gcutme;
  int ioprem, iprint, ispac, lrsud, lwsud, modpdf[2], nbtry, ncolo, nctry;
  int ndtry, netry, nflav, ngspl, nstru, nstry, nzbin, iop4jt[2], nprfmt;
  bool azsoft, azspin, cldir[2], hardme, nospac, prndec, prvtx, softme;
  bool zprime, prndef, prntex, prnweb;
} hwpram_;
#define hwpram hwpram_

const int nmxres = 500;

extern struct {
  double rltim[nmxres], rmass[nmxres], rspin[nmxres];
  int ichrg[nmxres], idpdg[nmxres], iflav[nmxres], nres;
  bool vtocdk[nmxres], vtordk[nmxres];
  bool qorqqb[nmxres], qborqq[nmxres];
} hwprop_;
#define hwprop hwprop_

#define hwcint hwcint_
#define hwcgse hwcgse_
#define hwctrm hwctrm_
#define hwuepr hwuepr_

extern "C" {
  void hwcint();
  void hwcgse();
  void hwctrm();
  void hwuepr();
}

inline void hwinit() { hwcint(); }

//////////////////
// Doxygen part //
//////////////////

/*!
  \file Herwig65Wrapper.H
  \brief Contains definitions of \b HERWIG common blocks
*/

#endif  
