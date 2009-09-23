#ifndef Lund_Wrapper_H
#define Lund_Wrapper_H

#include "Running_AlphaS.H"
#include "Lund_Interface.H"
#include "ISR_Handler.H"
#include "Exception.H"
#include "CXXFLAGS.H"
#include <cstring>

using namespace SHERPA;

extern "C" {
  void sphepc_(int*);
  void spinit_(const char*,const char*,const char*,double*,int,int,int);
  void spevnt_();
  void spexec_();
  void spprep_(int &);
  void spstrf_(int &);
  void spdecy_(int &);
  void spinki_(int*);
  void spinpr_();
  void splist_(int*);
  void spplst_();
  void spmult_(int*);
  double spalp_(double*);
  double shalps_(double*);
  int spcomp_(int*);
  void sppdpr_(double*,double*,double*);
  void shpdpr_(double*,double*,double*);
  void sherrm_(int*);
  void psherp_();
}

void sphepc(int mconv) { 
  sphepc_(&mconv); 
}

int spcomp(int kfcode) { 
  return spcomp_(&kfcode); 
}

void spinit(const char *frame,const char *beam,const char *target,double win) {
  spinit_(frame,beam,target,&win,strlen(frame),strlen(beam),strlen(target));
}

void spevnt() { 
  spevnt_(); 
}

void spexec() { 
  spexec_(); 
}

void spprep(int & callno) {
  spprep_(callno);
}

void spstrf(int & callno) {
  spstrf_(callno);
}

void spdecy(int & callno) {
  spdecy_(callno);
}

void spinki(int modki) { 
  spinki_(&modki); 
}

void spinpr() {
  spinpr_();
}

void splist(int dummy) { 
  splist_(&dummy); 
}

void spmult(int mmul) {
  spmult_(&mmul);
}

double shalps_(double *q2) {
  if (Lund_Interface::ExportAS()) return (*MODEL::as)(*q2);
  return spalp_(q2);
}

void shpdpr_(double *x,double *q2,double *xppr) {
  if (Lund_Interface::ExportPDF()) {
    for (int i=0;i<=10;++i) {
      Lund_Interface::ISRHandler()->PDF(0)->Calculate(*x,*q2);
      if (i!=3 && i!=4) {
	xppr[i+1]=Lund_Interface::ISRHandler()->PDF(0)->GetXPDF(ATOOLS::Flavour((kf_code)(i-5)));
      }
      else {
	xppr[i+1]=Lund_Interface::ISRHandler()->PDF(0)->GetXPDF(ATOOLS::Flavour((kf_code)(2-i)));
      }
    }
  }
  else {
    sppdpr_(x,q2,xppr);
  }
}

void sherrm_(int *error) {
  SHERPA::Lund_Interface::Error(*error);
}

void psherp_() {
  ATOOLS::exh->GenerateStackTrace(std::cout);
}

void ListLundParameters() { 
  spplst_(); 
}

const int nmaxjets  = 10000;
const int nmaxpart  = 500;
const int nmaxdecay = 8000;
const int nmaxpup   = 100;

const int nmxhep    = HEPEVT_CB_SIZE;

extern "C" {
extern struct {
  int mfudge[100];
  double pfudge[100];
} cfudge_;
#define cfudge cfudge_

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevt_;
#define hepevt hepevt_

extern struct {
  int n, npad;
  double k[5][nmaxjets], p[5][nmaxjets], v[5][nmaxjets];
} spjets_;
#define spjets spjets_

extern struct {
  int mstu[200];
  double paru[200];
  int mstj[200];
  double parj[200];
} spdat1_;
#define spdat1 spdat1_

extern struct {
  int kchg[4][nmaxpart];
  double pmas[4][nmaxpart];
  double parf[2000];
  double vckm[4][4];
} spdat2_;
#define spdat2 spdat2_

extern struct {
  int mdcy[3][nmaxpart], mdme[2][nmaxdecay];
  double brat[nmaxdecay];
  int kfdp[5][nmaxdecay];
} spdat3_;
#define spdat3 spdat3_

extern struct {
  char *chaf[2][nmaxpart];
} spdat4_;
#define spdat4 spdat4_

extern struct {
  int msel, mselpd, msub[500], kfin[81][2];
  double ckin[200];
} spsubs_;
#define spsubs spsubs_

extern struct {
  int mstp[200];
  double parp[200];
  int msti[200];
  double pari[200];
} sppars_;
#define sppars sppars_

extern struct {
  int mrsp[6];
  double rrsp[100];
} spdatr_;
#define spdatr spdatr_

extern struct {
  int mint[400];
  double vint[400];
} spint1_;
#define spint1 spint1_

extern struct {
  int iset[nmaxpart], kfpr[nmaxpart][2];
  double coef[20][nmaxpart];
  int icol[2][4][40];
} spint2_;
#define spint2 spint2_

extern struct {
  double xsfx[81][2];
  int isig[3][1000];
  double sigh[1000];
} spint3_;
#define spint3 spint3_

extern struct {
  int mwid[nmaxpart];
  double wids[5][nmaxpart];
} spint4_;
#define spint4 spint4_

extern struct {
  int ngenpd, ngen[3][501];
  double xsec[3][501];
} spint5_;
#define spint5 spint5_

extern struct {
  char *proc[501];
} spint6_;
#define spint6 spint6_

extern struct {
  double sigt[7][7][6];
} spint7_;
#define spint7 spint7_

extern struct {
  int imss[100];
  double rmss[100];
} spmssm_;
#define spmssm spmssm_

extern struct {
  double zmix[4][4], umix[2][2], vmix[2][2], smz[4], smw[2];
  double sfmix[4][16], zmixi[4][4], umixi[2][2], vmixi[2][2];
} spssmt_;
#define spssmt spssmt_

extern struct {
  double rvlam[3][3][3], rvlamp[3][3][3], rvlamb[3][3][3];
} spmsrv_;
#define spmsrv spmsrv_

extern struct {
  int itcm[100];
  double rtcm[100];
} sptcsm_;
#define sptcsm sptcsm_

extern struct {
  int ihist[4], indx[1000];
  double bin[20000];
} spbins_;

extern struct {
  int idbmup[2];
  double ebmup[2];
  int pdfgup[2], pdfsup[2], idwtup, nprup;
  double xsecup[nmaxpup], xerrup[nmaxpup], xmaxup[nmaxpup];
  int lprup[nmaxpup];
} heprup_;
#define heprup heprup_

}
#endif
