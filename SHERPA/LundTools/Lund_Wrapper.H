#ifndef SHERPA_LundTools_Lund_Wrapper_H
#define SHERPA_LundTools_Lund_Wrapper_H

#include "MODEL/Main/Running_AlphaS.H"
#include "SHERPA/LundTools/Lund_Interface.H"
#include "PDF/Main/ISR_Handler.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/CXXFLAGS.H"
#include <cstring>

using namespace SHERPA;

extern "C" {
  void pyhepc_(int*);
  void pyinit_(const char*,const char*,const char*,double*,int,int,int);
  void pyevnt_();
  void pyexec_();
  void pyprep_(int &);
  void pystrf_(int &);
  void pydecy_(int &);
  void pyinki_(int*);
  void pyinpr_();
  void pylist_(int*);
  void pyplst_();
  void pymult_(int*);
  double pyalp_(double*);
  double pyalps_(double*);
  int pycomp_(int*);
  void pypdpr_(double*,double*,double*);
  void pypdpr_(double*,double*,double*);
  void pyerrm_(int*);
}

void pyhepc(int mconv) { 
  pyhepc_(&mconv); 
}

int pycomp(int kfcode) { 
  return pycomp_(&kfcode); 
}

void pyinit(const char *frame,const char *beam,const char *target,double win) {
  pyinit_(frame,beam,target,&win,strlen(frame),strlen(beam),strlen(target));
}

void pyevnt() { 
  pyevnt_(); 
}

void pyexec() { 
  pyexec_(); 
}

void pyprep(int & callno) {
  pyprep_(callno);
}

void pystrf(int & callno) {
  pystrf_(callno);
}

void pydecy(int & callno) {
  pydecy_(callno);
}

void pyinki(int modki) { 
  pyinki_(&modki); 
}

void pyinpr() {
  pyinpr_();
}

void pylist(int dummy) { 
  pylist_(&dummy); 
}

void pymult(int mmul) {
  pymult_(&mmul);
}

double pyalps_(double *q2) {
  if (Lund_Interface::ExportAS()) return (*MODEL::as)(*q2);
  return pyalp_(q2);
}

void pyerrm_(int *error) {
  SHERPA::Lund_Interface::Error(*error);
}

void ListLundParameters() { 
  pyplst_(); 
}

const int nmaxjets  = 10000;
const int nmaxpart  = 500;
const int nmaxdecay = 8000;
const int nmaxpup   = 100;

const int nmxhep    = HEPEVT_CB_SIZE;

extern "C" {

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevt_;
#define hepevt hepevt_

extern struct {
  int n, npad;
  double k[5][nmaxjets], p[5][nmaxjets], v[5][nmaxjets];
} pyjets_;
#define pyjets pyjets_

extern struct {
  int mstu[200];
  double paru[200];
  int mstj[200];
  double parj[200];
} pydat1_;
#define pydat1 pydat1_

extern struct {
  int kchg[4][nmaxpart];
  double pmas[4][nmaxpart];
  double parf[2000];
  double vckm[4][4];
} pydat2_;
#define pydat2 pydat2_

extern struct {
  int mdcy[3][nmaxpart], mdme[2][nmaxdecay];
  double brat[nmaxdecay];
  int kfdp[5][nmaxdecay];
} pydat3_;
#define pydat3 pydat3_

extern struct {
  char *chaf[2][nmaxpart];
} pydat4_;
#define pydat4 pydat4_

extern struct {
  int msel, mselpd, msub[500], kfin[81][2];
  double ckin[200];
} pysubs_;
#define pysubs pysubs_

extern struct {
  int mstp[200];
  double parp[200];
  int msti[200];
  double pari[200];
} pypars_;
#define pypars pypars_

extern struct {
  int mrpy[6];
  double rrpy[100];
} pydatr_;
#define pydatr pydatr_

extern struct {
  int mint[400];
  double vint[400];
} pyint1_;
#define pyint1 pyint1_

extern struct {
  int iset[nmaxpart], kfpr[nmaxpart][2];
  double coef[20][nmaxpart];
  int icol[2][4][40];
} pyint2_;
#define pyint2 pyint2_

extern struct {
  double xsfx[81][2];
  int isig[3][1000];
  double sigh[1000];
} pyint3_;
#define pyint3 pyint3_

extern struct {
  int mwid[nmaxpart];
  double wids[5][nmaxpart];
} pyint4_;
#define pyint4 pyint4_

extern struct {
  int ngenpd, ngen[3][501];
  double xsec[3][501];
} pyint5_;
#define pyint5 pyint5_

extern struct {
  char *proc[501];
} pyint6_;
#define pyint6 pyint6_

extern struct {
  double sigt[7][7][6];
} pyint7_;
#define pyint7 pyint7_

extern struct {
  int imss[100];
  double rmss[100];
} pymssm_;
#define pymssm pymssm_

extern struct {
  double zmix[4][4], umix[2][2], vmix[2][2], smz[4], smw[2];
  double sfmix[4][16], zmixi[4][4], umixi[2][2], vmixi[2][2];
} pyssmt_;
#define pyssmt pyssmt_

extern struct {
  double rvlam[3][3][3], rvlamp[3][3][3], rvlamb[3][3][3];
} pymsrv_;
#define pymsrv pymsrv_

extern struct {
  int itcm[100];
  double rtcm[100];
} pytcsm_;
#define pytcsm pytcsm_

extern struct {
  int ihist[4], indx[1000];
  double bin[20000];
} pybins_;

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
