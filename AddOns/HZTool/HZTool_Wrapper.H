#ifndef AddOns_HZTool_HZTool_Wrapper_H
#define AddOns_HZTool_HZTool_Wrapper_H

#include "ATOOLS/Org/CXXFLAGS.H"
#include <cstring>

inline void MakeFortranString
(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) 
    output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) 
    output[j]=(char)input[j];
}

extern "C" {
  void hzinit_(double *xslo,double *nchavg);
  void hzevnt_(double *wgt);
  void hzfinl_();
}

void hzinit(double xslo,double nchavg) 
{ 
  hzinit_(&xslo,&nchavg); 
}

void hzevnt(double wgt) 
{
  hzevnt_(&wgt); 
}

void hzfinl() 
{
  hzfinl_(); 
}

#if HEPEVT_CB_SIZE != 4000
#error Wrong HepEvt common block size. Must be 4000 when using HZTool.
#error Use './configure --enable-hepevtsize=4000' to reconfigure.
#endif

const int nmxhep = HEPEVT_CB_SIZE;

extern "C" {

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevt_;
#define hepevt hepevt_

extern struct {
  int nevhep, nhep, isthep[nmxhep], idhep[nmxhep];
  int jmohep[nmxhep][2], jdahep[nmxhep][2];
  double phep[nmxhep][5], vhep[nmxhep][4];
} hepevtp_;
#define hepevtp hepevtp_

extern struct {
  double xsec;
  char gen[8];
  int ihchrg[nmxhep];
  float ntot, wtx;
  char beams[8][2];
} heracmn_;
#define heracmn heracmn_

extern struct {
  char hname[128];
} hzhname_;
#define hzhname hzhname_

}
#endif
