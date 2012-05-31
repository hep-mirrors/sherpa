#ifndef AddOns_MCFM_MCFM_Wrapper_H
#define AddOns_MCFM_MCFM_Wrapper_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include <cstring>

inline void MakeFortranString
(char *output,std::string input,unsigned int length)
{
  for (unsigned int i=0;i<length;++i) 
    output[i]=(char)32;
  for (size_t j=0;j<input.length();++j) 
    output[j]=(char)input[j];
}

#define MCFM_NF 5
#define MCFM_NMX 12

inline int mp(const int id,const int i)
{ return i*MCFM_NMX+id; }

inline int mr(const int i,const int j)
{ return (j+MCFM_NF)*(2*MCFM_NF+1)+(i+MCFM_NF); }

inline void GetMom(double *m,const int n,const ATOOLS::Vec4D &p)
{ m[mp(n,3)]=p[0]; for (int i(1);i<4;++i) m[mp(n,i-1)]=p[i]; }

extern "C" {
  void coupling_();
  void coupling2_();

  extern struct {
    bool verbose;
  } verbose_;
  extern struct {
    double epinv;
  } epinv_;
  extern struct {
    double epinv2;
  } epinv2_;
  extern struct {
    char scheme[4];
  } scheme_;
  extern struct {
    double deltar;
  } deltar_;
  extern struct {
    double scale,musq;
  } scale_;
  extern struct {
    double gsq,as,ason2pi,ason4pi;
  } qcdcouple_;
  extern struct {
    double gf_inp,aemmz_inp,xw_inp,wmass_inp,zmass_inp;
  } ewinput_;
  extern struct {
    int ewscheme;
  } ewscheme_;
  extern struct {
    double Gf,gw,xw,gwsq,esq,vevsq;
  } ewcouple_;
  extern struct {
    int nflav;
  } nflav_;
  extern struct {
    double md,mu,ms,mc,mb,mt,
      mel,mmu,mtau,
      hmass,hwidth,
      wmass,wwidth,
      zmass,zwidth,
      twidth,
      tauwidth,
      mtausq,mcsq,mbsq;
  } masses_;
  extern struct {
    int n2,n3;
    double mass2,width2,mass3,width3;
  } breit_;
  extern struct {
    double Vud,Vus,Vub,Vcd,Vcs,Vcb;
  } cabib_;
  extern struct {
    double susycoup;
  } susycoup_;
  extern struct {
    int nlooprun;
  } nlooprun_;
  extern struct {
    char pdlabel[7];
  } pdlabel_;
  extern struct {
    double amz;
  } couple_;
  extern struct {
    double mb_msbar;
  } mb_msbar_;
  extern struct {
    int nproc;
  } nproc_;
  extern struct {
    bool removebr;
  } removebr_;
  extern struct {
    bool zerowidth;
  } zerowidth_;
  extern struct {
    double bbsqmin,bbsqmax,wsqmin,wsqmax;
  } limits_;
  extern struct {
    Complex za[MCFM_NMX*MCFM_NMX], zb[MCFM_NMX*MCFM_NMX];
  } zprods_;
  extern struct{
    int nwz;
  }nwz_;
  extern struct{
    double s[MCFM_NMX*MCFM_NMX];
  }sprods_;
  extern struct{
    double delg1_z,delg1_g,lambda_g,lambda_z,
      h1Z,h2Z,h3Z,h4Z,h1gam,h2gam,h3gam,h4gam,
      delk_g,delk_z,tevscale;
    double xdelg1_z,xdelg1_g,xlambda_g,xlambda_z,
      h1tZ,h2tZ,h3tZ,h4tZ,h1tgam,h2tgam,h3tgam,h4tgam,
      xdelk_g,xdelk_z;
    bool anomtgc;
    }anomcoup_;
  
}

#endif
