#ifndef AddOns_HAWK_HAWK_Wrapper_H
#define AddOns_HAWK_HAWK_Wrapper_H

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/MyComplex.H"
#include <cstring>

#define HAWK_NF 5
#define HAWK_NP 6

extern "C" {
  extern struct {
    double mudim2;
  } uv_;
  extern struct {
    double lambda,lambda2,deltas;
  } ir_;
  extern struct {
    double gevfb;
  } gevfb_;
  //extern struct {
  //  double qu,qd,ql,qn,qf[4],mu,mu2,md,md2;
  //  Complex guu[3], gdd[3];
  //} qf_;
  // this common block does not need to be filled by us.
  extern struct {
    int dalpz,drbos,drferm,qalp;
  } rscheme_;
  extern struct {
    Complex ccw,ccw2,csw,csw2,cmw,cmw2,cmz,cmz2;
  } cparam_;
  extern struct {
    Complex xcw,xcw2,xsw,xsw2,xmw,xmw2,xmz,xmz2,xmh,xmh2,zero;
    Complex xml[3],xml2[3],xmqp[3],xmqp2[3],xmqm[3],xmqm2[3];
  } xparam_;
  extern struct {
    double  me2,rme2,mf2[5],rmf2[5];
    Complex cme2,cmf2[5];
  } fmass_;
  extern struct {
    double constff,constfi,constif,constii;
  } constsub_;
  extern struct {
    double chh2;
  } hh2_;
  extern struct {
    double pi,el,alpha,alpha0,alphaZ,GF,alphas;
    Complex v[9],cv[9];
    double cw,cw2,sw,sw2,mw,mw2,gw,mz,mz2,gz,mh,mh2;
    double ml[3],ml2[3],mqp[3],mqp2[3],mqm[3],mqm2[3];
  } param_;
  extern struct {
    int qborn,qw,qz,qschan,qtchan,qch2,qchint;
    int qbini,qbfin,qwidth,qfact,qbos,qferm,qsoft,qhh2;
    int qqcddiag,qqcdnondiag,qqcdgsplit,qqcdggfus,qcp;
  } rcoptions_;
}

#endif
