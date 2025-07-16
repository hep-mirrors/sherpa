#ifndef WS_MATH_UT
#define WS_MATH_UT
//#include <stdlib.h>

typedef double real_type;

//***********************************************************************
template <class T> int outrange(const T& x, const T& xmin, const T& xmax)
{return (x < xmin) ? -1 : ((x >= xmax)? 1 : 0);}

//-----------------------------------------------------------------------
class TLimits {
public:
  real_type lo, hi;
protected:
  int started;

public:
  TLimits(){ started = 0;}
  void Reset(){started = 0;}
  void Reset(real_type xl, real_type xh){
    lo = xl;
    hi = xh;
    started = 1;
  }

  real_type Lo(){
    if(!started) {
      cerr << "TLimits::Lo ERROR: not started." << endl;
      exit(101);
    }
    return lo;
  }

  real_type Hi(){
    if(!started) {
      cerr << "TLimits::Hi ERROR: not started." << endl;
      exit(101);
    }
    return hi;
  }

  void IncrLo(real_type x){
    if(x > lo) lo = x;
  }
  void DecrHi(real_type x){
    if(x < hi) hi = x;
  }

  void Intersect(real_type xl, real_type xh){
    if(!started) return;
    IncrLo(xl);
    DecrHi(xh);
  }

  void Intersect(const TLimits& lims){
    if(!started) return;
    IncrLo(lims.lo);
    DecrHi(lims.hi);
  }

  TLimits Scale(real_type s){
    if(!started) return *this;
    if(s >= 0.0) {lo = s*lo; hi = s*hi;}
    else {real_type lo1 = s*hi; hi = s*lo; lo = lo1;}
    return *this;
  }

  void Expand(real_type x){
    if(started){
      if(x < lo) lo = x;
      else if(x > hi) hi = x;
    }
    else {lo = hi = x; started = 1;}
  }

  int InRange(real_type x){
    return started && (lo < x) && (x < hi);
  }

  int InRangeR(real_type x){
    return started && (lo < x) && (x <= hi);
  }

  int InRangeL(real_type x){
    return started && (lo <= x) && (x < hi);
  }

  int InRangeLR(real_type x){
    return started && (lo <= x) && (x <= hi);
  }

  friend ostream& operator << (ostream& ops, TLimits rr){
    return ops << "[" << rr.Lo() << ", "
        << rr.Hi() << "] ";
  }

};

//-----------------------------------------------------------------------
struct TTableVar {
  real_type step, start, end, v;
  int logstep;

  TTableVar(real_type s, real_type e, int np, int logst=0){
    start = s;
    end = e;
    logstep = logst? 1 : 0;
    if(np > 1) step = logstep? log(end/start)/(np-1) : (end-start)/(np-1);
    else step = logstep;
  }

  TTableVar(const TLimits& lims, int np, int logst=0){
    *this = TTableVar(lims.lo, lims.hi, np, logst);
  }

  real_type Val(int ind){
    return v = logstep? start*exp(step*ind) : (start + ind*step);
  }
};

#endif
