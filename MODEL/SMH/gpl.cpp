#include "gpl.h"
#include <ginac/ginac.h>

extern "C" {
  struct dpx { double r, i;
    inline dpx(const double &_r=0.,const double &_i=0.):
      r(_r), i(_i) {} };
  dpx gpl1_(dpx *w1,dpx *x);
  dpx gpl2_(dpx *w1,dpx *w2,dpx *x);
}

using namespace GiNaC;

double GPL(double w1, double arg) {
  dpx _w1(w1), _x(arg);
  return gpl1_(&_w1,&_x).r;
};

double GPL(double w1, double w2, double arg) {
  dpx _w1(w1), _w2(w2), _x(arg);
  return gpl2_(&_w1,&_w2,&_x).r;
};
