#include "gpl.h"
#include <ginac/ginac.h>

using namespace GiNaC;

double GPL(double w1, double arg) {
  ex x1 = numeric(w1);
  ex xarg = numeric(arg);
  lst l = {x1};
  ex res = evalf(G( l, xarg ));
  return real(ex_to<numeric>(res)).to_double();
};

double GPL(double w1, double w2, double arg) {
  ex x1 = numeric(w1);
  ex x2 = numeric(w2);
  ex xarg = numeric(arg);
  lst l = {x1,x2};
  ex res = evalf(G( l, xarg ));
  return real(ex_to<numeric>(res)).to_double();
};
