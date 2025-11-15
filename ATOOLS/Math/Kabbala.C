#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Org/Exception.H"

using namespace ATOOLS;

Kabbala& Kabbala::operator+=(const Kabbala& k) {
  double a = abs(rishbon);
  double b = abs(k.Value());
  double max = ATOOLS::Max(a,b);
  if (max==0 || ATOOLS::IsZero(b/max)) {
    return *this;
  }
  if (ATOOLS::IsZero(a/max)) {
    rishbon = k.Value();
    shem = k.String();
    lambda = k.Lambda();
    return *this;
  }
  rishbon += k.Value();
  if (ATOOLS::IsZero(abs(rishbon)/max)) {
    return *this = Zero();
  }

  if (shem!=std::string("")) shem += std::string("+");
  shem    += k.String();  
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) + copy2(m);};
  return *this;
}

Kabbala& Kabbala::operator-=(const Kabbala& k) {
  double a = abs(rishbon);
  double b = abs(k.Value());
  double max = ATOOLS::Max(a,b);
  if (max==0 || ATOOLS::IsZero(b/max)) {
    return *this;
  }
  if (ATOOLS::IsZero(a/max)) {
    this->Negate();
    return *this;
  }
  rishbon -= k.Value();
  if (ATOOLS::IsZero(abs(rishbon)/max)) {
    return *this = Zero();
  }
  shem    += std::string("-(");
  shem    += k.String();
  shem    += std::string(")");
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) - copy2(m);};
  return *this;
}

Kabbala Kabbala::operator-() {
  this->Negate();
  return *this;
}

Kabbala Kabbala::operator+() {
  return Kabbala(shem,rishbon,lambda);
}

Kabbala& Kabbala::operator*=(const Kabbala& k) {
  if (rishbon==C_ZERO) return *this;
  if (k.Value()==C_ZERO) {
    return *this = Zero();
  }
  rishbon *= k.Value();
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
  shem += k.String();
  shem += std::string(")");  
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) * copy2(m);};
  return *this;
}

Kabbala& Kabbala::operator*=(const int& i) {
  if (i==0) return *this = Zero();
  rishbon *= i;
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
  MyStrStream sstr;
  sstr<<i;
  std::string istr;
  sstr>>istr;
  shem += istr;
  shem += std::string(")");
  Func copy(lambda);
  int f(i);
  lambda = [f, copy](Function_Argument m) {return Value_T(f, .0)*copy(m);};
  return *this;
}

Kabbala& Kabbala::operator*=(const double& d) {
  if (d==.0) return *this = Zero();
  rishbon *= d;
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
  MyStrStream sstr;  
  sstr<<d;
  std::string dstr;
  sstr>>dstr;
  shem += dstr;
  shem += std::string(")");
  Func copy(lambda);
  double f(d);
  lambda = [f, copy](Function_Argument m) {return Value_T(f, .0)*copy(m);};
  return *this;
}

Kabbala& Kabbala::operator*=(const Value_T& c) {
  if (c == C_ZERO) return *this = Zero();
  rishbon *= c;
  std::string save = shem;
  MyStrStream sstr;  
  sstr<<"("<<save<<")*("<<c.real()<<"+i*("<<c.imag()<<"))";
  sstr >> shem;
  Func copy(lambda);
  Value_T f(c);
  lambda = [f, copy](Function_Argument m) {return f*copy(m);};
  return *this;
}

Kabbala& Kabbala::operator/=(const Kabbala& k) {
  if (abs(k.Value()) == 0.0) THROW(fatal_error, "division by zero :(");
  if (abs(rishbon) == 0.0) return *this;
  rishbon /= k.Value();
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")/(");
  shem += k.String();
  shem += std::string(")");
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) / copy2(m);};
  return *this;
}
