#include "ATOOLS/Math/Kabbala.H"

using namespace ATOOLS;


Kabbala& Kabbala::operator+=(const Kabbala& k) {
  double a = abs(rishbon);
  double b = abs(k.Value());
  double max = ATOOLS::Max(a,b);
  if (ATOOLS::IsZero(b/max)) {
    return *this;
  }
  if (ATOOLS::IsZero(a/max)) {
    rishbon = k.Value();
    shem = k.String();
    return *this;
  }
  rishbon += k.Value();
  if (ATOOLS::IsZero(abs(rishbon)/max)) {
    shem = std::string("0");
    rishbon = Complex(0.,0.);
    return *this;
  }

  if (shem!=std::string("")) shem += std::string("+");
  shem    += k.String();
      
  return *this;
}
    
Kabbala& Kabbala::operator-=(const Kabbala& k) {
  double a = abs(rishbon);
  double b = abs(k.Value());
  double max = ATOOLS::Max(a,b);
  if (ATOOLS::IsZero(b/max)) {
    return *this;
  }
  if (ATOOLS::IsZero(a/max)) {
    rishbon = -k.Value();
    shem = std::string("-(")+k.String()+std::string(")");
    return *this;
  }
  rishbon -= k.Value();
  if (ATOOLS::IsZero(abs(rishbon)/max)) {
    shem = std::string("0");
    rishbon = Complex(0.,0.);
    return *this;
  }
  shem    += std::string("-(");
  shem    += k.String();
  shem    += std::string(")");
      
  return *this;
}
    
Kabbala Kabbala::operator-() {
  return Kabbala(std::string("-")+std::string("(")+shem+std::string(")"),-rishbon);
}

Kabbala Kabbala::operator+() {
  return Kabbala(shem,rishbon);
}
    
Kabbala& Kabbala::operator*=(const Kabbala& k) {
  if (rishbon==Complex(0.,0.)) return *this;
  if (k.Value()==Complex(0.,0.)) {
    shem = std::string("0");
    rishbon = Complex(0.,0.);
    return *this;
  }
  rishbon *= k.Value();
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
  shem += k.String();
  shem += std::string(")");
  return *this;
}
    
Kabbala& Kabbala::operator*=(const int& i) {
  rishbon *= i;
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
      
  MyStrStream sstr;
  sstr<<i;
  std::string istr;
  sstr>>istr;
      
  shem += istr;
  shem += std::string(")");
  return *this;
}

Kabbala& Kabbala::operator*=(const double& d) {
  rishbon *= d;
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")*(");
      
  MyStrStream sstr;  
  sstr<<d;
  std::string dstr;
  sstr>>dstr;
      
  shem += dstr;
  shem += std::string(")");
  return *this;
}

Kabbala& Kabbala::operator/=(const Kabbala& k) {
  rishbon /= k.Value();
  std::string save = shem;
  shem  = std::string("(") + save + std::string(")/(");
  shem += k.String();
  shem += std::string(")");
  return *this;
}
