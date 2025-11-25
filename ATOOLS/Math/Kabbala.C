#include "ATOOLS/Math/Kabbala.H"
#include "ATOOLS/Org/Exception.H"
#include "Kabbala.H"

namespace ATOOLS {

//Default constructor
Kabbala::Kabbala(){
  shem = "0";
  rishbon = C_ZERO;
  lambda = [](Function_Argument map) {return C_ZERO;};
}
//Constructor for fixed constant.
Kabbala::Kabbala(Complex c){
  MyStrStream ss;
  ss << c;
  ss >> shem;
  rishbon = c;
  Complex c1 (c);
  lambda = [c1](Function_Argument map) {return c1;};
}
// constructor without using initial value, that will be calculated from the function
Kabbala::Kabbala(std::string str, Func func, Function_Argument map){
  shem = str;
  lambda = func;
  Update(map);
}
// legacy constructor, reverts lambda to basic look up of the given String
Kabbala::Kabbala(std::string str ,Complex C) {
  shem    = str;
  rishbon = C;
  msg_Debugging() << "No proper function set for the Kabbala. Choosing Basic Lookup." << std::endl;
  lambda = BasicLookUpFunction();
}
// constructor with lambda and initial value, these may be different values, the value is NOT updated
Kabbala::Kabbala(std::string str, Complex C, Func func){
  shem = str;
  lambda = func;
  rishbon = C;
}
// copy constructor, this may be the reason the lambdas need to be copied explicitly edit: no
Kabbala::Kabbala(const Kabbala& k) {
  shem    = k.String();
  rishbon = k.Value();
  lambda = k.Lambda();
}


Complex Kabbala::Update(Function_Argument map){return rishbon = lambda(map);}

bool Kabbala::DependsOn(std::string param) const {return shem.find(param) != std::string::npos;}

Kabbala::Func Kabbala::BasicLookUpFunction(){
  std::string copy(shem);
  return [copy](Function_Argument map) {
    if (map->count(copy) == 0) return C_ZERO; 
    return Complex(map->at(copy), .0); 
  };
}

// Copy/Assignment
Kabbala& Kabbala::operator=(const Kabbala &k){
  if (this!=&k) {
    shem    = k.String();
    rishbon = k.Value();
    lambda = k.Lambda();
  } 
  return *this;
}

// in place add
Kabbala& Kabbala::operator+=(const Kabbala& k) {
  rishbon += k.Value();
  if (shem!=std::string("")) shem += std::string("+");
  shem    += k.String();  
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) + copy2(m);};
  return *this;
}

Kabbala& Kabbala::operator+=(const Complex& c){
  if (ATOOLS::IsZero(c)) return *this;
  rishbon += c;
  if (shem!=std::string("")) shem += std::string("+");
  MyStrStream ss;
  ss << "(" << c << ")";
  std::string s;
  ss >> s;
  shem += s;  
  Func copy1(lambda);
  Complex c2(c);
  lambda = [copy1, c2](Function_Argument m) {return copy1(m) + c2;};
  return *this;
}

// in place subtract
Kabbala Kabbala::operator-() {
  rishbon = -rishbon;
  shem = std::string("-(")+shem+std::string(")");
  Func copy(lambda);
  lambda = [copy](Function_Argument m) {return -copy(m);};
  return *this;
}

Kabbala& Kabbala::operator-=(const Kabbala& k) {
  rishbon -= k.Value();
  shem    += std::string("-(");
  shem    += k.String();
  shem    += std::string(")");
  Func copy1(lambda);
  Func copy2(k.Lambda());
  lambda = [copy1, copy2](Function_Argument m) {return copy1(m) - copy2(m);};
  return *this;
}

Kabbala& Kabbala::operator-=(const Complex& c) {
  if (ATOOLS::IsZero(c)) return *this;
  rishbon -= c;
  MyStrStream ss;
  ss << "-(" << c << ")";
  std::string s;
  ss >> s;
  shem    += s;
  Func copy1(lambda);
  Complex c2(c);
  lambda = [copy1, c2](Function_Argument m) {return copy1(m) - c2;};
  return *this;
}

// in place multiply
Kabbala& Kabbala::operator*=(const Kabbala& k) {
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

Kabbala& Kabbala::operator*=(const Complex& c) {
  if (ATOOLS::IsZero(c)) return *this = Kabbala();
  if (c == Complex(1., 0.)) return *this;
  rishbon *= c;
  std::string save = shem;
  MyStrStream sstr;  
  sstr<<"("<<save<<")*("<< c <<")";
  sstr >> shem;
  Func copy(lambda);
  Complex c2(c);
  lambda = [c2, copy](Function_Argument m) {return c2*copy(m);};
  return *this;
}

// Division in place, with constant
Kabbala& Kabbala::operator/=(const Kabbala& k) {
  // here the check is necessary, this is still unsafe if updated
  if (abs(k.Value()) == 0.0) THROW(fatal_error, "division by zero :(");
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

Kabbala& Kabbala::operator/=(const Complex& c) {
  if (ATOOLS::IsZero(c)) THROW(fatal_error, "division by zero :(");
  if (c == Complex(1., 0.)) return *this;
  rishbon *= c;
  std::string save = shem;
  MyStrStream sstr;  
  sstr<<"("<<save<<")/("<< c <<")";
  sstr >> shem;
  Func copy(lambda);
  Complex c2(c);
  lambda = [c2, copy](Function_Argument m) {return c2*copy(m);};
  return *this;
}

Kabbala operator/(const Complex& c, const Kabbala& k1) {
  if (abs(k1.Value()) == 0.0) THROW(fatal_error, "division by zero :(");
  if (ATOOLS::IsZero(c)) return Kabbala();
  Kabbala k(k1);
  k.SetValue(c/k1.Value());
  MyStrStream ss;
  ss << "(" << c << ")/(" << k1.String() << ")";
  std::string s;
  ss >> s;
  k.SetString(s);
  Kabbala::Func copy1(k1.Lambda());
  Complex c1(c);
  k.SetLambda([copy1, c1](Kabbala::Function_Argument m) {return c1/copy1(m);});
  return k;
}

// compare
bool operator==(const Kabbala& k1, const Kabbala& k2){
  if (k1.Value()!=k2.Value())   return false;
  if (k1.String()!=k2.String()) return false;
  return true;
}

bool operator<(const Kabbala& k1, const Kabbala& k2){
  if (k1.Value()==k2.Value()) return k1.String() < k2.String();
  return abs(k1.Value()) < abs(k2.Value());
}

// other operations
Kabbala exp(const Kabbala& k1) {
  Kabbala k(k1);
  k.SetValue(std::exp(k.Value()));
  k.SetString(std::string("exp(") + k.String() + std::string(")"));
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return std::exp(copy(m));});
  return k;
}

Kabbala sin(const Kabbala& k1) {
  Kabbala k(k1);
  k.SetValue(std::sin(k.Value()));
  k.SetString(std::string("sin(") + k.String() + std::string(")"));
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return std::sin(copy(m));});
  return k;
}

Kabbala cos(const Kabbala& k1) {
  Kabbala k(k1);
  k.SetValue(std::cos(k.Value()));
  k.SetString(std::string("cos(") + k.String() + std::string(")"));
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return std::cos(copy(m));});
  return k;
}

Kabbala sqrt(const Kabbala& k1) {
  Kabbala k(k1);
  k.SetValue(std::sqrt(k.Value()));
  k.SetString(std::string("sqrt(") + k.String() + std::string(")"));
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return std::sqrt(copy(m));});
  return k;
}

Kabbala pow(const Kabbala& k1, const Complex& c) {
  Kabbala k(k1);
  k.SetValue(std::pow(k.Value(), c));
  MyStrStream ss;
  ss << "(" << k.String() << ")^(" << c << ")";
  std::string s;
  ss >> s;
  k.SetString(s);
  Kabbala::Func copy(k.Lambda());
  Complex c1 (c);
  k.SetLambda([copy, c1](Kabbala::Function_Argument m){return std::pow(copy(m), c1);});
  return k;
}

Kabbala complexconjugate(const Kabbala& k1){
  Kabbala k(k1);
  k.SetValue(Complex(k.Value().real(), -k.Value().imag()));
  k.SetString("(" + k.String() +")*");
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return Complex(copy(m).real(), -copy(m).imag());});
}

/*Kabbala abs(const Kabbala& k1) {
  Kabbala k(k1);
  k.SetValue(std::abs(k.Value()));
  k.SetString(std::string("|") + k.String() + std::string("|"));
  Kabbala::Func copy(k.Lambda());
  k.SetLambda([copy](Kabbala::Function_Argument m){return std::abs(copy(m));});
  return k;
}*/
}