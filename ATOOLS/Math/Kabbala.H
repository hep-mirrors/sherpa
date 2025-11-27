#ifndef ATOOLS_Math_Kabbala_H
#define ATOOLS_Math_Kabbala_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H" 
#include "ATOOLS/Org/MyStrStream.H"
#include <functional>
#include <map>

#define C_ZERO Complex(0., 0.)

namespace ATOOLS {
  /*
  Class for storing (complex) constants which can be updated later on.

  Each object can update based on external real constants, which need to be passed.
  */
  class Kabbala {
    public:
      // types for the Lambda Functions for better readability
      using Function_Argument = std::map<std::string, double_t>*;
      using Func = std::function<Complex(Function_Argument)>;
    public:
      // members are a string, the current value and the lambda to update this value
      std::string  shem;
      Complex rishbon;
      Func lambda = [](Function_Argument map) {return C_ZERO;}; //needs initial value to prevent bad function call
    public:
      //constructors
      Kabbala();
      Kabbala(Complex C);
      Kabbala(std::string str, Func func, Function_Argument map);
      Kabbala(std::string str ,Complex C);
      Kabbala(std::string str, Complex C, Func func);
      Kabbala(const Kabbala& k);

      // Stuff for the lambda functions
      Complex Update(Function_Argument map);
      bool DependsOn(std::string param) const;
      Func BasicLookUpFunction();

      // Member access
      Complex Value() const {return rishbon;} 
      std::string String() const {return shem;}
      Func Lambda() const {return lambda;}
    
      // Member setters
      void SetValue(Complex val) {rishbon = val;}
      void SetString(std::string str) {shem = str;} 
      void SetLambda(Func func) {lambda = func;}

      /*
      a LOT of useful operators :)
      */
      Kabbala& operator=(const Kabbala& k);

      Kabbala& operator+=(const Kabbala& k);
      Kabbala& operator+=(const Complex& c);
      Kabbala& operator+=(const int& i) {return *this += Complex(i);}
      Kabbala& operator+=(const double& d) {return *this += Complex(d);}

      Kabbala& operator-=(const Kabbala& k);
      Kabbala& operator-=(const Complex& c);
      Kabbala& operator-=(const int& i) {return *this -= Complex(i);}
      Kabbala& operator-=(const double& d) {return *this -= Complex(d);}

      Kabbala& operator*=(const Kabbala& k);
      Kabbala& operator*=(const Complex& c);
      Kabbala& operator*=(const int& i) {return *this *= Complex(i);}
      Kabbala& operator*=(const double& d) {return *this *= Complex(d);}

      Kabbala& operator/=(const Kabbala& k);
      Kabbala& operator/=(const Complex& c);
      Kabbala& operator/=(const int& i) {return *this /= Complex(i);}
      Kabbala& operator/=(const double& d) {return *this /= Complex(d);}
      
      Kabbala operator+() {return Kabbala(shem,rishbon,lambda);}
      Kabbala operator+(const Kabbala& k) const {Kabbala K(*this); K+=k; return K;}
      Kabbala operator+(const Complex& c) const {Kabbala K(*this); K+=c; return K;}
      Kabbala operator+(const int& i)     const {Kabbala K(*this); K+=i; return K;}
      Kabbala operator+(const double& d)  const {Kabbala K(*this); K+=d; return K;}

      Kabbala operator-();
      Kabbala operator-(const Kabbala& k) const {Kabbala K(*this); K-=k; return K;}
      Kabbala operator-(const Complex& c) const {Kabbala K(*this); K-=c; return K;}
      Kabbala operator-(const int& i)     const {Kabbala K(*this); K-=i; return K;}
      Kabbala operator-(const double& d)  const {Kabbala K(*this); K-=d; return K;}

      Kabbala operator*(const Kabbala& k) const {Kabbala K(*this); K*=k; return K;}
      Kabbala operator*(const Complex& c) const {Kabbala K(*this); K*=c; return K;}
      Kabbala operator*(const int& i)     const {Kabbala K(*this); K*=i; return K;}
      Kabbala operator*(const double& d)  const {Kabbala K(*this); K*=d; return K;}

      Kabbala operator/(const Kabbala& k) const {Kabbala K(*this); K/=k; return K;}
      Kabbala operator/(const Complex& c) const {Kabbala K(*this); K/=c; return K;}
      Kabbala operator/(const int& i)     const {Kabbala K(*this); K/=i; return K;}
      Kabbala operator/(const double& d)  const {Kabbala K(*this); K/=d; return K;}
  };
    
  // Easy Output
  inline std::ostream &operator<<(std::ostream& s, const Kabbala& k){s << k.String() << " with current value: " << k.Value(); return s;} 
  
  // this is deprecated, rather use operator ==
  inline bool IsZero(const Kabbala& a) {return ATOOLS::IsZero(a.Value());}
  inline int IsEqual(const Kabbala& a,const Complex& b) {return ATOOLS::IsEqual(a.Value(),b);}

  /*
  Some Calculations that are useful
  */
  Kabbala exp(const Kabbala& k1);
  Kabbala sin(const Kabbala& k1);
  Kabbala cos(const Kabbala& k1);
  Kabbala sqrt(const Kabbala& k1);
  Kabbala pow(const Kabbala& k1, const Complex& c);
  inline Kabbala sqr(const Kabbala& k1) {return pow(k1, Complex(2.));}
  Kabbala complexconjugate(const Kabbala& k1);
  //Kabbala abs(const Kabbala& k1);

  /*
  Some tools for comparison
  */
         bool operator==(const Kabbala& k1, const Kabbala& k2);
  inline bool operator==(const Kabbala& k1, const Complex& c) {return k1.Value() == c;}
  inline bool operator!=(const Kabbala& k1, const Kabbala& k2) {return !(k1==k2);}
  inline bool operator!=(const Kabbala& k1, const Complex& c) {return k1.Value() != c;}
         bool operator<(const Kabbala& k1, const Kabbala& k2);
  inline bool operator<(const Kabbala& k1, const Complex& c) {return abs(k1.Value()) < abs(c);}

  /*
  if the operands are the other way around
  */
  inline Kabbala operator+(const Complex& c,  const Kabbala& k) {return k+c;}
  inline Kabbala operator+(const int& i,      const Kabbala& k) {return k+i;}
  inline Kabbala operator+(const double& d,   const Kabbala& k) {return k+d;}

  inline Kabbala operator-(const Complex& c,  const Kabbala& k) {return -(k-c);}
  inline Kabbala operator-(const int& i,      const Kabbala& k) {return -(k-i);}
  inline Kabbala operator-(const double& d,   const Kabbala& k) {return -(k-d);}

  inline Kabbala operator*(const Complex& c,  const Kabbala& k) {return k*c;}
  inline Kabbala operator*(const int& i,      const Kabbala& k) {return k*i;}
  inline Kabbala operator*(const double& d,   const Kabbala& k) {return k*d;}

  Kabbala operator/(const Complex& c,  const Kabbala& k);
  inline Kabbala operator/(const int& i,      const Kabbala& k) {return Complex(i)/k;}
  inline Kabbala operator/(const double& d,   const Kabbala& k) {return Complex(d)/k;}

  /*
  some conversion; TODO only here because model not written properly yet
  */
  inline double ToDouble(const Kabbala& k) {return k.Value().real();}
  inline std::string ToString(const Kabbala& k) {return k.String();}
}

#endif
