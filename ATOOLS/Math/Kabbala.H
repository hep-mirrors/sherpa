#ifndef ATOOLS_Math_Kabbala_H
#define ATOOLS_Math_Kabbala_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H" 
#include "ATOOLS/Org/MyStrStream.H"

namespace ATOOLS {
  class Kabbala {
    std::string  shem;
    Complex rishbon;
  public:
    Kabbala(const std::string str = std::string(""),const Complex C = Complex(0.,0.)) {
      shem    = str;
      rishbon = C;
    } 
    Kabbala(const Kabbala& k) {
      shem    = k.String();
      rishbon = k.Value();
    }
    Kabbala& operator=(const Kabbala& k) {
      if (this!=&k) {
	shem    = k.String();
	rishbon = k.Value();
      } 
      return *this;
    }
    Kabbala& operator+=(const Kabbala& k);
    Kabbala& operator-=(const Kabbala& k);
    
    Kabbala operator-();
    Kabbala operator+();

    Kabbala& operator*=(const int& i);
    Kabbala& operator*=(const Kabbala& k);
    Kabbala& operator*=(const double& d);
    Kabbala& operator*=(const Complex& c);
    Kabbala& operator/=(const Kabbala& k);
    
    const Complex &Value() const {return rishbon;} 
    const std::string &String() const {return shem;}
   
    void SetValue(Complex val) {rishbon = val;}
    void SetString(std::string str) {shem    = str;} 

    friend int IsZero(const Kabbala& a);
    friend int IsEqual(const Kabbala& a,const Complex& b);
  };
  
  inline int IsZero(const Kabbala& a) {
    return ATOOLS::IsZero(a.Value());
  }

  inline int IsEqual(const Kabbala& a,const Complex& b) {
    return ATOOLS::IsEqual(a.Value(),b);
  }

  inline Kabbala operator+(const Kabbala& k1,const Kabbala& k2)
  {
    Kabbala k(k1);
    k += k2;
    return k;
  }

  inline Kabbala operator-(const Kabbala& k1,const Kabbala& k2)
  {
    Kabbala k(k1);
    k -= k2;
    return k;
  }

  inline Kabbala operator* (const Kabbala& k1, const Kabbala& k2) 
  {
    Kabbala k(k1);
    k *= k2;
    return k;
  }

  inline Kabbala operator* (const int& i1, const Kabbala& k2) 
  {
    Kabbala k(k2); 
    k *= i1;
    return k;
  }

  inline Kabbala operator* (const Kabbala& k1,const int& i2) 
  {
    return i2*k1;
  }

  
  inline Kabbala operator* (const double d1, const Kabbala& k2) 
  {
    Kabbala k(k2); 
    k *= d1;
    return k;
  }
  
  inline Kabbala operator* (const Kabbala& k1,const double d2) 
  {
    return d2*k1;
  }
  
    inline Kabbala operator* (const Complex c1, const Kabbala& k2) 
  {
    Kabbala k(k2); 
    k *= c1;
    return k;
  }
  
  inline Kabbala operator* (const Kabbala& k1,const Complex c2) 
  {
    Kabbala k(k1); 
    k *= c2;
    return k;
  }
  
  inline Kabbala operator/ (const Kabbala& k1, const Kabbala& k2) 
  {
    Kabbala k(k1);
    k /= k2;
    return k;
  }

  inline bool operator==(const Kabbala& k1,const Kabbala& k2) 
  {
    if (k1.Value()!=k2.Value())   return false;
    if (k1.String()!=k2.String()) return false;
    return true;
  }

  inline bool operator!=(const Kabbala& k1,const Kabbala& k2)
  {
    return !(k1==k2);
  }

  inline Kabbala exp(const Kabbala& k1)
  {
    Kabbala k(k1);

    k.SetValue(std::exp(k.Value()));
    k.SetString(std::string("exp(") + k.String() + std::string(")"));

    return k;
  }

  inline Kabbala sin(const Kabbala& k1)
  {
    Kabbala k(k1);
    return k;
  }
}

#endif


