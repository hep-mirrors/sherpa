#ifndef Kabbala_h
#define Kabbala_h

#include "MyString.H"
#include "MyComplex.H"
#include "MathTools.H" 
#include "MyStrStream.H"

namespace AMEGIC {
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
    Kabbala& operator+=(const Kabbala& k) {
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
    
    Kabbala& operator-=(const Kabbala& k) {
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
    
    Kabbala operator-() {
      return Kabbala(std::string("-")+std::string("(")+shem+std::string(")"),-rishbon);
    }

    Kabbala operator+() {
      return Kabbala(shem,rishbon);
    }
    
    Kabbala& operator*=(const Kabbala& k) {
      rishbon *= k.Value();
      std::string save = shem;
      shem  = std::string("(") + save + std::string(")*(");
      shem += k.String();
      shem += std::string(")");
      return *this;
    }
    
    Kabbala& operator*=(const int& i) {
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

    Kabbala& operator*=(const double& d) {
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

    Kabbala& operator/=(const Kabbala& k) {
      rishbon /= k.Value();
      std::string save = shem;
      shem  = std::string("(") + save + std::string(")/(");
      shem += k.String();
      shem += std::string(")");
      return *this;
    }
    
    Complex Value()  const {return rishbon;} 
    std::string  String() const {return shem;}
   
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


