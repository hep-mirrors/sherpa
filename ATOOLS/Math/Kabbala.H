#ifndef ATOOLS_Math_Kabbala_H
#define ATOOLS_Math_Kabbala_H

#include "ATOOLS/Math/MyComplex.H"
#include "ATOOLS/Math/MathTools.H" 
#include "ATOOLS/Org/MyStrStream.H"
#include <functional>
#include <map>

#define C_ZERO Complex(0., 0.)

namespace ATOOLS {
  class Kabbala {
    public:
      using Function_Argument = std::map<std::string, double_t>*;
      using Value_T = Complex;
      using Func = std::function<Value_T(Function_Argument)>;
    protected:
      std::string  shem;
      Value_T rishbon;
      Func lambda = [](Function_Argument map) {return C_ZERO;}; //needs initial value to prevent bad function call
    public:
      //Default constructor
      Kabbala(){
        shem = "";
        rishbon = C_ZERO;
        lambda = [](Function_Argument map) {return C_ZERO;};
      }
      // constructor without using initial value, that will be calculated from the function
      Kabbala(const std::string str, const Func func, Function_Argument map){
        shem = str;
        lambda = func;
        Update(map);
      }
      // legacy constructor, reverts lambda to basic look up of the given String
      Kabbala(const std::string str ,const Value_T C) {
        shem    = str;
        rishbon = C;
        msg_Debugging() << "No proper function set for the Kabbala. Choosing Basic Lookup." << std::endl;
        lambda = BasicLookUpFunction();
      }
      // constructor with lambda and initial value, these may be different values, the value is NOT updated
      Kabbala(const std::string str, const Value_T C, const Func func){
        shem = str;
        lambda = func;
        rishbon = C;
      }
      // copy constructor, this may be the reason the lambdas need tobe copied explicitly
      Kabbala(const Kabbala& k) {
        shem    = k.String();
        rishbon = k.Value();
        lambda = k.Lambda();
      }

      // Zero Kabbala
      static Kabbala Zero(){
        return Kabbala(std::string("0"), C_ZERO, [](Function_Argument map) {return C_ZERO;});
      }

      // in-place negation of the kabbala
      void Negate(){
        rishbon = -rishbon;
        shem = std::string("-(")+shem+std::string(")");
        Func copy(lambda);
        lambda = [copy](Function_Argument m) {return -copy(m);};
      }

      // update the value using the lambda and a map of external constants
      Value_T Update(Function_Argument map){
        return rishbon = lambda(map);
      }

      // look for occurence of the parameter in the string
      bool DependsOn(std::string param){
        return shem.find(param) != std::string::npos;
      }

      // Basis of recursive lambdas, uses copy of const string for lookup
      Func BasicLookUpFunction(){
        std::string copy(shem);
        return [=](Function_Argument map) {
          if (map->count(copy) == 0) return C_ZERO; 
          return Value_T(map->at(copy), .0); 
        };
      }

      Kabbala& operator=(const Kabbala& k) {
        if (this!=&k) {
          shem    = k.String();
          rishbon = k.Value();
          lambda = k.Lambda();
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
      Kabbala& operator*=(const Value_T& c);
      Kabbala& operator/=(const Kabbala& k);
      
      const Value_T &Value() const {return rishbon;} 
      const std::string &String() const {return shem;}
      const Func &Lambda() const {return lambda;}
    
      void SetValue(Value_T val) {rishbon = val;}
      void SetString(std::string str) {shem = str;} 
      void SetLambda(Func func) {lambda = func;}

      friend int IsZero(const Kabbala& a);
      friend int IsEqual(const Kabbala& a,const Value_T& b);
  };
    
    inline std::ostream &operator<<(std::ostream& s, const Kabbala& k)
    {
      s << k.String() << " with current value: " << k.Value();
      return s;
    } 
    
    inline int IsZero(const Kabbala& a) {
      return ATOOLS::IsZero(a.Value());
    }

    inline int IsEqual(const Kabbala& a,const Kabbala::Value_T& b) {
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
      // why is this different
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
    
      inline Kabbala operator* (const Kabbala::Value_T c1, const Kabbala& k2) 
    {
      Kabbala k(k2); 
      k *= c1;
      return k;
    }
    
    inline Kabbala operator* (const Kabbala& k1,const Kabbala::Value_T c2) 
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
      Kabbala::Func copy(k.Lambda());
      k.SetLambda([copy](Kabbala::Function_Argument m){return std::exp(copy(m));});
      return k;
    }

    inline Kabbala sin(const Kabbala& k1)
    {
      // why
      Kabbala k(k1);
      return k;
    }

}

#endif


