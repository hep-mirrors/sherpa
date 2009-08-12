#include "ATOOLS/Math/Function_Base.H"

#include "ATOOLS/Math/Algebra_Interpreter.H"
#include <math.h>

using namespace ATOOLS;

Function_Base::~Function_Base() {}

void Function_Base::SetDefault(double _m_defval)    
{ m_defval=_m_defval; }

void Function_Base::SetType(std::string _m_type) 
{ m_type=_m_type; }

void Function_Base::SetParameters(double *parameters)    
{ return; }

double Function_Base::GetValue(double x)        
{ return (*this)(x); }   

double Function_Base::operator()(double x)         
{ return m_defval; }

double Function_Base::operator()()               
{ return m_defval; }

std::string Function_Base::Type()                     
{ return m_type; }


double Function_Base::FindZero(double x_min, double x_max,int MAX_ITR,double precision) {

  double root = 0.5 * (x_min + x_max);
  double x_lower = x_min;
  double x_upper = x_max;

  bool SUCCESS = false;

  for(int i=0;i<MAX_ITR;i++){
    SUCCESS=this->IterateBisection(root,x_lower,x_upper,precision);
    if(SUCCESS) break;
  }

  return root;

}


bool Function_Base::IterateBisection(double& root, double& x_lower, double& x_upper, double precision){
  // Bisection algorithm inspired by gsl-1.9

  double x_bisect,f_bisect;

  const double f_lower = this->GetValue(x_lower);
  const double f_upper = this->GetValue(x_upper);


  if ( fabs(f_lower-0.0)<precision )
    {
      root = x_lower;
      x_upper = x_lower;
      return true;
    }

  if ( fabs(f_upper-0.0)<precision)
    {
      root = x_upper ;
      x_lower = x_upper;
      return true;
    }

  x_bisect = ( x_lower + x_upper) / 2.0;
  f_bisect = this->GetValue(x_bisect);

  if ((f_lower > 0.0 && f_bisect < 0.0) || (f_lower < 0.0 && f_bisect > 0.0))
    {
      root = 0.5 * (x_lower + x_bisect) ;
      x_upper = x_bisect;
    }
  else
    {
      root = 0.5 * (x_bisect + x_upper) ;
      x_lower = x_bisect;
    }

  return false; // keep searching!!!

}

namespace ATOOLS {

  class Function_Wrapper: public Function {
  private:

    Function_Base *p_f;

  public:

    inline Function_Wrapper(Function_Base *const f):
      Function(f->Name()), p_f(f) {}

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      Term *res(Term::New((*p_f)(args[0]->Get<double>())));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class Function_Wrapper

}// end of namespace ATOOLS

Function *Function_Base::GetAIFunction()
{
  return new Function_Wrapper(this);
}

