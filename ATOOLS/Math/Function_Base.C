#include "ATOOLS/Math/Function_Base.H"

#include "ATOOLS/Org/My_Limits.H"
#include "ATOOLS/Math/Algebra_Interpreter.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Org/Message.H"
#include <math.h>

using namespace ATOOLS;

Function_Base::~Function_Base() {}

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

/* (C) Copr. 1986-92 Numerical Recipes Software */

double Function_Base::WDBSolve
(const double &y,const double &xmin,const double &xmax,
 const double &precision,const int maxit)
{
  double eps=std::numeric_limits<double>::epsilon();
  double a=xmin, b=xmax, c=b;
  double fa=(*this)(a)-y, fb=(*this)(b)-y, fc=fb;
  double d=0.0, e=0.0, p=0.0, q=0.0, r=0.0, s=0.0, x=0.0, t=0.0;
  for (int it=0;it<maxit;++it) {
    if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (dabs(fc)<dabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    t=2.0*eps*dabs(b)+0.5*precision;
    x=0.5*(c-b);
    if (dabs(x)<=t || fb==0.0) return b;
    if (dabs(e)>=t && dabs(fa)>dabs(fb)) {
      s=fb/fa;
      if (a==c) {
	p=2.0*x*s;
	q=1.0-s;
      }
      else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*x*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p>0.0) q=-q;
      p=dabs(p);
      double min1=3.0*x*q-dabs(t*q);
      double min2=dabs(e*q);
      if (2.0*p<(min1<min2?min1:min2)) {
	e=d;
	d=p/q;
      }
      else {
	d=x;
	e=d;
      }
    }
    else {
      d=x;
      e=d;
    }
    a=b;
    fa=fb;
    if (dabs(d)>t) b+=d;
    else b+=x>=0.0?dabs(t):-dabs(t);
    fb=(*this)(b)-y;
  }
  msg_Error()<<METHOD<<"(): No solution for f(x) = "
	     <<y<<" in ["<<xmin<<","<<xmax<<"]"<<std::endl;
  return 0.5*(xmin+xmax);
}

double Function_Base::FindZero(double x_min, double x_max,int MAX_ITR,double precision)
{
  return WDBSolve(0.0,x_min,x_max,precision,MAX_ITR);
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

  class GMean_Function_Wrapper: public Function {
  private:

    Function_Base *p_f;

  public:

    inline GMean_Function_Wrapper(Function_Base *const f):
      Function("GMean_"+f->Name()), p_f(f) {}

    Term *Evaluate(const std::vector<Term*> &args) const
    {
      msg_Debugging()<<"GMean_"<<p_f->Name()<<"(): {\n";
      double ym(1.0), xm(1.0);
      double xmax(std::numeric_limits<double>::max()), xmin(-xmax);
      for (size_t i(0);i<args.size();++i) {
	double cx(args[i]->Get<double>()), cy((*p_f)(cx));
	msg_Debugging()<<"  x_{"<<i<<"} = "<<cx
		       <<", y_{"<<i<<"} = "<<cy<<"\n";
	xm*=cx;
	ym*=cy;
	if (cx<xmax) xmax=cx;
	if (cx>xmin) xmin=cx;
      }
      ym=pow(ym,1.0/args.size());
      xm=p_f->WDBSolve(ym,xmin,xmax);
      if (!IsEqual((*p_f)(xm),ym)) msg_Error()<<"GMean_"<<
	p_f->Name()<<"(): Could not solve for x."<<std::endl; 
      msg_Debugging()<<"} -> y = "<<ym<<" -> x = "<<xm<<"\n";
      Term *res(Term::New(xm));
      p_interpreter->AddTerm(res);
      return res;
    }

  };// end of class GMean_Function_Wrapper

}// end of namespace ATOOLS

Function *Function_Base::GetAIFunction()
{
  return new Function_Wrapper(this);
}

Function *Function_Base::GetAIGMeanFunction()
{
  return new GMean_Function_Wrapper(this);
}

