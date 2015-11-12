#ifndef SHRIMPS_Tools_Special_Functions_H
#define SHRIMPS_Tools_Special_Functions_H

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Gauss_Integrator.H"

using namespace ATOOLS;

namespace SHRIMPS {
  /*!
    \file This file contains an assortment of special functions
    neccessary for the SHRIMPS::Form_Factor.  
    
    \todo Further expand and move to ATOOLS::MathTools
  */
  class Special_Functions {
  public:
    double LnGamma(const double & x) const
    {
      static const double coeff[6] = 
	{ 76.18009172947146, -86.50532032941677, 24.01409824083091, 
	  -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5 };
      double y(x), tmp(x+5.5);
      tmp -= (x+0.5)*log(tmp);
      double arg(1.000000000190015);
      for (short int j=0;j<6;j++) arg += coeff[j]/++y;
      return -tmp+log(2.5066282746310005*arg/x);
    }

    double IncompleteGamma(const double & a, const double & x) const
    {
      double incgam(0.);
      if (x<0. || a<0.) {
	msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<x<<"):\n"
		   <<"   Out of bounds. "
		   <<"Will return 0 and hope for the best.\n";
	return incgam;
      }
      if (a==0.) {
	incgam = -0.577215664902-log(x);
	double pref(1.);
	for (int i=1;i<21;i++) {
	  incgam += pref*pow(x,i);
	  pref *= -double(i)/sqr(double(i+1));
	}
	return incgam;
      }
      double lngamma(LnGamma(a));
      if (x<a+1.) {
	if (x==0.) return incgam;
	else {
	  double aprime(a), del(1./a), sum(1./a);
	  for (short int i=0;i<100;i++) {
	    ++aprime;
	    del *= x/aprime;
	    sum += del;
	    if (dabs(del)<dabs(sum)*1.e-12) {
	      incgam = 1.-sum*exp(-x+a*log(x)-lngamma);
	      break;
	    }
	  }
	}
      }
      else {
	msg_Error()<<"Error in "<<METHOD<<"("<<a<<", "<<x<<") :\n"
		   <<"   Not implemented yet. "
		   <<"Will return 0 and hope for the best.\n";
      }
      return incgam;
    }

    double Jn(const int order,const double & arg) const 
    {
      double jn(0.),darg(dabs(arg));
      if (order==0) {
	if (darg<=1.e-12) return 1.;
	if (darg<8.) {
	  double x2(darg*darg);
	  double num = 
	    ((((-184.9052456*x2+77392.33017)*x2-
	       11214424.18)*x2+651619640.7)*x2-13362590354.)*x2+57568490574.;
	  double den =
	    ((((x2+267.8532712)*x2+59272.64853)*
	      x2+9494680.718)*x2+1029532985)*x2+57568490411.;
	  jn = num/den;
	}
	else {
	  double x2(64.0/sqr(darg)), xx(darg-0.785398164);
	  double pref1 = 
	    (((0.2993887211e-6*x2-0.2073370639e-5)*x2+
	      0.2734510407e-4)*x2-0.1098628627e-2)*x2+1.;
	  double pref2 =
	    (((-0.934945152e-7*x2+0.7621095161e-6)*x2-
	      0.6911147651e-5)*x2+0.1430488765e-3)*x2-0.1562499995e-1;
	  jn = sqrt(0.636619772/darg)*(cos(xx)*pref1-8./darg*sin(xx)*pref2);
	}
      }
      else {
	msg_Error()<<"Error in "<<METHOD<<"("<<order<<", "<<arg<<"):\n"
		   <<"   Not implemented yet.  Exit the run.\n";
	exit(1);
      }
      return jn;
    }

    double In(const int order,const double & arg) const 
    {
      double in(0.),darg(dabs(arg)), y;
      switch (order) {
      case 1:
	if (darg<3.75) {
	  y  = sqr(arg/3.75);
	  in = darg*((((((0.32411e-3*y+0.301532e-2)*y+0.2658733e-1)*y+
			0.15084934)*y+0.51498869)*y+0.87890594)*y+0.5);
	}
	else {
	  y  = 3.75/darg;
	  in = ((((((((-0.420059e-2*y+0.1787654e-1)*y-0.2895312e-1)*y+
		     0.2282967)*y-0.1031555e-1)*y+0.163801e-2)*y-
		  0.362018e-2)*y-0.3988024e-1)*y+0.39894228) * 
	    exp(darg)/sqrt(darg);
	}
	break;
      case 0:
	if (darg<3.75) {
	  y  = sqr(arg/3.75);
	  in = (((((0.45813e-2*y+0.360768e-1)*y+0.2659723)*y+1.2067492)*y
		 +3.0899424)*y+3.5156229)*y+1.0;
	}
	else {
	  y  = 3.75/darg;
	  in = ((((((((0.392377e-2*y-0.1647633e-1)*y+0.2635537e-1)*y-
		     0.2057706)*y+0.916281e-2)*y-0.157565e-2)*y+
		  0.225319e-2)*y+0.1328592e-1)*y+0.39894228) * 
	    exp(darg)/sqrt(darg);
	}    
	break;
      default:
	msg_Error()<<"Error in "<<METHOD<<"("<<order<<", "<<arg<<"):\n"
		   <<"   Not implemented yet.  Exit the run.\n";
	exit(1);
      }
      return in;
    }


    double I0(const double & arg)  const { return In(0,arg); }

    double Kn(const int order,const double & arg) const
    {
      if (arg<=0.) return 0.;
      double kn(0.),darg(dabs(arg)), y;
      if (order==1) {
	if (darg<=2.0) {
	  y   = sqr(darg)/4.;
	  kn  = log(darg/2.)*In(1,darg);
	  kn += ((((((-0.4686e-4*y-0.110404e-2)*y-0.1919402e-1)*y-0.18156897)*y
		   -0.67278579)*y+0.15443144)*y+1.)/darg;
	}
	else {
	  y   = 2./darg;
	  kn  = exp(-darg)/sqrt(darg);
	  kn *= ((((((-0.68245e-3*y+0.325614e-2)*y-0.780353e-2)*y+
		    0.1504268e-1)*y-0.3655620e-1)*y+0.23498619)*y+1.25331414);
	}
      }
      else {
	msg_Error()<<"Error in "<<METHOD<<"("<<order<<", "<<arg<<")\n:"
		   <<"   Not implemented yet.  Exit the run.\n";
	exit(1);
      }
      return kn;
    }

    void TestBessel(const std::string & dirname) const
    {
      double arg(0.);
      std::ofstream was;
      std::string filename = dirname+std::string("/BesselJ0I0I1K1.dat");
      was.open(filename.c_str());
      long int cnt(0);
      while (arg<10.) {
	was<<arg<<"   "<<Jn(0,arg)<<"   "
	   <<In(0,arg)<<"   "<<In(1,arg)<<"   "<<Kn(1,arg)<<std::endl;
	//if (cnt==0 || !(cnt%500)) 
	//  msg_Out()<<"  J_0("<<arg<<") = "<<Jn(0,arg)
	//	   <<"  K_1("<<arg<<") = "<<Kn(1,arg)<<std::endl;
	arg += 0.001;
	cnt++;
      }
      was.close();
    }
  };

  extern Special_Functions SF;


  class ExpInt : public ATOOLS::Function_Base {
  public:
    virtual double GetValue(double T) { return -exp(-T)/T; }
    virtual double operator()(double T) { return GetValue(T); }
    virtual double GetExpInt(double X){
	if(X>0.) exit(1);
  	ATOOLS::Gauss_Integrator integrator(this);
  	return integrator.Integrate(-X,100000.,1.e-2,1);
        }
  };

}
#endif
