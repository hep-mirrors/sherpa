#ifndef ATOOLS_Math_Special_Functions_H
#define ATOOLS_Math_Special_Functions_H

#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include <limits>

// TOOD check this file for correct implementations on this branch
namespace ATOOLS {
  /*!
    \file
    This file contains an assortment of special functions

    \todo We should use the C++23 libraries there and match e.g. the
          various Bessel functions or the exponential integral here.
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
	//return -std::expint(-x);
	incgam = -0.577215664902-log(x);
	double prev = 0., fac(1.);
	size_t i(1);
	msg_Debugging()<<METHOD<<"("<<a<<", "<<x<<"), inc = "<<incgam<<"\n";
	while (dabs(1.-incgam/prev)>1.e-12) {
	  prev    = incgam;
	  incgam += fac*pow(x,i);
	  fac    *= -double(i)/sqr(double(i+1));
	  i++;
	}
        msg_Debugging()<<"---> inc = "<<incgam<<"\n";
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
      switch (order) {
      case 0: return J0(arg);
      case 1: return J1(arg);
      default:
	THROW(critical_error,"J("+ToString(order)+") not implemented yet.");
      }
    }

    double J0(const double & arg) const
    {
      double jn(0.),darg(dabs(arg));
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
      return jn;
    }

    double J1(const double & arg) const
    {
      if (dabs(arg)<=1.e-12) return 0.;
      double x2 = arg*arg;
      if (dabs(arg)>=100.) {
	return (sqrt(2./(M_PI*arg)) *
		(1.    + 15./(128.*x2))   * cos(arg-3.*M_PI/4.) +
		(3./8. - 105./(1024.*x2)) * cos(arg-3.*M_PI/4.)/arg);
      }
      // using expansion from Maass & Martin, Results in Physics 8 (2018) 1234
      // https://www.sciencedirect.com/science/article/pii/S2211379718300111
      double l  = 0.1, l2 = l*l, l2term = sqrt(1.+l2*x2);
      double p0 = 1.776322448,   p1 = 0.2250803518;
      double P0 = -0.7763224930, P1 = -0.03147133771;
      double q1 = 0.4120981204,  q2 = 0.006571619275;
      double p2 = 2.*sqrt(l/M_PI)*q2, P2 = -l*p2;
      return ( 1./(2.*sqrt(l2term)) *
	       ( (p0+p1*x2+p2*x2*x2)/(1.+q1*x2+q2*x2*x2)*sin(arg) +
		 arg/l2term * (P0+P1*x2+P2*x2*x2)/(1.+q1*x2+q2*x2*x2)*cos(arg) ) );
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

    double Kn(const int order,const double & x) const
    {
      if (x<=0.) return 0.;
      double kn(0.),darg(dabs(x)), y;
      if (order==0) {
        // The implementation is based on the formulas from the "Handbook of Mathematical
        // Functions" by Abramowitz and Stegun (A&S).
        if (x <= 0.0) {
          throw std::domain_error("Argument x for bessel_K0 must be positive.");
        }

        // --- Crossover point between the two approximation regimes ---
        const double X_CROSSOVER = 2.0;

        if (x <= X_CROSSOVER) {
          // --- Series Expansion for small x (0 < x <= 2.0) ---
          // This regime handles the logarithmic singularity near x=0.
          // A&S 9.6.11: K₀(x) = -[ln(x/2) + γ]I₀(x) + Σ...
          // We compute I₀(x) and the summation part simultaneously.

          const double EULER_GAMMA = 0.57721566490153286060; // γ

          double term = 1.0;      // First term of I₀(x) series
          double i0_sum = 1.0;    // Sum for I₀(x)
          double k0_sum = 0.0;    // Sum for the second part of K₀(x)
          double psi_k = 0.0;     // Digamma approximation: H_k - γ, where H_k is harmonic number
          double k_fact = 1.0;
          double x_pow_2k = 1.0;

          const double x_half_sq = 0.25 * x * x;
          const double log_x_half = log(0.5 * x);

          for (int k = 1; k <= 15; ++k) { // 15 terms are sufficient for double precision up to x=2
            // Update terms for the k-th iteration
            psi_k += 1.0 / k;
            k_fact *= k;
            x_pow_2k *= x_half_sq;

            term = x_pow_2k / (k_fact * k_fact);

            // Accumulate sums
            i0_sum += term;
            k0_sum += term * psi_k;

            // Stop if the term becomes negligible to avoid unnecessary computation
            if (term < std::numeric_limits<double>::epsilon() * i0_sum) {
              break;
            }
          }

          return k0_sum - (log_x_half + EULER_GAMMA) * i0_sum;

        } else {
          // --- Asymptotic Expansion for large x (x > 2.0) ---
          // A&S 9.7.2: K₀(x) ~ sqrt(π/(2x)) * e^(-x) * [ 1 - 1/(8x) + 9/(128x^2) - ... ]

          const double inv_x = 1.0 / x;
          const double inv_x_sq = inv_x * inv_x;

          double series = 1.0;
          series -= 0.125 * inv_x;                   // -1/(8x)
          series += 0.0703125 * inv_x_sq;            // +9/(128x^2)
          series -= 0.0732421875 * inv_x_sq * inv_x; // -75/(1024x^3)

          return exp(-x) * sqrt(M_PI / (2.0 * x)) * series;
        }
      }
      else if (order==1) {
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
        msg_Error()<<"Error in "<<METHOD<<"("<<order<<", "<<x<<")\n:"
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
      // Test BesselK0
      std::cout << "------------------------------------------------------------------" << std::endl;
      std::cout << std::left << std::setw(12) << "x"
                << std::setw(22) << "bessel_K0(x)"
                << std::setw(22) << "Reference"
                << std::setw(15) << "Rel. Error" << std::endl;
      std::cout << "------------------------------------------------------------------" << std::endl;

      // Test cases: struct containing {x_value, reference_K0(x_value)}
      // Reference values obtained from WolframAlpha.
      const std::vector<std::pair<double, double>> test_cases = {
              {0.01, 4.721244730161094},
              {0.1,  2.427069024702016},
              {0.5,  0.924419071227666},
              {1.0,  0.42102443824070834},
              {1.5,  0.2138055626475258},
              {2.0,  0.1138938727495334}, // The crossover point
              {2.1,  0.10078374088996696},
              {5.0,  0.0036910983340425947},
              {10.0, 0.00001778006231616765},
              {20.0, 5.741237815336525e-10},
              {50.0, 3.4101677497894956e-23},
              {100.0, 4.656628229175903e-45}
      };

      bool all_tests_passed = true;
      for (const auto& test : test_cases) {
        double x = test.first;
        double ref_val = test.second;
        double calc_val = Kn(0, x);
        double rel_error = (ref_val != 0) ? std::abs((calc_val - ref_val) / ref_val) : std::abs(calc_val);

        std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(2) << x
                  << std::scientific << std::setprecision(12)
                  << std::setw(22) << calc_val
                  << std::setw(22) << ref_val
                  << std::setw(15) << rel_error;

        if (rel_error > 1e-6) {
          all_tests_passed = false;
          std::cout << "  <-- FAILED";
        }
        std::cout << std::endl;
      }
      std::cout << "------------------------------------------------------------------" << std::endl;
      // Test edge cases
      try {
        std::cout << "\nTesting edge case x=0..." << std::endl;
        Kn(0, 0.0);
      } catch (const std::domain_error& e) {
        std::cout << "Caught expected exception: " << e.what() << std::endl;
      }
      std::cout << "\nTest summary: " << (all_tests_passed ? "All tests passed." : "Some tests FAILED.") << std::endl;

      // Test BesselK1
      std::cout << "------------------------------------------------------------------" << std::endl;
      std::cout << std::left << std::setw(12) << "x"
                << std::setw(22) << "bessel_K1(x)"
                << std::setw(22) << "Reference"
                << std::setw(15) << "Rel. Error" << std::endl;
      std::cout << "------------------------------------------------------------------" << std::endl;

      // Test cases: struct containing {x_value, reference_K1(x_value)}
      // Reference values obtained from WolframAlpha or trusted high-precision sources.
      const std::vector<std::pair<double, double>> test_cases_K1 = {
              {0.01, 99.97389411829624},
              {0.1, 9.853844780870606},
              {0.5, 1.656441120003301},
              {1.0, 0.6019072301972347},
              {1.5, 0.27738780045684386},
              {2.0, 0.13986588181652243},
              {2.1, 0.12274641153350789},
              {5.0, 0.004044613445452164},
              {10.0, 0.00001864877345382558},
              {20.0, 5.883057969557039e-10},
              {50.0, 3.4441022267175555e-23},
              {100.0, 4.67985373563691e-45}
      };

      bool all_tests_passed_K1 = true;
      for (const auto& test : test_cases_K1) {
        double x = test.first;
        double ref_val = test.second;
        double calc_val = Kn(1, x);
        double rel_error = (ref_val != 0) ? std::abs((calc_val - ref_val) / ref_val) : std::abs(calc_val);

        std::cout << std::left << std::setw(12) << std::fixed << std::setprecision(2) << x
                  << std::scientific << std::setprecision(12)
                  << std::setw(22) << calc_val
                  << std::setw(22) << ref_val
                  << std::setw(15) << rel_error;

        if (rel_error > 1e-6) {
          all_tests_passed_K1 = false;
          std::cout << "  <-- FAILED";
        }
        std::cout << std::endl;
      }
      std::cout << "------------------------------------------------------------------" << std::endl;

      // Test edge cases
      try {
        std::cout << "\nTesting edge case x=0..." << std::endl;
        Kn(1, 0.0);
      } catch (const std::domain_error& e) {
        std::cout << "Caught expected exception: " << e.what() << std::endl;
      }

      std::cout << "\nTest summary: " << (all_tests_passed_K1 ? "All tests passed." : "Some tests FAILED.") << std::endl;
    }
  };
  extern Special_Functions SF;

  class ExpInt : public ATOOLS::Function_Base {
  public:
    double GetValue(double T) override { return -exp(-T)/T; }
    double operator()(double T) override { return GetValue(T); }
    double GetExpInt(double X){
      if(X>0.) exit(1);
      ATOOLS::Gauss_Integrator integrator(this);
      return integrator.Integrate(-X,100000.,1.e-2,1);
    }
  };

}
#endif
