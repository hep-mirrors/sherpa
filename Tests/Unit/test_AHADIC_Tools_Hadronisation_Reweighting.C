#include "AHADIC++/Tools/Hadronisation_Reweighting.H"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// The quadrature under test, a nested class of the reweighting.
typedef AHADIC::Hadronisation_Reweighting::Frag_Norm Frag_Norm;

// Independent oracle for
//   log int_zmin^zmax dz z^alpha (1-z)^beta exp(-c/z),
// deliberately not sharing any code with Frag_Norm: a brute-force composite
// 8-point Gauss-Legendre rule on 2*NH panels graded geometrically in z below
// the centre of the interval and in 1-z above it, accumulated with
// log-sum-exp. Converged to ~1e-11 for the integrands tested here.
namespace {

  double reference(const double zmin, const double zmax, const double alpha,
                   const double beta, const double c) {
    static const double x[8] = {
      -9.60289856497536176e-01, -7.96666477413626728e-01,
      -5.25532409916328991e-01, -1.83434642495649780e-01,
       1.83434642495649780e-01,  5.25532409916328991e-01,
       7.96666477413626728e-01,  9.60289856497536176e-01};
    static const double w[8] = {
      1.01228536290376689e-01, 2.22381034453374343e-01,
      3.13706645877887047e-01, 3.62683783378361768e-01,
      3.62683783378361768e-01, 3.13706645877887047e-01,
      2.22381034453374343e-01, 1.01228536290376689e-01};
    const int NH = 20000;
    const double zc = 0.5*(zmin+zmax);
    std::vector<double> edge(2*NH+1);
    edge[0] = zmin; edge[NH] = zc; edge[2*NH] = zmax;
    for (int k=1; k<NH; ++k) {
      edge[k]    = zmin*std::pow(zc/zmin, double(k)/NH);
      edge[NH+k] = 1.-(1.-zmax)*std::pow((1.-zc)/(1.-zmax), double(NH-k)/NH);
    }
    std::vector<double> e, ww;
    e.reserve(2*NH*8); ww.reserve(2*NH*8);
    double emax = -std::numeric_limits<double>::max();
    for (int ip=0; ip<2*NH; ++ip) {
      const double a = edge[ip], b = edge[ip+1];
      const double mid = 0.5*(a+b), half = 0.5*(b-a);
      if (half<=0.) continue;
      for (int i=0; i<8; ++i) {
        const double z = mid + half*x[i];
        const double v = alpha*std::log(z) + beta*std::log1p(-z) - c/z;
        e.push_back(v); ww.push_back(half*w[i]);
        if (v>emax) emax = v;
      }
    }
    double sum = 0.;
    for (size_t i=0; i<e.size(); ++i) sum += ww[i]*std::exp(e[i]-emax);
    return emax + std::log(sum);
  }

  // What the reweighting actually uses: the quadrature error that survives in
  // log(w) = log(f_var/N_var) - log(f_nom/N_nom). The f's are exact, so this
  // is the difference of the two normalisation errors.
  double log_weight_bias(const double zmin, const double zmax,
                         const double alpha, const double beta,
                         const double cnom, const double cvar,
                         const double cmin, const double cmax) {
    Frag_Norm fn;
    fn.SetRange(zmin,zmax,alpha,beta,cmin,cmax);
    const double enom = fn(alpha,beta,cnom) - reference(zmin,zmax,alpha,beta,cnom);
    const double evar = fn(alpha,beta,cvar) - reference(zmin,zmax,alpha,beta,cvar);
    return enom - evar;
  }
}

TEST_CASE("Frag_Norm reproduces a known integral", "[AHADIC][Frag_Norm]") {
  Frag_Norm fn;
  // c = 0 reduces to an incomplete Beta function; check against the oracle.
  fn.SetRange(1.e-3,0.999,2.5,0.13,0.,0.);
  CHECK_THAT(fn(2.5,0.13,0.),
             Catch::Matchers::WithinAbs(reference(1.e-3,0.999,2.5,0.13,0.),1.e-9));
  // int_0^1 z^1 (1-z)^0 = 1/2 over almost the full range
  fn.SetRange(1.e-9,1.-1.e-9,1.,0.,0.,0.);
  CHECK_THAT(std::exp(fn(1.,0.,0.)),
             Catch::Matchers::WithinRel(0.5,1.e-8));
}

TEST_CASE("Frag_Norm peak location", "[AHADIC][Frag_Norm]") {
  // Stationary point of alpha ln z + beta ln(1-z) - c/z.
  const double alpha = 2.5, beta = 0.13, c = 3.;
  const double zp = Frag_Norm::Peak(alpha,beta,c,1.e-4,0.9999);
  const double d  = 1.e-6;
  const auto f = [&](double z){
    return alpha*std::log(z)+beta*std::log1p(-z)-c/z; };
  CHECK(f(zp) > f(zp+d));
  CHECK(f(zp) > f(zp-d));
  // c = 0 gives the Beta mode alpha/(alpha+beta)
  CHECK_THAT(Frag_Norm::Peak(alpha,beta,0.,1.e-4,0.9999),
             Catch::Matchers::WithinRel(alpha/(alpha+beta),1.e-10));
}

// The regression this guards: a rule whose panel widths do not depend on c
// stops resolving exp(-c/z) once the peak, of width ~z*^2/c, falls inside a
// single panel. The error then lands undamped in log(w) and compounds over
// the splittings of an event.
TEST_CASE("Frag_Norm resolves strongly peaked integrands",
          "[AHADIC][Frag_Norm]") {
  // Narrow z window away from z=1 with a large exponent: the configuration
  // that used to be underestimated by many e-foldings.
  CHECK(std::abs(log_weight_bias(0.05,0.20,3.26,0.11,200.,800.,200.,800.))
        < 1.e-4);
  CHECK(std::abs(log_weight_bias(0.001,0.20,2.50,0.25,200.,800.,200.,800.))
        < 1.e-4);
  CHECK(std::abs(log_weight_bias(0.01,0.30,3.26,0.11,180.,721.,180.,721.))
        < 1.e-4);
  // ... and the same at the other end, where the peak sits against zmin.
  CHECK(std::abs(log_weight_bias(1.e-5,0.999,2.50,0.13,0.324,0.072,0.072,0.324))
        < 1.e-4);
}

TEST_CASE("Frag_Norm is accurate over the AHADIC variation ranges",
          "[AHADIC][Frag_Norm]") {
  // alpha, beta, gamma of the four cluster types at their defaults, scanned
  // over generous variation ranges of gamma and KT_0.
  struct Type { double alpha, beta, gamma; };
  const std::vector<Type> types = {
    {2.50,0.13,0.45}, {1.26,0.98,0.05}, {3.26,0.11,0.39}, {2.50,0.25,0.50}};
  const std::vector<double> gfac  = {0.2,1.0,1.8};
  const std::vector<double> kt0   = {0.5,1.0,1.5};
  const std::vector<double> zmins = {1.e-5,1.e-2,0.2};
  const std::vector<double> zmaxs = {0.2,0.9,0.9999};
  const std::vector<double> scale = {1.,400.,2500.};

  double worst = 0.;
  Frag_Norm fn;
  for (const Type& t : types) {
    for (const double zmin : zmins) {
      for (const double zmax : zmaxs) {
        if (zmax<=3.*zmin) continue;
        for (const double s : scale) {
          std::vector<double> cs;
          for (const double f : gfac) {
            const double g = t.gamma*f;
            cs.push_back(std::abs(g)>5.e-3 ? g*s : 0.);
          }
          for (const double k : kt0)
            cs.push_back(std::abs(t.gamma)>5.e-3 ? t.gamma*s/(k*k) : 0.);
          const double cnom = t.gamma*s;
          cs.push_back(cnom);
          const double cmin = *std::min_element(cs.begin(),cs.end());
          const double cmax = *std::max_element(cs.begin(),cs.end());
          fn.SetRange(zmin,zmax,t.alpha,t.beta,cmin,cmax);
          // one reference evaluation per exponent, then all pairs against the
          // nominal
          std::vector<double> err(cs.size());
          for (size_t i=0; i<cs.size(); ++i)
            err[i] = fn(t.alpha,t.beta,cs[i])
                   - reference(zmin,zmax,t.alpha,t.beta,cs[i]);
          const double enom = err.back();
          for (size_t i=0; i<cs.size(); ++i)
            worst = std::max(worst, std::abs(enom-err[i]));
        }
      }
    }
  }
  INFO("worst log-weight bias over the scan: " << worst);
  CHECK(worst < 1.e-4);
}
