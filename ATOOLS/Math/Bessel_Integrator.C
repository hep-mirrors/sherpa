#include "ATOOLS/Math/Bessel_Integrator.H"
#include "ATOOLS/Math/Special_Functions.H"
#include <iomanip>

using namespace ATOOLS;
using namespace std;

Bessel_Integrator::Bessel_Integrator(ATOOLS::Function_Base * f,const size_t & order) :
  m_kernel(f,order), m_order(order),
  m_maxbins(50), m_depth(10), m_iterator(1), m_xmin(0.) {
}

double Bessel_Integrator::operator()(const double & xmin, const double & xmax) {
  if (xmin!=0. || xmax!=-1) exit(1);
  m_xmin = xmin; m_xmax = xmax;
  if (FillBins(false)) return (m_M[m_maxdepth-1][0]/m_N[m_maxdepth-1][0]);
  return 0.;
}

bool Bessel_Integrator::FillBins(const bool & output) {
  if (output) msg_Out()<<"=== "<<METHOD<<":\n";
  FixBins(output);
  m_F.resize(m_maxbins+1,0.);
  m_Psi.resize(m_maxbins+1,0.);
  m_M.resize(m_depth);
  m_N.resize(m_depth);
  for (size_t i=0;i<m_depth;i++) {
    m_M[i].resize(m_maxbins+1,0.);
    m_N[i].resize(m_maxbins+1,0.);
  }
  m_maxdepth = m_depth;

  if (output) msg_Out()<<"=== "<<METHOD<<" start filling the supports, "
		       <<"max depth = "<<m_maxdepth<<":\n";
  Gauss_Integrator gauss(&m_kernel);
  double F = 0.;
  for (size_t i=1;i<m_maxbins+1;i++) {
    double xmin, xmax;
    if (m_iterator==0)      { xmin = m_zeroes[i-1];  xmax = m_zeroes[i]; }
    else if (m_iterator==1) { xmin = m_extrema[i-1]; xmax = m_extrema[i]; }
    else exit(1);
    F += m_Psi[i-1] = gauss.Integrate(xmin,xmax,1.e-6);
    if (dabs(m_Psi[i-1])<1.e-99 && i<m_maxdepth+1) {
      if (output) msg_Out()<<"   found a zero in integral over "
			   <<"x in ["<<xmin<<", "<<xmax<<"]: "<<m_Psi[i-1]<<" for "
			   <<"i = "<<i<<" --> new max depth = "<<(i-1)<<"\n";
      m_maxdepth = i-1;
      return false;
    }
    m_F[i-1]    = F;
    m_M[0][i-1] = m_F[i-1]/m_Psi[i-1];
    m_N[0][i-1] = 1./m_Psi[i-1];
    if (output) {
      msg_Out()<<"  bin i = "<<std::setw(2)<<i<<": "
	       <<"Psi["<<std::setw(8)<<xmin<<", "<<std::setw(8)<<xmax<<"] = "
	       <<std::setw(12)<<std::setprecision(6)<<m_Psi[i-1]<<", "
	       <<"F["<<std::setw(8)<<xmin<<", "<<std::setw(8)<<xmax<<"] = "
	       <<std::setw(12)<<std::setprecision(6)<<m_F[i-1]<<" --> "
	       <<"M_{-1}^("<<std::setw(2)<<i<<") = "
	       <<std::setw(12)<<std::setprecision(6)<<m_M[0][i-1]<<", and "
	       <<"N_{-1}^("<<std::setw(2)<<i<<") = "
	       <<std::setw(12)<<std::setprecision(6)<<m_N[0][i-1]<<"\n";
    }
  }
  for (size_t d=1;d<m_maxdepth;d++) {
    for (size_t i=1;i<m_maxbins+1-d;i++) {
      double xmin, xmax;
      if (m_iterator==0)      { xmin = m_zeroes[i];  xmax = m_zeroes[i+1+d]; }
      else if (m_iterator==1) { xmin = m_extrema[i]; xmax = m_extrema[i+1+d]; }
      else exit(1);
      double norm = 1./xmin - 1./xmax;
      m_M[d][i-1] = (m_M[d-1][i-1]-m_M[d-1][i])/norm;
      m_N[d][i-1] = (m_N[d-1][i-1]-m_N[d-1][i])/norm;
    }
    if (output && d!=m_maxdepth-1) {
      msg_Out()<<"=== "<<METHOD<<"(depth = "<<std::setw(2)<<d<<"): "
	       <<std::setw(12)<<std::setprecision(6)<<(m_M[d][0]/m_N[d][0])<<"\n";
    }
  }
  if (output)
    msg_Out()<<"=== "<<METHOD<<": Integral value("
	     <<"depth = "<<std::setw(2)<<m_maxdepth<<"): "
	     <<std::setw(12)<<std::setprecision(6)
	     <<(m_M[m_maxdepth-1][0]/m_N[m_maxdepth-1][0])<<" "
	     <<"for f(x) * J_{"<<m_order<<"}(x)\n";
  return true;
}

void Bessel_Integrator::FixBins(const bool & output) {
  m_zeroes.resize(m_maxbins+1);
  m_extrema.resize(m_maxbins+1);
  m_zeroes[0] = m_extrema[0] = m_xmin;
  for (size_t i=1;i<m_maxbins+1;i++) {
    switch (m_order) {
    case 0:
      if (i==1)  { m_zeroes[i] =  2.404825; break; }
      if (i==2)  { m_zeroes[i] =  5.520078; break; }
      if (i==3)  { m_zeroes[i] =  8.653727; break; }
      if (i==4)  { m_zeroes[i] = 11.791534; break; }
      if (i==5)  { m_zeroes[i] = 14.930917; break; }
      if (i==6)  { m_zeroes[i] = 18.071063; break; }
      if (i==7)  { m_zeroes[i] = 21.211636; break; }
      if (i==8)  { m_zeroes[i] = 24.352471; break; }
      if (i==9)  { m_zeroes[i] = 27.493479; break; }
      if (i==10) { m_zeroes[i] = 30.634606; break; }
    case 1:
      if (i==1)  { m_zeroes[i] =  3.831705; break; }
      if (i==2)  { m_zeroes[i] =  7.015586; break; }
      if (i==3)  { m_zeroes[i] = 10.173468; break; }
      if (i==4)  { m_zeroes[i] = 13.323691; break; }
      if (i==5)  { m_zeroes[i] = 16.470630; break; }
      if (i==6)  { m_zeroes[i] = 19.615858; break; }
      if (i==7)  { m_zeroes[i] = 22.760084; break; }
      if (i==8)  { m_zeroes[i] = 25.903672; break; }
      if (i==9)  { m_zeroes[i] = 29.046828; break; }
      if (i==10) { m_zeroes[i] = 32.189679; break; }
    case 2:
      if (i==1) { m_zeroes[i] =  5.1356; break; }
      if (i==2) { m_zeroes[i] =  8.4172; break; }
      if (i==3) { m_zeroes[i] = 11.6198; break; }
      if (i==4) { m_zeroes[i] = 14.7960; break; }
      if (i==5) { m_zeroes[i] = 17.9598; break; }
    default:
      if (i==1) {
	double n = double(i);
	m_zeroes[i] = ( n + 1.8557571*pow(n,1./3.) + 1.033150*pow(n,-1./3.) -
			0.00397/n - 0.0908*pow(n,-5./3.) + 0.043*pow(n,-7./3.));
      }
      else if (i==2) {
	double n = double(i);
	m_zeroes[i] = ( n + 3.2446076*pow(n,1./3.) + 3.158244*pow(n,-1./3.) -
			0.08331/n - 0.8437*pow(n,-5./3.) + 0.864*pow(n,-7./3.));
      }
      else m_zeroes[i] = 2.*m_zeroes[i-1]-m_zeroes[i-2];
      break;
    }
    m_extrema[i] = (m_zeroes[i]+m_zeroes[i-1])/2.;
  }
  if (output) {
    msg_Out()<<"=== "<<METHOD<<" yields the first "<<m_maxbins<<" "
	     <<"zeroes and extrema of Bessel J_{"<<m_order<<"}:\n";
    for (size_t i=0;i<m_maxbins+1;i++) {
      msg_Out()<<"  x_{"<<m_order<<", "<<std::setw(2)<<i<<"} = "
	       <<std::setw(12)<<std::setprecision(8)<<m_zeroes[i]<<" --> "
	       <<"y_{"<<m_order<<", "<<std::setw(2)<<i<<"} = "
	       <<std::setw(12)<<std::setprecision(8)<<m_extrema[i]<<".\n";
    }
  }
}

double Bessel_Integrator::Kernel::operator()(double x)  {
  return (*p_func)(x)*ATOOLS::SF.Jn(m_order,x);
}

double Bessel_Integrator::TestFunction::operator()(double x) {
  switch (m_mode) {
  case 4: return (1.-exp(-x))/(x*log(1.+sqrt(2.)));
  case 3: return x/(1.+x*x);
  case 2: return log(1.+x*x)/2.;
  case 1: return x*x;
  case 0:
  default:
    break;
  }
  return 1;
}

void Bessel_Integrator::Test() {
  Function_Base * f;
  for (size_t i=0;i<5;i++) {
    m_kernel.SetFunc(f = new Bessel_Integrator::TestFunction(i));
    m_kernel.SetOrder((i==2 || i==0) ? 1 : 0);
    FillBins(true);
    delete f;
  }
}

