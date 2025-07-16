#include "SHRiMPS/Tools/DEQ_Solver.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/MathTools.H"

using namespace SHRIMPS;
using namespace ATOOLS;

DEQ_Solver::DEQ_Solver(DEQ_Kernel_Base * kernel,const size_t & dim,
		       const deqmode::code & deq,const int & test) :
  p_kernel(kernel), m_dim(dim), 
  m_x(std::vector<std::vector<double> >(m_dim)),
  m_xsave(std::vector<std::vector<double> >(m_dim)),
  m_deqmode(deq), m_test(test)
{ }

void DEQ_Solver::InitIteration(const std::vector<double> & x0,
			       const int & steps) {
  m_stepsize = (m_yend-m_ystart)/double(steps);
  for (size_t i=0;i<m_dim;i++) {
    m_x[i].clear();
    m_x[i].resize(steps+1);
    m_x[i][0] = x0[i];
  }
  msg_Tracking()<<METHOD<<"(steps = "<<steps<<": "
		<<"x0[0] = "<<m_x[0][0]<<", x0[1] = "<<m_x[1][0]<<".\n";
}

void DEQ_Solver::RunIteration(const int & steps) {
  switch (m_deqmode) {
  case deqmode::SimpleEuler:
    SimpleEuler(steps);
    break;
  case deqmode::RungeKutta2:
    RungeKutta2(steps);
    break;
  case deqmode::RungeKutta4:
  default:
    RungeKutta4(steps);
    break;
  }
}

bool DEQ_Solver::CheckAccuracy(const double & accu,double & diffmax) {
  // m_x is the grid with twice as many steps as m_xsave
  msg_Tracking()<<METHOD<<"("<<m_x[0].size()<<"): ";
  double diff(0.);
  diffmax = -1.;
  for (size_t j=0;j<m_xsave[0].size();j++) {
    for (size_t i=0;i<m_dim;i++) {
      diff = dabs(2.*(m_xsave[i][j]-m_x[i][2*j])/(m_xsave[i][j]+m_x[i][2*j]));
      if (diff>accu) {
	msg_Tracking()<<" --> diff = "<<diff<<" from ["<<i<<", "<<j<<"] --> "
		      <<m_xsave[i][j]<<" and "<<m_x[i][2*j]<<".\n";
	return false;
      }
      if (diff>diffmax) diffmax = diff;
    }
  }
  msg_Tracking()<<" --> ok.\n";
  return true;
}

void DEQ_Solver::SaveResult() {
  msg_Tracking()<<METHOD<<" for size = "<<m_x.size()<<"("<<m_x[0].size()<<").\n";
  for (size_t i=0;i<m_dim;i++) {
    m_xsave[i].clear();
    m_xsave[i].resize(m_x[i].size());
    //copy(m_x[i].begin(),m_x[i].end(),back_inserter(m_xsave[i]));
    copy(m_x[i].begin(),m_x[i].end(),m_xsave[i].begin());
  }
}

void DEQ_Solver::RestoreResult() {
  for (size_t i=0;i<m_dim;i++) {
    m_x[i].clear();
    copy(m_xsave[i].begin(),m_xsave[i].end(),back_inserter(m_x[i]));
  }
}

void DEQ_Solver::IncreaseAccuracy(int & steps) {
  steps      *= 2;
  m_stepsize /= 2.;
}

void DEQ_Solver::SolveSystem(const std::vector<double> & x0,int & steps,
			     const double & accu)
{
  msg_Tracking()<<"In "<<METHOD<<"(steps = "<<steps<<", accu = "<<accu<<").\n";
  if (x0.size()!=m_dim) exit(1);
  double diffmax(-1.);
  int iterations(0);
  bool run(true);
  // Initiate and store reference with half the steps.
  InitIteration(x0,int(steps/2));
  RunIteration(int(steps/2));
  SaveResult();
  do {
    InitIteration(x0,steps);
    RunIteration(steps);
    if (CheckAccuracy(accu,diffmax)) {
      run = false;
    }
    else {
      // update reference and increase number of steps by two.
      SaveResult();
      IncreaseAccuracy(steps);
    }
  } while (run && iterations++<5);
  msg_Tracking()<<"Out "<<METHOD<<"(steps = "<<steps<<", accu = "<<accu<<") "
		<<"yields final accuracy "<<diffmax<<".\n";
}


void DEQ_Solver::SimpleEuler(const int & steps) {
  std::vector<double> x1, f1;
  x1.resize(m_dim);
  f1.resize(m_dim);
  for (size_t i(0);i<m_dim;i++) x1[i] = m_x[i][0];

  for (int step(0);step<steps;step++) {
    f1  = (*p_kernel)(x1);
    for (size_t i(0);i<m_dim;i++) x1[i] += m_stepsize*f1[i];
    for (size_t i(0);i<m_dim;i++) m_x[i][step+1] = x1[i];
  }
}

void DEQ_Solver::RungeKutta2(const int & steps) {
  std::vector<double> x1, x2, f1, f2;
  x1.resize(m_dim);
  f1.resize(m_dim);
  x2.resize(m_dim);
  f2.resize(m_dim);
  for (size_t i(0);i<m_dim;i++) x1[i] = m_x[i][0];

  for (int step(0);step<steps;step++) {
    f1  = (*p_kernel)(x1);
    for (size_t i(0);i<m_dim;i++) x2[i]  = x1[i] + m_stepsize/2.*f1[i];
    f2  = (*p_kernel)(x2);
    for (size_t i(0);i<m_dim;i++) x1[i] += m_stepsize*f2[i];
    for (size_t i(0);i<m_dim;i++) m_x[i][step+1] = x1[i];
  }
}

void DEQ_Solver::RungeKutta4(const int & steps) {
  std::vector<double> x1, x2, x3, x4, f1, f2, f3, f4;
  x1.resize(m_dim);
  f1.resize(m_dim);
  x2.resize(m_dim);
  f2.resize(m_dim);
  x3.resize(m_dim);
  f3.resize(m_dim);
  x4.resize(m_dim);
  f4.resize(m_dim);
  for (size_t i(0);i<m_dim;i++) x1[i] = m_x[i][0];

  for (int step(0);step<steps;step++) {
    f1  = (*p_kernel)(x1);
    for (size_t i(0);i<m_dim;i++) x2[i]  = x1[i] + m_stepsize/2.*f1[i];
    f2  = (*p_kernel)(x2);
    for (size_t i(0);i<m_dim;i++) x3[i]  = x1[i] + m_stepsize/2.*f2[i];
    f3  = (*p_kernel)(x3);
    for (size_t i(0);i<m_dim;i++) x4[i]  = x1[i] + m_stepsize*f3[i];
    f4  = (*p_kernel)(x4);


    for (size_t i(0);i<m_dim;i++) 
      x1[i] += m_stepsize*(f1[i]+2.*f2[i]+2.*f3[i]+f4[i])/6.;
    for (size_t i(0);i<m_dim;i++) 
      m_x[i][step+1] = x1[i];
  }
}
