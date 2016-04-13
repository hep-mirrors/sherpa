#include "SHRiMPS/Eikonals/Eikonal_Creator.H"
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"
#include "SHRiMPS/Tools/Kernels.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;
using namespace std;

Eikonal_Creator::Eikonal_Creator() :
  p_ff1(NULL), p_ff2(NULL),
  m_params(MBpars.GetEikonalParameters()),
  m_Bsteps(100),m_ff1steps(50), m_ff2steps(50) // was 400 & 100
{ }


void Eikonal_Creator::SetFormFactors(Form_Factor * ff1,Form_Factor * ff2) {
  p_ff1 = ff1; p_ff2 = ff2;
}

Omega_ik * Eikonal_Creator::InitialiseEikonal()
{
  Eikonal_Contributor * omegai = new Eikonal_Contributor(p_ff1,p_ff2,m_params);
  Eikonal_Contributor * omegak = new Eikonal_Contributor(p_ff1,p_ff2,m_params);
  FillBYGrids(omegai,omegak);

  Omega_ik * eikonal = new Omega_ik(m_params);
  eikonal->SetContributors(omegai,omegak);
  double prefactor(ATOOLS::sqr(p_ff1->Prefactor()*p_ff2->Prefactor()));
  eikonal->SetPrefactor(prefactor); 
  CreateImpactParameterGrid(eikonal);
  return eikonal;
}

void Eikonal_Creator::
FillBYGrids(Eikonal_Contributor * omegai,Eikonal_Contributor * omegak)
{
  msg_Info()<<METHOD<<" with lambda = "<<m_params.lambda<<", "
	    <<"Delta = "<<m_params.Delta<<") in "
	    <<"["<<(-m_params.Ymax)<<", "<<m_params.Ymax<<"].\n";

  omegai->PrepareGrid(m_ff1steps+1,m_ff2steps+1);
  omegak->PrepareGrid(m_ff1steps+1,m_ff2steps+1);
  double ff1max(p_ff1->Maximum());
  double ff2max(p_ff2->Maximum());
  double deltaff1(ff1max/double(m_ff1steps));
  double deltaff2(ff2max/double(m_ff2steps));
  double ff1, ff2;

  int ysteps(200);
  DEQ_Kernel_Base * deqkernel = 
    new DEQ_Kernel_NoKT(m_params.lambda,m_params.Delta,m_params.absorp);
  DEQ_Solver solver(deqkernel,2,deqmode::RungeKutta4);
  solver.SetInterval(-m_params.Ymax,m_params.Ymax);

  for (int i=0;i<m_ff1steps+1;i++) {
    for (int j=0;j<m_ff2steps+1;j++) {
      ff1 = Max(0.,ff1max-i*deltaff1);
      ff2 = Max(0.,ff2max-j*deltaff2);
      FixGridAndBorders(&solver,ysteps,ff1,ff2);
      omegai->InsertValues(i,j,solver.X()[0]);
      omegak->InsertValues(i,j,solver.X()[1]);
    }
  }
  delete deqkernel;
}


void Eikonal_Creator::FixGridAndBorders(DEQ_Solver * solver,int & ysteps,
					const double & ff1, const double & ff2)
{
  // these are the mock boundary conditions for the two terms
  // omega_i(k) and omega_k(i).  for y = -ymax, omega_{i(k)} = ff1, but
  // we have to guess the value of omega_{k(i)} there - we only know it
  // at y = ymax: there it is ff2.  we therefore have to iteratively
  // reconstruct a good starting point for oemga_{k(i)} for y = -ymax,
  // until we arrive a the correct boundary conditions.
  std::vector<double> x0(2,0.);
  x0[0] = ff1;
  x0[1] = ff2*(2.*m_params.Ymax)*
    exp(exp(-m_params.lambda/2.*(ff1+ff2))*m_params.Delta);

  int    n(0);
  double f_i(0.), x_i(0.), f_im1(f_i), x_im1(x_i), accu(0.005);
  std::vector<std::vector<double> > res; 
  do {
    solver->SolveSystem(x0,ysteps,accu);
    x_i = solver->X()[1][0];
    f_i = solver->X()[1][ysteps];
    if (n==0) x0[1] = ff2;
    else x0[1] = x_i-(f_i-ff2) * (x_i-x_im1)/(f_i-f_im1);
    x_im1 = x_i;
    f_im1 = f_i;
    n++;
  } while (dabs((f_i-ff2)/(f_i+ff2))>m_params.accu);
}

void Eikonal_Creator::CreateImpactParameterGrid(Omega_ik * eikonal)
{
  double B(0.), Bmax(m_params.bmax), deltaB(Bmax/double(m_Bsteps)), yref(0.);
  msg_Info()<<METHOD<<" up to B = "<<Bmax<<" in "<<m_Bsteps<<" steps, "
	    <<"delta = "<<deltaB<<".\n";
  eikonal->SetDeltaB(deltaB);

  double beta02(m_params.beta02), value(0.), accu(0.01);

  std::vector<double> * gridB(eikonal->GetImpactParameterGrid());
  std::vector<double> * gridBmax(eikonal->GetImpactParameterMaximumGrid());

  Integration_Kernel_B2 intkernel(eikonal->GetSingleTerm(0),
				  eikonal->GetSingleTerm(1));
  Gauss_Integrator integrator(&intkernel); 
  while (B<=Bmax) {
    intkernel.SetB(B);
    intkernel.SetYref(yref);
    intkernel.ResetMax();
    value = integrator.Integrate(0.,Bmax,accu,1)/beta02;
    if (dabs(value)<1.e-12) value  = 0.;
    gridB->push_back(value);
    gridBmax->push_back(intkernel.Max());
    B += deltaB;
  }
}
