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

//double Analytic_Eikonal::operator()(const double & B) const {
//  if (B<0.) return 0.;
//  double value(exp(-B*B*m_Lambda2/(4.*(2.+m_kappa_i+m_kappa_k))));
//  return m_prefactor*value;
//}

Eikonal_Creator::Eikonal_Creator(const Eikonal_Parameters & params) :
  p_ff1(NULL), p_ff2(NULL),
  m_params(params),
  m_Bsteps(400),m_ff1steps(100), m_ff2steps(100)
{  
  msg_Info()<<METHOD<<"("<<m_params.Ymax<<").\n";
}


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
  //CreateImpactParameterGrid(eikonal);
  return eikonal;
}

void Eikonal_Creator::
FillBYGrids(Eikonal_Contributor * omegai,Eikonal_Contributor * omegak)
{
  omegai->PrepareGrid(m_ff1steps+1,m_ff2steps+1);
  omegak->PrepareGrid(m_ff1steps+1,m_ff2steps+1);
  double ff1max(p_ff1->Maximum());
  double ff2max(p_ff2->Maximum());
  double deltaff1(ff1max/double(m_ff1steps));
  double deltaff2(ff2max/double(m_ff2steps));
  double ff1, ff2;

  msg_Out()<<METHOD<<" with lambda = "<<m_params.lambda<<", "
	   <<"Delta = "<<m_params.Delta<<") in "
	   <<"["<<(-m_params.Ymax)<<", "<<m_params.Ymax<<"].\n";
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
  msg_Tracking()<<METHOD<<"("<<ff1<<", "<<ff2<<"): new round.\n";
  std::vector<double> x0(2,0.);
  x0[0] = ff1;
  x0[1] = ff2*exp(exp(-m_lambda/2.*(ff1+ff2))*m_Delta)*(2.*m_Y);

  int    n(0);
  double f_i(0.), x_i(0.), f_im1(f_i), x_im1(x_i), accu(0.01);
  //accu(m_params.accu);
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
    msg_Tracking()<<"  n = "<<n<<" x_iml = "<<x_im1<<", f_iml = "<<f_im1
		  <<" --> "<<x0[1]<<" vs. "<<solver->X()[0][ysteps-1]<<".\n";
  } while (dabs((f_i-ff2)/(f_i+ff2))>m_params.accu);
  msg_Tracking()<<METHOD<<" reaches accuracy goal after "<<n<<" iterations, "
		<<"FormFactors = "<<ff1<<" & "<<ff2<<".\n";
}

void Eikonal_Creator::CreateImpactParameterGrid(Omega_ik * eikonal)
{
  m_b1max = m_b2max = m_params.bmax;
  std::vector<double> * gridB(eikonal->GetImpactParameterGrid());
  std::vector<double> * gridBmax(eikonal->GetImpactParameterMaximumGrid());

  Integration_Kernel_B2 intkernel(eikonal->GetSingleTerm(0),
				  eikonal->GetSingleTerm(1));
  Gauss_Integrator integrator(&intkernel); 

  double B(0.), deltaB(m_Bmax/double(m_Bsteps)), yref=0.;
  double value(0.);

  msg_Tracking()<<"   "<<METHOD<<" : "
		<<"Start producing impact parameter grid for "
		<<"{ik} = {"<<p_ff1->Number()<<" "<<p_ff2->Number()<<"}, "
		<<std::endl
		<<"   y = "<<yref<<", b_max = "<<m_Bmax<<" with "
		<<" b1max = "<<m_b1max<<"."<<std::endl;

  gridB->clear();
  gridBmax->clear();
  while (B<=1.0001*m_Bmax) {
    intkernel.SetB(B);
    intkernel.SetYref(yref);
    intkernel.ResetMax();
    value = integrator.Integrate(0.,m_b1max,m_accu,1)/m_beta02;
    if (dabs(value)<1.e-12) value  = 0.;
    gridB->push_back(value);
    gridBmax->push_back(intkernel.Max());
    msg_Tracking()<<"   B = "<<B
		  <<" Omega_{"<<p_ff1->Number()<<" "<<p_ff2->Number()<<"}"
		  <<" = "<<value<<" (max = "<<intkernel.Max()<<")."<<std::endl;
    B += deltaB;
  }
  //intkernel.PrintErrors();
  msg_Tracking()<<"   "<<METHOD<<" : Produced impact parameter grid.of size "
		<<gridB->size()<<std::endl
		<<"   and maximum grid of size "<<gridBmax->size()<<"."
		<<std::endl;
}

/*
void Eikonal_Creator::TestEikonal(Omega_ik * omegaik) const
{  
  Analytic_Eikonal eikonal(m_Delta,m_Y,p_ff1->Kappa(),p_ff2->Kappa(),
			   p_ff1->Beta0()*p_ff2->Beta0(),p_ff1->Lambda2());
  double B;
  for (int j=0;j<20;j++) {
    B = j*0.5;
    msg_Out()<<"  Omega_{ik}("<<B<<") : ana = "
	     <<eikonal(B)<<", num = "<<(*omegaik)(B)<<std::endl;
  }
  std::string filename("InclusiveQuantities/eikonals-ana.dat");
  std::ofstream was;
  was.open(filename.c_str());
  ysteps=100;
  was<<"# Delta = "<<m_Delta<<" Y = "<<m_Y<<" kappa_0 = "<<p_ff1->Kappa()<<" kappa_1 = "<<p_ff2->Kappa()
     <<" beta0_0 = "<<p_ff1->Beta0()<<" beta0_1 = "<<p_ff2->Beta0()<<" Lambda^2 = "<<p_ff1->Lambda2()<<std::endl;
  was<<"#  b1=b2    y     Omega_{1(2)}: num       ana     Omega_{(1)2}: num       ana "<<std::endl;
  for (int i=0;i<80;i++) {
    b1 = b2 = i*0.1;
    for (int j=0;j<ysteps;j++) {
      y        = -m_Y+j*(2.*m_Y)/double(ysteps);
      value12  = (*omegaik->GetSingleTerm(0))(b1,b2,y);
      value12a = ana12(b1,y);
      value21  = (*omegaik->GetSingleTerm(1))(b1,b2,y);
      value21a = ana21(b2,y);
      was<<b1<<"   "<<y<<"   "<<value12<<"   "<<value12a<<"  "
	 <<value21<<"   "<<value21a<<std::endl;
    }
    was<<std::endl<<std::endl;
  }
  was<<std::endl<<std::endl;
  was<<"# B    Omega_{ik}(B) : ana     num  "<<std::endl;
  for (int j=0;j<200;j++) {
    B = j*0.05;
    was<<B<<"   "<<eikonal(B)<<"   "<<(*omegaik)(B)<<std::endl;
  }
  was.close();
}
*/
