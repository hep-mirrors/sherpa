#include "SHRiMPS/Eikonals/Eikonal_Creator.H"
#include "SHRiMPS/Eikonals/Eikonal_Contributor.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "SHRiMPS/Tools/DEQ_Solver.H"
#include "SHRiMPS/Tools/Kernels.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

double Analytic_Contributor::
operator()(const double & b,const double & y) const {
  if (y<-m_Y || y>m_Y || b>p_ff->Bmax()) return 0.;
  return p_ff->FourierTransform(b)*exp(m_Delta*(m_Y+m_sign*y));
}

double Analytic_Eikonal::operator()(const double & B) const {
  if (B<0.) return 0.;
  double value(exp(-B*B*m_Lambda2/(4.*(2.+m_kappa_i+m_kappa_k))));
  return m_prefactor*value;
}

Eikonal_Creator::Eikonal_Creator(const int & test) :
  p_ff1(NULL), p_ff2(NULL),
  m_lambda((test==0)?MBpars("lambda"):0.), m_Delta(MBpars("Delta")), 
  m_beta02(sqr(MBpars("beta0"))), 
  m_absorp(MBpars.Absorption()),
  m_originalY(MBpars("originalY")), m_cutoffY(MBpars("deltaY")), 
  m_Y(m_originalY-m_cutoffY),
  m_Bmin(MBpars("bmin")), m_Bmax(MBpars("bmax")), m_Bsteps(400),
  m_test(test), m_accu(MBpars("accu")),
  m_ff1steps(100), m_ff2steps(100)
{  }

Omega_ik * Eikonal_Creator::CreateEikonal(Form_Factor * ff1,Form_Factor * ff2)
{
  p_ff1 = ff1;  p_ff2 = ff2; 

  msg_Tracking()
    <<METHOD<<"(lambda = "<<m_lambda<<", Delta = "<<m_Delta<<") "
    <<"in Y = "<<m_Y<<" "
    <<"(from "<<m_originalY<<" - "<<m_cutoffY<<")."<<std::endl
    <<"   Will now produce initial grids for FF = "
    <<p_ff1->FourierTransform(0.)<<" and "
    <<p_ff2->FourierTransform(0.)<<"."<<std::endl;

  Omega_ik * eikonal = new Omega_ik(p_ff1,p_ff2,m_Bsteps);
  CreateEikonalTerms(eikonal);
  CreateImpactParameterGrid(eikonal);
  if (m_test) TestEikonal(eikonal);

  return eikonal;
}

void Eikonal_Creator::CreateEikonalTerms(Omega_ik * eikonal)
{
  m_b1min   = m_b2min = 0.;
  m_b1max   = m_Bmax;
  m_b2max   = m_Bmax;
  Eikonal_Contributor * Omegai(eikonal->GetSingleTerm(0));
  Eikonal_Contributor * Omegak(eikonal->GetSingleTerm(1));

  m_ff1max  = p_ff1->FourierTransform(0.);
  m_ff2max  = p_ff2->FourierTransform(0.);
  double deltaff1(m_ff1max/double(m_ff1steps));
  double deltaff2(m_ff2max/double(m_ff2steps));
  double ff1, ff2;
  int ysteps(-1);
  Omegai->PrepareGrid(m_ff1steps+1,m_ff2steps+1);
  Omegak->PrepareGrid(m_ff1steps+1,m_ff2steps+1);

  DEQ_Kernel_Base * deqkernel(new DEQ_Kernel_NoKT(m_lambda,m_Delta,m_absorp));
  DEQ_Solver solver(deqkernel,2,deqmode::RungeKutta4);
  solver.SetInterval(-m_Y,m_Y);

  for (int i=0;i<m_ff1steps+1;i++) {
    ff1 = Max(0.,m_ff1max-i*deltaff1);
    for (int j=0;j<m_ff2steps+1;j++) {
      ff2 = Max(0.,m_ff2max-j*deltaff2);
      ysteps = FixBorders(&solver, ff1, ff2, ysteps);
      Omegai->InsertValues(i,j,solver.X()[0]);
      Omegak->InsertValues(i,j,solver.X()[1]);
    }
  }

  delete deqkernel;
}


int Eikonal_Creator::FixBorders(DEQ_Solver * solver,
				const double & ff1, const double & ff2, 
				const int & steps)
{
  std::vector<double> x0(2,0.);
  x0[0]  = ff1;
  x0[1]  = ff2*exp(exp(-m_lambda/2.*(ff1+ff2))*m_Delta)*(2.*m_Y);

  std::vector<std::vector<double> > res; 
  double accu(steps>0?-1.:m_accu), f_i(0.), x_i(0.), f_im1(f_i), x_im1(x_i);
  int    ysteps, n(0);
  do {
    ysteps = (accu>0.?64:steps);
    solver->SolveSystem(x0,ysteps,accu);
    res    = solver->X();
    x_i    = res[1][0];
    f_i    = res[1][ysteps];
    if (n==0) x0[1] = ff2;
    else {
      if (dabs((f_i-f_im1)/(f_i+f_im1))<1.e-12 ||
	  (dabs(f_i)<1.e-10 && dabs(f_im1)<1.e-10)) {
	break;
      }
      else x0[1] = x_i-(f_i-ff2) * (x_i-x_im1)/(f_i-f_im1);
    }
    x_im1 = x_i;
    f_im1 = f_i;
    n++;
    msg_Debugging()<<"   Done with the "<<n<<"th round, "
		   <<"size = "<<res[0].size()<<" : "
		   <<res[0][0]<<" --> "<<res[0][ysteps]<<", "
		   <<res[1][0]<<" --> "<<f_i<<"."<<std::endl;
  } while(dabs((f_i-ff2)/(f_i+ff2))>m_accu ||
	  f_i<ff2);

  return ysteps;
}

void Eikonal_Creator::CreateImpactParameterGrid(Omega_ik * eikonal)
{
  std::vector<double> * gridB(eikonal->GetImpactParameterGrid());
  std::vector<double> * gridBmax(eikonal->GetImpactParameterMaximumGrid());

  Integration_Kernel_B2 intkernel(eikonal->GetSingleTerm(0),
				  eikonal->GetSingleTerm(1));
  Gauss_Integrator integrator(&intkernel); 

  double B(m_Bmin), deltaB((m_Bmax-m_Bmin)/double(m_Bsteps)), yref=0.;
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


void Eikonal_Creator::TestEikonal(Omega_ik * omegaik) const
{
  if (m_test==1) {
    msg_Out()<<"In "<<METHOD<<":"<<std::endl
	     <<"   Check accuracy of DEQ solution vs. analytical result."
	     <<std::endl
	     <<"   To this end, set lambda = 0 ("<<m_lambda<<")."<<std::endl;
    Analytic_Contributor ana12(p_ff1,m_Delta,m_Y,+1);
    Analytic_Contributor ana21(p_ff2,m_Delta,m_Y,-1);

    double b1,b2,y,value12,value12a,value21,value21a;
    int ysteps(11);
    for (int i=0;i<8;i++) {
      b1 = b2 = i*1.;
      msg_Out()<<"  "<<" ff1("<<b1<<") = "<<p_ff1->FourierTransform(b1)<<","
	       <<" ff2("<<b2<<") = "<<p_ff2->FourierTransform(b2)<<std::endl; 
      for (int j=0;j<ysteps;j++) {
	y        = -m_Y+j*(2.*m_Y)/double(ysteps-1);
	value12  = (*omegaik->GetSingleTerm(0))(b1,b2,y);
	value12a = ana12(b1,y);
	value21  = (*omegaik->GetSingleTerm(1))(b1,b2,y);
	value21a = ana21(b2,y);
	msg_Out()<<"   y = "<<y<<"   "
		 <<"Omega_{1(2)}: num = "<<value12<<" (ana = "<<value12a<<"), "
		 <<"Omega_{(1)2}: num = "<<value21<<" (ana = "<<value21a<<")."
		 <<std::endl;
      }
    }
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
  //exit(1);
}
