#include "SHRiMPS/Eikonals/Omega_ik.H"
#include "SHRiMPS/Tools/MinBias_Parameters.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Message.H"

using namespace SHRIMPS;
using namespace ATOOLS;

Omega_ik::Omega_ik(const Eikonal_Parameters & params) :
  m_weights(params), m_bmax(params.bmax),
  p_Omegaik(0), p_Omegaki(0) 
{ 
  m_gridB.clear();
  m_gridBmax.clear();
  m_gridD.clear();
}

Omega_ik::~Omega_ik() {
  if (p_Omegaik) delete p_Omegaik;
  if (p_Omegaki) delete p_Omegaki;
}

void Omega_ik::SetContributors(Eikonal_Contributor * Omegaik,
			       Eikonal_Contributor * Omegaki) { 
  p_Omegaik = Omegaik; p_Omegaki = Omegaki; 
  m_weights.SetSingleOmegaTerms(p_Omegaik,p_Omegaki);
}

void Omega_ik::SetDeltaB(const double & deltaB) {
  m_deltaB = deltaB;
}

void Omega_ik::SetPrefactor(const double & prefactor) { 
  m_prefactor = prefactor; 
}

Eikonal_Contributor * Omega_ik::GetSingleTerm(const int & i) {
  if (i==0) return p_Omegaik;
  if (i==1) return p_Omegaki;
  msg_Error()<<"Error in "<<METHOD<<"("<<i<<"):"<<std::endl
	     <<"   Out of range.  Will exit the run."<<std::endl;
  exit(1);  
}

Eikonal_Weights * Omega_ik::GetWeights() {
  return &m_weights;
}

std::vector<double> * Omega_ik::GetImpactParameterGrid() { 
  return &m_gridB; 
}

std::vector<double> * Omega_ik::GetImpactParameterMaximumGrid() { 
  return &m_gridBmax; 
}

double Omega_ik::operator()(const double & B) const {
  if (B<0. || B>=m_bmax) return 0.;
  size_t Bbin(int(B/m_deltaB));
  return ((m_gridB[Bbin]*((Bbin+1)*m_deltaB-B)+
	   m_gridB[Bbin+1]*(B-Bbin*m_deltaB))/m_deltaB);
}

ATOOLS::Vec4D Omega_ik::SelectB1B2(double & b1,double & b2,const double & B) {
  double maxvalue(1.1*Maximum(B));
  double theta(0.),b1max(p_Omegaik->B1max()),value(0.);  
  bool   accept(false);
  while (!accept) {
    theta = 2.*M_PI*ATOOLS::ran->Get();
    b1    = b1max * ATOOLS::ran->Get();
    b2    = sqrt(B*B+b1*b1-2.*B*b1*cos(theta));
    if (b1>m_bmax || b2>m_bmax) continue;
    value = b1*(*p_Omegaik)(b1,b2,0.)*(*p_Omegaki)(b1,b2,0.);
    if (value>maxvalue) 
      msg_Error()<<"Warning in "<<METHOD<<"("<<b1<<", "<<b2<<", "<<B<<"):"
		 <<std::endl
		 <<"   Value = "<<value
		 <<" > maxvalue = "<<maxvalue<<"."<<std::endl;
    if (value/maxvalue>ATOOLS::ran->Get()) accept=true;
  }
  return Vec4D(0.,b1*cos(theta),b1*sin(theta),0.);
}

double Omega_ik::Maximum(const double & B) const {
  if (B<0. || B>=m_bmax) return 0.;
  size_t Bbin(int(B/m_deltaB));
  return ((m_gridBmax[Bbin]*((Bbin+1)*m_deltaB-B)+
	   m_gridBmax[Bbin+1]*(B-Bbin*m_deltaB))/m_deltaB);
}

void Omega_ik::PrepareQT(const double & b1,const double & b2) {  
  double D1,D2,invD,y;
  p_Omegaik->SetB1B2(b1,b2);
  p_Omegaki->SetB1B2(b1,b2);
  Gauss_Integrator inti(p_Omegaik), intk(p_Omegaki);
  m_gridD.clear();
  for (int i=0;i<=m_Ysteps;i++) {
    y    = m_Y*(1.-2.*double(i)/m_Ysteps);
    D1   = inti.Integrate(-m_Y,y,2.e-2,1);
    D1  += intk.Integrate(-m_Y,y,2.e-2,1);
    D2   = inti.Integrate(y,m_Y,2.e-2,1);
    D2  += intk.Integrate(y,m_Y,2.e-2,1);
    invD = 1./D1+1./D2;
    m_gridD.push_back(invD);
  }
}

double Omega_ik::
EffectiveIntercept(double b1,double b2,const double & y) {
  return m_weights.EffectiveIntercept(b1,b2,y);
}


void Omega_ik::TestEikonal(Analytic_Eikonal * anaeik,
			   const std::string & dirname) const
{
  std::string filename(dirname+"eikonals-ana.dat");
  std::ofstream was;
  was.open(filename.c_str());
  was<<"# B    Omega_{ik}(B) : ana     num  "<<std::endl;
  double errmax(0.), B, ana, num;
  for (int j=0;j<200;j++) {
    B   = j*0.05;
    ana = (*anaeik)(B);
    num = (*this)(B);
    if (2.*(ana-num)/(ana+num)>errmax) errmax = 2.*(ana-num)/(ana+num);
    was<<B<<"   "<<ana<<"   "<<num<<"\n";
  }
  was.close();
  msg_Info()<<METHOD<<" with maximal error: "<<(100.*errmax)<<" %.\n";
}

void Omega_ik::TestIndividualGrids(Analytic_Contributor * ana12,
				   Analytic_Contributor * ana21,
				   const double & Ymax,
				   const std::string & dirname) const
{
  std::ofstream was;
  std::string filename = dirname+std::string("/SingleTerms_b1_0_b2_0.dat");
  was.open(filename.c_str());
  msg_Out()<<"In "<<METHOD<<":"<<std::endl
	   <<"   Check accuracy of DEQ solution vs. analytical result "
	   <<"in Y-range ["<<-Ymax<<", "<<Ymax<<"].\n"
	   <<"   To this end, we have set lambda = 0.\n";
  double b1,b2,y,value12,value12a,value21,value21a,maxerr(0.);
  int ysteps(20);
  for (int i=0;i<10;i++) {
    for (int j=i;j<11;j++) {
      b1 = i*0.5;
      b2 = j*0.5;
      for (int k=0;k<ysteps;k++) {
	y        = -Ymax+k/double(ysteps-1)*(2.*Ymax);
	value12  = (*p_Omegaik)(b1,b2,y);
	value12a = (*ana12)(b1,y);
	value21  = (*p_Omegaki)(b1,b2,y);
	value21a = (*ana21)(b2,y);
	if (i==0 && j==0) {
	  was<<y
	     <<" "<<value12<<" "<<value21<<" "
	     <<value12a<<" "<<value21a<<".\n";
	}
	if (value12>1.e-6 && value12a>1.e-6 &&
	    (value12-value12a)/(value12+value12a)>maxerr)
	  maxerr = (value12-value12a)/(value12+value12a);
	if (value21>1.e-6 && value21a>1.e-6 &&
	    (value21-value21a)/(value21+value21a)>maxerr)
	  maxerr = (value21-value21a)/(value21+value21a);
	msg_Out()<<"   y = "<<y<<"   "
		 <<"Omega_{1(2)} = "<<value12<<" (ana = "<<value12a<<"), "
		 <<"Omega_{(1)2} = "<<value21<<" (ana = "<<value21a<<")\n.";
      }
    }   
  }
  msg_Out()<<"Maximal error: "<<(100.*maxerr)<<" %.\n";
  was.close();
}
