#include "SINIC++/Eikonals/Single_Channel_Eikonal.H"
#include "ATOOLS/Math/Gauss_Integrator.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/MyStrStream.H"

using namespace SINIC;
using namespace ATOOLS;

Single_Channel_Eikonal::
Single_Channel_Eikonal(const deqmode::code & deq,const int & test) :
  m_test(test), m_deqmode(deq), p_ff1(NULL), p_ff2(NULL), 
  p_convolution2D(NULL), p_integrator(NULL)
{ 
}

Single_Channel_Eikonal::~Single_Channel_Eikonal() {
  if (p_ff1) { delete p_ff1; p_ff1 = NULL; }
  if (p_ff2) { delete p_ff2; p_ff2 = NULL; }
  if (p_convolution2D) { delete p_convolution2D; p_convolution2D = NULL; }
  if (p_integrator)    { delete p_integrator;    p_integrator    = NULL; }
}

void Single_Channel_Eikonal::
Initialise(Form_Factor * ff1,Form_Factor * ff2,
	   const double & lambda,const double & alpha,
	   const double & Y,const double & ycutoff)
{
  p_ff1       = ff1; 
  m_ff1max    = p_ff1->FFmax();
  m_ff1bins   = p_ff1->Bbins();
  m_deltaff1  = m_ff1max/m_ff1bins; 
  m_b1max     = p_ff1->Bmax();
  p_ff2       = ff2; 
  m_ff2max    = p_ff2->FFmax();
  m_ff2bins   = p_ff2->Bbins();
  m_deltaff2  = m_ff2max/m_ff2bins; 
  m_b2max     = p_ff2->Bmax();
  m_Bmax      = Max(m_b1max,m_b2max);
  m_Bbins     = Max(m_ff1bins,m_ff2bins);
  m_deltaB    = m_Bmax/m_Bbins;

  m_beta2     = p_ff1->Beta0()*p_ff2->Beta0();

  m_lambda    = lambda;
  m_alpha     = alpha;
  m_expfactor = 1./2.;

  m_Y         = Y;
  m_ycutoff   = ycutoff;
  m_yshift    = m_Y-m_ycutoff;
  m_ybins     = 20;
  m_deltay    = 2.*m_yshift/m_ybins;

  if (m_test==10) {
    m_ycutoff = 0.;
    m_yshift  = m_Y;
    m_lambda  = 0.;
  }

  m_accu      = 1.e-2;
  msg_Tracking()<<METHOD<<"(lambda = "<<m_lambda<<", alpha = "<<m_alpha<<") "
		<<"in Y = "<<m_yshift<<":"<<std::endl
		<<"   Form factors: ff1,2max = "<<m_ff1max<<", "<<m_ff2max
		<<" in b up to "<<m_b1max<<", "<<m_b2max<<"."<<std::endl
		<<"   Will now produce initial grids."<<std::endl;
  if (m_test==0 || m_test==1 || m_test==10) ProduceInitialGrids();
  if (m_test==10) TestSingleEikonal();
}

void Single_Channel_Eikonal::ProduceInitialGrids()
{
  msg_Tracking()<<"In "<<METHOD<<"("<<m_ff1bins<<" * "<<m_ff2bins<<")."
		<<std::endl;
  
  p_convolution2D = new Convolution2D(this,m_yshift,m_test);
  p_integrator    = new ATOOLS::Gauss_Integrator(p_convolution2D);
  m_grid1 = 
    std::vector<std::vector<std::vector<double> > >
    (m_ff1bins+1,std::vector<std::vector<double> >
     (m_ff2bins+1,std::vector<double>(2,0.)));
  m_grid2 =
    std::vector<std::vector<std::vector<double> > >
    (m_ff1bins+1,std::vector<std::vector<double> >
     (m_ff2bins+1,std::vector<double>(2,0.)));

  double value1, value2;
  for (int i=0;i<=m_ff1bins;i++) {
    for (int j=0;j<=m_ff2bins;j++) {
      InitialiseBoundaries(i,j,value1,value2);
      if (i==0 && j==0) m_ybins = AdjustGrid(i,j,value1,value2);
                   else SolveSystem(i,j,value1,value2,m_ybins);
    }
  }
  msg_Tracking()<<METHOD<<": "
		<<"Produced initial grids of DEq solutions."<<std::endl;
}

void Single_Channel_Eikonal::InitialiseBoundaries(const int & i,const int & j,
						  double & val1,double & val2) 
{
  val1 = m_ff1max-i*m_deltaff1;
  val2 = m_ff2max-j*m_deltaff2;
}

int Single_Channel_Eikonal::
AdjustGrid(const int & i,const int & j,double & val1,double & val2) {
  int ysteps(m_ybins);
  SolveSystem(i,j,val1,val2,ysteps);
  std::vector<double> comp1, comp2;
  do {
    comp1 = m_grid1[i][j];
    comp2 = m_grid2[i][j];

    ysteps   *= 2;
    m_deltay /= 2.;
    SolveSystem(i,j,val1,val2,ysteps);
  } while (CheckAccuracy(i,j,ysteps,comp1,comp2));
  return ysteps;
}


void Single_Channel_Eikonal::
SolveSystem(const int & i,const int & j,double & val1,
	    double & val2,const int & steps)
{
  switch (m_deqmode) {
  case deqmode::RungeKutta2:
    return RungeKutta2(i,j,val1,val2,steps);
  case deqmode::RungeKutta4:
    return RungeKutta4(i,j,val1,val2,steps);
  case deqmode::RungeKutta4Transformed:
  default:
    return RungeKutta4Transformed(i,j,val1,val2,steps);
  }
}


bool Single_Channel_Eikonal::
CheckAccuracy(const int & i,const int & j,const int & ysteps,
	      const std::vector<double> & comp1,
	      const std::vector<double> & comp2)
{
  double diff1(0.), diff2(0.), ave1(0.), ave2(0.);
  diff1 = diff2 = 0.;
  for (int step=2;step<ysteps;step+=2) {
    double test((comp1[step/2]+comp1[step/2+1])/2.),  x1(m_grid1[i][j][step]);
    double testp((comp2[step/2]+comp2[step/2+1])/2.), xp1(m_grid2[i][j][step]);
    
    ave1 += dabs(test/x1-1.);
    ave2 += dabs(testp/xp1-1.);
    if (dabs(test/x1-1.)  > diff1) {
      diff1 = dabs(test/x1-1.);
    } 
    if (dabs(testp/xp1-1.)> diff2) {
      diff2 = dabs(testp/xp1-1.); 
    }
  }
  if (diff1>m_accu || diff2>m_accu) {
    return true;
  }
  return false;
}  

void Single_Channel_Eikonal::
RungeKutta4Transformed(const int & i,const int & j,
		       double & val1,double & val2,const int & steps) {
  double norm_ik(val1),norm_ki(val2),Omega_ik(norm_ik),Omega_ki(norm_ki);
  double deltay(m_deltay),y(0.);
  double barDelta(m_alpha*exp(-m_lambda*m_expfactor*(norm_ik+norm_ki)));
  double x1_ik(1.),x1_ki(1.),exp1,x2_ik,x2_ki,exp2;
  double x3_ik,x3_ki,exp3,x4_ik,x4_ki,exp4;
  double f1_ik,f1_ki,f2_ik,f2_ki,f3_ik,f3_ki,f4_ik,f4_ki;
  m_grid1[i][j].clear();
  m_grid2[i][j].clear();
  m_grid1[i][j].push_back(Omega_ik);
  m_grid2[i][j].push_back(Omega_ki);

  for (int step=0;step<steps;step++) {
    // x_1 = x(y_1), y_1=y, f(x_1,y_1) = x_1*barDelta*exp1
    exp1   = exp(-m_lambda*m_expfactor*norm_ik*(exp(barDelta*y)*x1_ik-1.)
		 -m_lambda*m_expfactor*norm_ki*(exp(barDelta*y)*x1_ki-1.))-1.;
    f1_ik  = x1_ik*barDelta*exp1;
    f1_ki  = x1_ki*barDelta*exp1; 
    // x_2 = x_1 + delta y/2 * f(x_1,y_1), y_2 = y+delta y/2, 
    // f(x_2,y_2) = x_2*barDelta*exp2
    x2_ik  = x1_ik + deltay/2. * f1_ik;
    x2_ki  = x1_ki + deltay/2. * f1_ki;
    exp2   = exp(-m_lambda*m_expfactor*norm_ik*(exp(barDelta*(y+deltay/2.))*x2_ik-1.)
		 -m_lambda*m_expfactor*norm_ki*(exp(barDelta*(y+deltay/2.))*x2_ki-1.))-1.;
    f2_ik  = x2_ik*barDelta*exp2;
    f2_ki  = x2_ki*barDelta*exp2; 
    // x_3 = x_1 +delta y/2 * f(x_2,y_2), y_3 = y+delta y/2, 
    // f(x_3,y_3) = x_3*barDelta*exp3
    x3_ik  = x1_ik + deltay/2. * f2_ik;
    x3_ki  = x1_ki + deltay/2. * f2_ki;
    exp3   = exp(-m_lambda*m_expfactor*norm_ik*(exp(barDelta*(y+deltay/2.))*x3_ik-1.)
		 -m_lambda*m_expfactor*norm_ki*(exp(barDelta*(y+deltay/2.))*x3_ki-1.))-1.;
    f3_ik  = x3_ik*barDelta*exp3;
    f3_ki  = x3_ki*barDelta*exp3;
    // x_4 = x_1 +delta y * f(x_3,y_3), y_4 = y+delta y/, 
    // f(x_4,y_4) = x_4*barDelta*exp4 
    x4_ik  = x1_ik + deltay * f3_ik;
    x4_ki  = x1_ki + deltay * f3_ki;
    exp4   = exp(-m_lambda*m_expfactor*norm_ik*(exp(barDelta*(y+deltay))*x4_ik-1.)
		 -m_lambda*m_expfactor*norm_ki*(exp(barDelta*(y+deltay))*x4_ki-1.))-1.;
    f4_ik  = x4_ik*barDelta*exp4;
    f4_ki  = x4_ki*barDelta*exp4;


    x1_ik += deltay*(f1_ik+2.*f2_ik+2.*f3_ik+f4_ik)/6.;
    x1_ki += deltay*(f1_ki+2.*f2_ki+2.*f3_ki+f4_ki)/6.;
    
    y     += deltay;

    Omega_ik = norm_ik*exp(barDelta*y)*x1_ik;
    Omega_ki = norm_ki*exp(barDelta*y)*x1_ki;
    m_grid1[i][j].push_back(Omega_ik);
    m_grid2[i][j].push_back(Omega_ki);
  }
}

void Single_Channel_Eikonal::RungeKutta4(const int & i,const int & j,
					 double & val1,double & val2,
					 const int & steps) {
  double x1_ik(val1),x1_ki(val2),exp1,x2_ik,x2_ki,exp2;
  double x3_ik,x3_ki,exp3,x4_ik,x4_ki,exp4;
  double f1_ik,f1_ki,f2_ik,f2_ki,f3_ik,f3_ki,f4_ik,f4_ki,deltay(m_deltay);

  m_grid1[i][j].clear();
  m_grid2[i][j].clear();
  m_grid1[i][j].push_back(x1_ik);
  m_grid2[i][j].push_back(x1_ki);

  msg_Tracking()<<" y = "<<(-m_yshift)<<": "
		<<"Omega_ik = "<<x1_ik<<", Omega_ki = "<<x1_ki<<" "
		<<"(expterm = "<<exp(m_alpha*m_yshift)<<")."<<std::endl;

  for (int step=0;step<steps;step++) {
    exp1   = exp(-m_lambda*m_expfactor*(x1_ik+x1_ki));
    f1_ik  = x1_ik*m_alpha*exp1; 
    f1_ki  = x1_ki*m_alpha*exp1; 
    x2_ik  = x1_ik+deltay/2.*f1_ik;
    x2_ki  = x1_ki+deltay/2.*f1_ki;
    exp2   = exp(-m_lambda*m_expfactor*(x2_ik+x2_ki));
    f2_ik  = x2_ik*m_alpha*exp2; 
    f2_ki  = x2_ki*m_alpha*exp2; 
    x3_ik  = x1_ik+deltay/2.*f2_ik;
    x3_ki  = x1_ki+deltay/2.*f2_ki;
    exp3   = exp(-m_lambda*m_expfactor*(x3_ik+x3_ki));
    f3_ik  = x3_ik*m_alpha*exp3; 
    f3_ki  = x3_ki*m_alpha*exp3; 
    x4_ik  = x1_ik+deltay*f3_ik;
    x4_ki  = x1_ki+deltay*f3_ki;
    exp4   = exp(-m_lambda*m_expfactor*(x4_ik+x4_ki));
    f4_ik  = x4_ik*m_alpha*exp4; 
    f4_ki  = x4_ki*m_alpha*exp4; 

    x1_ik += deltay*(f1_ik+2.*f2_ik+2.*f3_ik+f4_ik)/6.;
    x1_ki += deltay*(f1_ki+2.*f2_ki+2.*f3_ki+f4_ki)/6.;

    m_grid1[i][j].push_back(x1_ik);
    m_grid2[i][j].push_back(x1_ki);
  }
}



void Single_Channel_Eikonal::RungeKutta2(const int & i,const int & j,
					 double & val1,double & val2,const int & steps) {
  double x1_ik(val1),x1_ki(val2),x2_ik,x2_ki,exp1,exp2,f1_ik,f1_ki,f2_ik,f2_ki;

  m_grid1[i][j].clear();
  m_grid2[i][j].clear();
  m_grid1[i][j].push_back(x1_ik);
  m_grid2[i][j].push_back(x1_ki);

  for (int step=0;step<steps;step++) {
    exp1   = exp(-m_lambda*m_expfactor*(x1_ik+x1_ki));
    f1_ik  = x1_ik * m_alpha * exp1;
    f1_ki  = x1_ki * m_alpha * exp1;
    x2_ik  = x1_ik + m_deltay/2.*f1_ik;
    x2_ki  = x1_ki + m_deltay/2.*f1_ki;
    exp2   = exp(-m_lambda*m_expfactor*(x2_ik+x2_ki));
    f2_ik  = x2_ik * m_alpha * exp2;
    f2_ki  = x2_ki * m_alpha * exp2;
    x1_ik += m_deltay*f2_ik;
    x1_ki += m_deltay*f2_ki;

    m_grid1[i][j].push_back(x1_ik);
    m_grid2[i][j].push_back(x1_ki);
  }
}






void Single_Channel_Eikonal::ProduceImpactParameterGrid(const double & y) {
  bool   run(true);
  double B(0.), value, valmax, valmin;
  msg_Tracking()<<METHOD<<" : Start producing impact parameter grid for "
		<<"y = "<<y<<", b_max = "<<m_Bmax<<"."<<std::endl;
  while (run) {
    m_gridB.clear();
    valmax = -1.;
    valmin = -1.;
    while (B<=m_Bmax) {
      value = IntegrateOutImpactParameters(B,y);
      if (dabs(value)<1.e-12) value  = 0.;
      if (value>0.) {
	if (value<valmin) valmin = value;
	if (value>valmax) valmax = value;
      }
      m_gridB.push_back(value);
      B += m_deltaB;
    }
    run = false;
    if (valmin/valmax>m_accu) {
      m_Bbins *= 2;
      m_Bmax  *= 2.;
      msg_Tracking()<<METHOD<<" does not meet accuracy goal in B_max = "
		    <<m_Bmax<<std::endl;
      run = true;
    }
    if (!run) {
      for (size_t i=0;i<m_gridB.size();i++) {
	double Btest((i+0.5)*m_deltaB);
	double fit((*this)(Btest));
	double exact(IntegrateOutImpactParameters(Btest,y));
	if (dabs(fit/exact-1.)>m_accu && sqr(exact/valmax)>m_accu) {
	  msg_Tracking()<<METHOD<<" does not meet accuracy goal "
			<<"("<<(dabs(fit/exact-1.)>0.01)<<") "
			<<"in "<<m_Bbins<<" steps:"<<std::endl
			<<" i = "<<i<<", B = "<<Btest<<": "
			<<dabs(fit/exact-1.)<<" from :"
			<<exact<<" vs. "<<fit<<"."<<std::endl
			<<"   Use now "<<m_Bbins<<" steps --> "
			<<"delta_B = "<<m_deltaB<<"."<<std::endl;
	  B         = 0.;
	  m_Bbins  *= 2;
	  m_deltaB /= 2.;
	  run = true;
	  break;
	}
      }
    }
  }
  msg_Tracking()<<METHOD<<" : Produced impact parameter grid."<<std::endl;
  p_convolution2D->PrintErrors();
} 


double Single_Channel_Eikonal::
IntegrateOutImpactParameters(const double & B,const double & Y) {
  p_convolution2D->SetB(B);
  p_convolution2D->SetY(Y);
  return p_integrator->Integrate(0.,m_b1max,m_accu,1)/m_beta2;
}

bool Single_Channel_Eikonal::
GeneratePositions(const double & B,double & b1,double & theta1) {
  double b2, value, y(0.);
  do {
    b1     = ran->Get()*m_b1max;
    theta1 = 2.*M_PI*ran->Get();
    b2     = sqrt(B*B+b1*b1-2.*B*b1*cos(theta1));
    value  = Omega12(b1,b2,y)*Omega21(b1,b2,y);
  } while (value<m_maxconv*ran->Get());
  return true;
}

double Single_Channel_Eikonal::
Omega12(const double & b1,const double & b2,const double & y,
	const bool & plot) const 
{
  if (b1>m_b1max || b1<0. || 
      b2>m_b2max || b2<0. || 
      y>m_yshift || y<-m_yshift) return 0.;
  double ff1(p_ff1->FourierTransform(b1)), ff2(p_ff2->FourierTransform(b2));
  int    ff1bin(int((m_ff1max-ff1)/m_deltaff1));
  int    ff2bin(int((m_ff2max-ff2)/m_deltaff2)); 
  double yactual(y+m_yshift);
  int    ybin(int(yactual/m_deltay));
  if (ff1bin<0 || ff1bin>m_ff1bins || 
      ff2bin<0 || ff2bin>m_ff2bins || 
      ybin<0 || ybin>m_ybins) {
    msg_Error()<<"Error in "<<METHOD<<": bins out of bounds."<<std::endl
	       <<"   b1 = "<<b1<<", b2 = "<<b2<<" --> "
	       <<"ff1 = "<<ff1<<", ff2 = "<<ff2<<", y = "<<y<<";"<<std::endl
	       <<"   ==> ff1bin = "<<ff1bin<<"("<<m_ff1bins<<"), "
	       <<"ff2bin = "<<ff2bin<<"("<<m_ff2bins<<"), "
	       <<"ybin = "<<ybin<<"("<<m_ybins<<")."<<std::endl;
    return 0.;
  }
  double ff1up(m_ff1max-ff1bin*m_deltaff1);
  double ff1low(m_ff1max-(ff1bin+1)*m_deltaff1);
  double ff2up(m_ff2max-ff2bin*m_deltaff2);
  double ff2low(m_ff2max-(ff2bin+1)*m_deltaff2);
  double yup((ybin+1)*m_deltay), ylow(ybin*m_deltay);
  double numer = 
    m_grid1[ff1bin+1][ff2bin+1][ybin+0]*
    (ff1up-ff1)*(ff2up-ff2)*(yup-yactual)   +
    m_grid1[ff1bin+1][ff2bin+0][ybin+0]*
    (ff1up-ff1)*(ff2-ff2low)*(yup-yactual)  +
    m_grid1[ff1bin+0][ff2bin+1][ybin+0]*
    (ff1-ff1low)*(ff2up-ff2)*(yup-yactual)  +
    m_grid1[ff1bin+0][ff2bin+0][ybin+0]*
    (ff1-ff1low)*(ff2-ff2low)*(yup-yactual) +
    m_grid1[ff1bin+1][ff2bin+1][ybin+1]*
    (ff1up-ff1)*(ff2up-ff2)*(yactual-ylow)  +
    m_grid1[ff1bin+1][ff2bin+0][ybin+1]*
    (ff1up-ff1)*(ff2-ff2low)*(yactual-ylow) +
    m_grid1[ff1bin+0][ff2bin+1][ybin+1]*
    (ff1-ff1low)*(ff2up-ff2)*(yactual-ylow) +
    m_grid1[ff1bin+0][ff2bin+0][ybin+1]*
    (ff1-ff1low)*(ff2-ff2low)*(yactual-ylow);
  double denom = m_deltay * m_deltaff1 * m_deltaff2;
  return numer/denom;
}

double Single_Channel_Eikonal::
Omega21(const double & b1,const double & b2,const double & y,
	const bool & plot) const 
{
  if (b1>m_b1max || b1<0. || 
      b2>m_b2max || b2<0. || 
      y>m_yshift || y<-m_yshift) return 0.;

  double ff1(p_ff1->FourierTransform(b1)), ff2(p_ff2->FourierTransform(b2));
  int    ff1bin(int((m_ff1max-ff1)/m_deltaff1));
  int    ff2bin(int((m_ff2max-ff2)/m_deltaff2)); 
  double yactual(m_yshift-y);
  int    ybin(int(yactual/m_deltay));
  if (ff1bin<0 || ff1bin>m_ff1bins || 
      ff2bin<0 || ff2bin>m_ff2bins || 
      ybin<0 || ybin>m_ybins) {
    msg_Error()<<"Error in "<<METHOD<<": bins out of bounds."<<std::endl
	       <<"   b1 = "<<b1<<", b2 = "<<b2<<" --> "
	       <<"ff1 = "<<ff1<<", ff2 = "<<ff2<<", y = "<<y<<";"<<std::endl
	       <<"   ==> ff1bin = "<<ff1bin<<"("<<m_ff1bins<<"), "
	       <<"ff2bin = "<<ff2bin<<"("<<m_ff2bins<<"), "
	       <<"ybin = "<<ybin<<")."<<std::endl;
    return 0.;
  }
  double ff1up(m_ff1max-ff1bin*m_deltaff1);
  double ff1low(m_ff1max-(ff1bin+1)*m_deltaff1);
  double ff2up(m_ff2max-ff2bin*m_deltaff2);
  double ff2low(m_ff2max-(ff2bin+1)*m_deltaff2);
  double yup((ybin+1)*m_deltay), ylow(ybin*m_deltay);
  double numer = 
    m_grid2[ff1bin+1][ff2bin+1][ybin+0]*
    (ff1up-ff1)*(ff2up-ff2)*(yup-yactual)   +
    m_grid2[ff1bin+1][ff2bin+0][ybin+0]*
    (ff1up-ff1)*(ff2-ff2low)*(yup-yactual)  +
    m_grid2[ff1bin+0][ff2bin+1][ybin+0]*
    (ff1-ff1low)*(ff2up-ff2)*(yup-yactual)  +
    m_grid2[ff1bin+0][ff2bin+0][ybin+0]*
    (ff1-ff1low)*(ff2-ff2low)*(yup-yactual) +
    m_grid2[ff1bin+1][ff2bin+1][ybin+1]*
    (ff1up-ff1)*(ff2up-ff2)*(yactual-ylow)  +
    m_grid2[ff1bin+1][ff2bin+0][ybin+1]*
    (ff1up-ff1)*(ff2-ff2low)*(yactual-ylow) +
    m_grid2[ff1bin+0][ff2bin+1][ybin+1]*
    (ff1-ff1low)*(ff2up-ff2)*(yactual-ylow) +
    m_grid2[ff1bin+0][ff2bin+0][ybin+1]*
    (ff1-ff1low)*(ff2-ff2low)*(yactual-ylow);
  double denom = m_deltay * m_deltaff1 * m_deltaff2;

  return numer/denom;
}

double Single_Channel_Eikonal::operator()(double B) {
  if (B<0. || B>=m_Bmax) return 0.;
  int Bbin(int(B/m_deltaB));
  return (m_gridB[Bbin]*((Bbin+1)*m_deltaB-B)+
	  m_gridB[Bbin+1]*(B-Bbin*m_deltaB))/m_deltaB;
}


double Single_Channel_Eikonal::Convolution2D::operator()(double b1) {
  p_convolution1D->Setb1(b1);
  return 2.*b1*p_integrator->Integrate(0.,M_PI,m_accu,1);
}

double Single_Channel_Eikonal::Convolution2D::Convolution1D::
operator()(double theta1) {
  double b2((m_b==0.)?m_b1:sqrt(m_b*m_b+m_b1*m_b1-2.*m_b*m_b1*cos(theta1)));
  double eik12(p_eikonal->Omega12(m_b1,b2,m_y));
  double eik21(p_eikonal->Omega21(m_b1,b2,m_y));
  double value(eik12*eik21);
  if (value>p_eikonal->MaxConv()) p_eikonal->SetMaxConv(value);
  return value;
}


void Single_Channel_Eikonal::PrintOmega_ik() {
  double b1(0.), b2(0.), y(-m_Y);
  double deltab1(3.), deltab2(3.), deltay(0.1);
  while (b1<3.) {
    while (b2<3.) {
      std::cout<<"Omega_ik for b1 = "<<b1<<" b2 = "<<b2<<"."<<std::endl;
      y = -m_Y;
      while (y<m_Y) {
	std::cout<<" "<<y<<"  "<<Omega12(b1,b2,y)<<std::endl;
	y += deltay;
      }
      b2 += deltab2;
    }
    b1 += deltab1;
  }
}

void Single_Channel_Eikonal::TestDEQSolvers() {
  double deltab1(0.1), deltab2(0.1);
  if (m_deqmode!=deqmode::RungeKutta4Transformed) {
    msg_Out()<<"In "<<METHOD<<"(RK2,4) : Y = "<<m_Y<<", "
	     <<"prefactor = "<<Prefactor()<<"."<<std::endl
	     <<"   Omega_{1(2)}(0,b2max,0) = "
	     <<Omega12(0.,0.999*m_b2max,-m_yshift)<<", "
	     <<"Omega_{2(1)}(0,b2max,Y) = "
	     <<Omega21(0.,0.999*m_b2max,m_yshift)<<"."<<std::endl;
    std::ofstream was;
    int i(0), j(0);
    std::string filename;
    double b1, b2, omega12, omega21, y;
    for (j=0;j<=15;j+=3) {
      //i = j;
      filename = std::string("Omega_"+ATOOLS::ToString(i)+"(")+
	ATOOLS::ToString(j)+std::string(").dat");
      was.open(filename.c_str());
      b1 = i*deltab1;
      b2 = j*deltab2;
      msg_Out()<<"   Init system for {i, j} = {"<<i<<", "<<j<<"}"
	       <<" --> b1 = "<<(i*deltab1)<<" b2 = "<<(j*deltab2)<<"."
	       <<std::endl;
      for (int ystep=0;ystep<m_ybins;ystep++) {
	y = -m_yshift+ystep*m_deltay;
	omega12 = Omega12(b1,b2,y);
	omega21 = Omega21(b1,b2,y);
	msg_Out()<<y
		 <<"   "<<omega12<<"  "<<m_grid1[i][j][ystep]
		 <<"   "<<omega21<<"  "<<m_grid2[i][j][m_ybins-ystep]<<std::endl;
	was<<y
	   <<"   "<<omega12<<"  "<<m_grid1[i][j][ystep]
	   <<"   "<<omega21<<"  "<<m_grid2[i][j][m_ybins-ystep]<<std::endl;
      }
      msg_Out()<<"   wrote out Omega_{1(2)}(0,0,y) and Omega_{2(1)}(0,0,y) "
	       <<"in Omega.dat."<<std::endl;
      was.close();
    }    
  }
  else {
    msg_Out()<<METHOD<<"(RK4T) : Y = "<<m_Y<<", "
	     <<"prefactor = "<<Prefactor()<<"."<<std::endl
	     <<"   Omega_{1(2)}(0,b2max,0) = "
	     <<Omega12(0.,0.999*m_b2max,-m_yshift)<<", "
	     <<"Omega_{2(1)}(0,b2max,Y) = "
	     <<Omega21(0.,0.999*m_b2max,m_yshift)<<"."<<std::endl;
    std::ofstream was;
    int i(0), j(0), step(1);
    std::string filename;
    double b1, b2, y;
    for (j=0;j<=15;j+=step) {
      i = j;
      filename = 
	std::string("RK4Prime_Omega_"+ATOOLS::ToString(i)+"(")+
	ATOOLS::ToString(j)+std::string(").dat");
      was.open(filename.c_str());
      b1 = i*deltab1;
      b2 = j*deltab2;
      msg_Out()<<"   Init system for {i, j} = {"<<i<<", "<<j<<"}"
	       <<" --> b1 = "<<(i*deltab1)<<" b2 = "<<(j*deltab2)<<"."
	       <<std::endl;
      for (int ystep=0;ystep<m_ybins;ystep++) {
	y = -m_yshift+ystep*m_deltay;
	was<<y<<"   "<<Omega12(b1,b2,y)<<"  "<<m_grid1[i][j][ystep]<<std::endl;
      }
      msg_Out()<<"   wrote out barOmega_{1(2)}(0,0,y) in Omega.dat."<<std::endl;
      was.close();
      if (i==3) step = 3;
    }
  }
}


void Single_Channel_Eikonal::
TestSingleEikonal(const double & b1,const double & b2) {
  msg_Tracking()<<METHOD<<": "
		<<"check grid with exact result in each ybin "
		<<"and between the bins;"<<std::endl
		<<"   fix b1 = "<<b1<<" and b2 = "<<b2<<"."<<std::endl;
  double beta02(m_beta2), Lambda2(p_ff1->Lambda2()), kappa(p_ff1->Kappa());
  double pref(beta02*Lambda2/(4.*M_PI)*exp(-sqr(b1)*Lambda2/(4.*(1.+kappa))));
  double ana, num, y;
  msg_Tracking()<<"   prefactor = "<<pref<<" from beta^2 = "<<beta02<<", "
		<<"Lambda^2 = "<<Lambda2<<", "
		<<"and kappa = "<<kappa<<std::endl
		<<"   (Note: In the testing approximation lambda=0 "
		<<"and the two Omega's decouple.)"<<std::endl; 
  std::ofstream was;
  was.open("Omega_Analytical.dat");
  for (int ystep=0;ystep<2*m_ybins;ystep++) {
    y   = -m_yshift+ystep*m_deltay/2.;
    ana = pref*exp(m_alpha*(y+m_yshift));
    num = Omega12(b1,b2,y);
    msg_Tracking()<<"y = "<<y<<"  ana = "<<ana<<", num = "<<num<<", "
		  <<"delta = "<<dabs(1.-num/ana)<<"."<<std::endl; 
    was<<y<<"  "<<ana<<"  "<<num<<"  "<<dabs(1.-num/ana)<<std::endl; 
  }
  was.close();
  exit(1);
}

