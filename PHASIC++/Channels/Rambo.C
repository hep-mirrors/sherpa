#include "PHASIC++/Channels/Rambo.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace PHASIC;
using namespace ATOOLS;


Rambo::Rambo(int _nin,int _nout,const Flavour * fl, const Mass_Selector* _ms) :
  p_ms(_ms)
{
  nin=_nin; nout=_nout;
  xm2 = new double[nin+nout+1];
  p2  = new double[nin+nout+1];  
  E   = new double[nin+nout+1];
  ms  = new double[nin+nout+1];
  rans= 0;
  rannum=0;
  massflag = 0;
  for (short int i=0;i<nin+nout;i++) {
    ms[i]=sqr(p_ms->Mass(fl[i]));
    if (!ATOOLS::IsZero(ms[i])) massflag = 1;
  } 

  double   pi2log = log(M_PI/2.);
  double * Z      = new double[nout+1];
  Z[2] = pi2log;
  for (short int k=3;k<=nout;k++) Z[k] = Z[k-1]+pi2log-2.*log(double(k-2));
  for (short int k=3;k<=nout;k++) Z[k] = Z[k]-log(double(k-1));
  Z_N  = Z[nout];
  delete[] Z;
}

Rambo::~Rambo() 
{
  if (xm2) { delete [] xm2; xm2 = 0; }
  if (p2)  { delete [] p2;  p2  = 0; }
  if (E)   { delete [] E;   E   = 0; }
}



void Rambo::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<nin;i++) sump += p[i];
  double ET = sqrt(sump.Abs2());
  weight    = 1.;
  if (massflag) MassiveWeight(p,ET);
  weight   *= exp((2.*nout-4.)*log(ET)+Z_N)/pow(2.*M_PI,nout*3.-4.);
}

void Rambo::GeneratePoint(Vec4D * p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<nin;i++) sump += p[i];

  double ET = sqrt(sump.Abs2());
  
  double Q, S, C, F, G, A, X, RMAS, BQ, e;
  short int i;
  Vec4D R;
  Vec3D B;
  
  for(i=nin;i<nin+nout;i++) {
    C     = 2*ran->Get()-1;
    S     = sqrt(1-C*C);
    F     = 2*M_PI*ran->Get();
    Q     = -log( std::min( 1.0-1.e-10, std::max(1.e-10,ran->Get()*ran->Get()) ) );
    p[i]  = Vec4D(Q, Q*S*::sin(F), Q*S*cos(F), Q*C);
    R    += p[i]; 
  }

  RMAS = sqrt(R.Abs2());
  B    = (-1)*Vec3D(R)/RMAS;
  G    = R[0]/RMAS;
  A    = 1.0/(1.0+G);
  X    = ET/RMAS;
  
  for(i=nin;i<nin+nout;i++) {
    e     = p[i][0];
    BQ    = B*Vec3D(p[i]);
    p[i]  = X*Vec4D((G*e+BQ),Vec3D(p[i])+B*(e+A*BQ));
  }

  weight = 1.;
  //if (massflag)
  MassivePoint(p,ET); // The boost is numerically not very precise, MassivePoint is always called for momentum conservation
}

void Rambo::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran) {
  GeneratePoint(p,cuts);
}

void Rambo::MassiveWeight(Vec4D* p,double ET)
{
  itmax = 6;
  accu  = ET * pow(10.,-14.);

  double xmt = 0.; 
  for (short int i=nin;i<nin+nout;i++) {
    xm2[i]   = 0.;
    xmt     += sqrt(ms[i]);
    p2[i]    = sqr(Vec3D(p[i]).Abs());
  }
  double x   = 1./sqrt(1.-sqr(xmt/ET));
  xmt        = 0.;

  // Massive particles : Rescale their momenta by a common factor x

  // Loop to calculate x
  double f0,g0,x2;    
  short int iter = 0; 
  for (;;) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=nin;i<nin+nout;i++) {
      E[i] = sqrt(xm2[i]+x2*p2[i]);
      f0  += E[i];
      g0  += p2[i]/E[i];
    }
    if (dabs(f0)<accu) break; 
    iter++;
    if (iter>itmax) break;
    x -= f0/(x*g0);  
  }
  
  double wt2 = 1.;
  double wt3 = 0.;
  double v;
  
  // Calculate Momenta + Weight 
  for (short int i=nin;i<nin+nout;i++) {
    v    = Vec3D(p[i]).Abs();
    wt2 *= v/p[i][0];
    wt3 += v*v/p[i][0];
  }  
  x      = 1./x;
  weight = exp((2.*nout-3.)*log(x)+log(wt2/wt3*ET));
}

void Rambo::MassivePoint(Vec4D* p,double ET)
{
  itmax = 6;
  accu  = ET * 1.e-14; //pow(10.,-14.);


  double xmt = 0.;
  double x;
 
  for (short int i=nin;i<nin+nout;i++) {
    xmt   += sqrt(ms[i]);
    xm2[i] = ms[i];
    p2[i]  = sqr(p[i][0]);
  }

  x = sqrt(1.-sqr(xmt/ET));

  // Massive particles : Rescale their momenta by a common factor x
    
  // Loop to calculate x

  double f0,g0,x2;
  short int iter = 0; 
  for (;;) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=nin;i<nin+nout;i++) {
      E[i] = sqrt(xm2[i]+x2*p2[i]);
      f0  += E[i];
      g0  += p2[i]/E[i];
    }
    if (dabs(f0)<accu) break; 
    iter++;
    if (iter>itmax) break;
    x -= f0/(x*g0);  
  }
  
  // Construct Momenta
  for (short int i=nin;i<nin+nout;i++) p[i] = Vec4D(E[i],x*Vec3D(p[i]));
}







