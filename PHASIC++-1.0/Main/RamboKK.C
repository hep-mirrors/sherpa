#include "RamboKK.H"
#include "Message.H"
#include "Random.H"
#include "Run_Parameter.H"
#include "Couplings_LED.H"

using namespace PHASIC;
using namespace AORGTOOLS;
using namespace APHYTOOLS;
using namespace AMATOOLS;
using namespace AMEGIC;

RamboKK::RamboKK(int _Nin,int _Nout,Flavour* fl) : Nin(_Nin), Nout(_Nout) 
{
  double pi2log =log(M_PI/2.);
  
  double* Z = new double[Nout+1];

  Z[2]=pi2log;
  for (short int k=3;k<=Nout;k++) Z[k] = Z[k-1]+pi2log-2.*log(double(k-2));
  for (short int k=3;k<=Nout;k++) Z[k] = Z[k]-log(double(k-1));

  Z_N = Z[Nout];
  delete[] Z;
  
  xm2 = new double[Nin+Nout];
  p2  = new double[Nin+Nout];  
  E   = new double[Nin+Nout];
  ms_out = new double[Nin+Nout];

  massflag = 0;
  
  for (short int i=0;i<Nin+Nout;i++) {
    ms_out[i] = sqr(fl[i].Mass());
    //cout<<"Masses in Rambo: "<<sqrt(ms_out[i])<<endl;
    if (!AMATOOLS::IsZero(ms_out[i])) massflag = 1;
  } 

  kkp=-1;mpss=1.;
  for (short int i=Nin;i<Nin+Nout;i++) {
    if(fl[i].IsKK()){
      if(AMATOOLS::IsZero(ms_out[i])){
	cout<<"Please initialize with nonzero particle mass ("<<fl[i]<<") !"<<endl;
	abort();
      }
      kkp=i;
      Couplings_LED  CplLED;
      CplLED.Init();
      ed=CplLED.Ned();
      R2=sqr(CplLED.R());
      double mm=rpa.gen.Ecms();
      for(short int j=Nin;j<Nin+Nout;j++)
	if(j!=i)mm-=sqrt(ms_out[j]);
      maxM2=sqr(mm);
      maxN=sqrt(maxM2*R2/4./sqr(M_PI));
      mpss=1./ed*pow(maxM2,0.5*(double(ed)))/pow(CplLED.Ms(),2.+(double(ed)))/CplLED.Gn();
      break;
    }
  }

}

void RamboKK::Set_KKmass()
{
  if(kkp==-1)return;
  double *nv=new double[ed];
  double ms2;
  do{
    ms2=0;
    //cout<<"Set_KKmass: "<<kkp<<","<<maxM2<<","<<maxN<<":";
    for (short int i=0;i<ed;i++) {
      nv[i]=ran.Get()*maxN;
      //cout<<nv[i]<<",";
      ms2+=sqr(nv[i]);
    }
    ms2*=4*sqr(M_PI)/R2;
    //cout<<ms2<<endl;
  }while (ms2>maxM2);
  ms_out[kkp]=ms2;
  delete[] nv;
}

void RamboKK::Generate_Weight(Vec4D* p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<Nin;i++) sump += p[i];
  double ET = sqrt(sump.Abs2());

  weight = 1.;
  if (massflag) Massive_Weight(p,ET);
  weight *= ::exp((2.*Nout-4.)*log(ET)+Z_N)/pow(2.*M_PI,Nout*3.-4.);
  weight *=mpss;
}

void RamboKK::Generate_Point(Vec4D* p,Cut_Data * cuts)
{

  Set_KKmass();
  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<Nin;i++) sump += p[i];
  double ET = sqrt(sump.Abs2());
  
  double Q, S, C, F, G, A, X, RMAS, BQ, e;
  short int i;
  Vec4D R;
  Vec3D B;
  
  for(i=Nin;i<Nin+Nout;i++) {
    C     = 2*ran.Get()-1;
    S     = sqrt(1-C*C);
    F     = 2*M_PI*ran.Get();
    Q     = -log(ran.Get()*ran.Get());
    p[i]  = Vec4D(Q, Q*S*::sin(F), Q*S*cos(F), Q*C);
    R    += p[i]; 
  }

  RMAS = sqrt(R.Abs2());
  B    = (-1)*Vec3D(R)/RMAS;
  G    = R[0]/RMAS;
  A    = 1.0/(1.0+G);
  X    = ET/RMAS;
  
  for(i=Nin;i<Nin+Nout;i++) {
    e     = p[i][0];
    BQ    = B*Vec3D(p[i]);
    p[i]  = X*Vec4D((G*e+BQ),Vec3D(p[i])+B*(e+A*BQ));
  }

  weight = 1.;
  if (massflag) Massive_Point(p,ET);
  //for (short int i=0;i<Nin+Nout;i++) cout<<i<<". Momentum: "<<p[i]<<";"<<sqrt(p[i].Abs2())<<endl;
}

void RamboKK::Generate_Point(Vec4D* p,Cut_Data * cuts,double* ran)
{
  msg.Error()<<"In Rambo a dummy routine!"<<endl;
  Generate_Point(p,cuts);
}

void RamboKK::Massive_Weight(Vec4D* p,double ET)
{
  //Maxnumber of Iterations
  short int itmax = 6;
  //Accuracy
  double accu = ET*pow(10.,-14.);

  double xmt = 0.;
  double x;
 
  for (short int i=Nin;i<Nin+Nout;i++) {
    xm2[i] = 0.;
    xmt   += sqrt(ms_out[i]);
    p2[i]  = sqr(Vec3D(p[i]).Abs());
  }
  x = 1./sqrt(1.-sqr(xmt/ET));
  xmt = 0.;
  //
  // MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X    
  double f0,g0,x2;
    
  short int iter = 0; 
  for (;;) {
    f0=-ET;
    g0=0.;
    x2=x*x;
    for (short int i=Nin;i<Nin+Nout;i++) {
      E[i] = sqrt(xm2[i]+x2*p2[i]);
      f0 += E[i];
      g0 += p2[i]/E[i];
    }
    if (dabs(f0)<accu) break; 
    iter++;
    if (iter>itmax) break;
    x-=f0/(x*g0);  
  }
  
  double wt2=1.;
  double wt3=0.;
  double v;
  
  // Calculate Momenta + Weight 
  for (short int i=Nin;i<Nin+Nout;i++) {
    v    = Vec3D(p[i]).Abs();
    wt2 *= v/p[i][0];
    wt3 += v*v/p[i][0];
  }  
  x = 1./x;
  weight = ::exp((2.*Nout-3.)*log(x)+log(wt2/wt3*ET));
}

void RamboKK::Massive_Point(Vec4D* p,double ET)
{
  //Maxnumber of Iterations
  short int itmax = 6;
  //Accuracy
  double accu = ET*pow(10.,-14.);

  double xmt = 0.;
  double x;
 
  for (short int i=Nin;i<Nin+Nout;i++) {
    xmt   += sqrt(ms_out[i]);
    xm2[i] = ms_out[i];
    p2[i]  = sqr(p[i][0]);
  }

  x = sqrt(1.-sqr(xmt/ET));
  //
  // MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X    
  double f0,g0,x2;
    
  short int iter = 0; 
  for (;;) {
    f0=-ET;
    g0=0.;
    x2=x*x;
    for (short int i=Nin;i<Nin+Nout;i++) {
      E[i] = sqrt(xm2[i]+x2*p2[i]);
      f0 += E[i];
      g0 += p2[i]/E[i];
    }
    if (dabs(f0)<accu) break; 
    iter++;
    if (iter>itmax) break;
    x-=f0/(x*g0);  
  }

  //double wt2=1.;
  //double wt3=0.;
  //double v;
  
  // Calculate Momenta + Weight 
  for (short int i=Nin;i<Nin+Nout;i++) {
    //v = x*p[i][0];
    //wt2 *= v/E[i];
    //wt3 += v*v/E[i];
    p[i] = Vec4D(E[i],x*Vec3D(p[i]));
  }  

  //weight = ::exp((2.*Nout-3.)*log(x)+log(wt2/wt3*ET));
}







