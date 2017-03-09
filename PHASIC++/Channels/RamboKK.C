#include "PHASIC++/Channels/RamboKK.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Process/Process_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/My_MPI.H"
#include "MODEL/Main/Model_Base.H"

using namespace PHASIC;
using namespace ATOOLS;
using namespace std;

RamboKK::RamboKK(int _nin,int _nout,const Flavour * fl)// : nin(_nin), nout(_nout)
{
  nin=_nin;nout=_nout;
  xm2 = new double[nin+nout];
  p2  = new double[nin+nout];  
  E   = new double[nin+nout];
  ms  = new double[nin+nout];
  rans= 0;
  rannum=0;
  massflag = 0;
  for (short int i=0;i<nin+nout;i++) {
    ms[i] = ATOOLS::sqr(fl[i].Mass());
    if (!ATOOLS::IsZero(ms[i])) massflag = 1;
  } 

  double   pi2log = log(M_PI/2.);
  double * Z      = new double[nout+1];
  Z[2] = pi2log;
  for (short int k=3;k<=nout;k++) Z[k] = Z[k-1]+pi2log-2.*log(double(k-2));
  for (short int k=3;k<=nout;k++) Z[k] = Z[k]-log(double(k-1));
  Z_N  = Z[nout];
  delete[] Z;

  kkp=-1;mpss=1.;
  int mode = MODEL::s_model->ScalarNumber(std::string("KK_mode"));
  for (int i=nin;i<nin+nout;i++) {
    if(fl[i].IsKK() && (mode==1 || mode==2 || mode==5)){
      if(ATOOLS::IsZero(ms[i])){
	msg_Error()<<"Error in RamboKK: "<<endl
		   <<"   Please initialize with nonzero particle mass ("<<fl[i]<<") !"<<std::endl;
	Abort();
      }
      kkp=i;
      
      ed  = MODEL::s_model->ScalarNumber(std::string("ED"));
      r2  = sqr(MODEL::s_model->ScalarConstant(std::string("Radius")));
      gn  = MODEL::s_model->ScalarConstant(std::string("G_Newton"));

      //Calculation of Gamma(ed/2)
      if(ed%2==0) gam=1.;
      else gam=sqrt(M_PI);
      for(int k=2-ed%2;k<ed;k+=2)gam*=0.5*k;

      double mm=rpa->gen.Ecms();
      prevET = mm;
      for(int j=nin;j<nin+nout;j++)
	if(j!=i) mm -= sqrt(ms[j]);
      maxm2=sqr(mm);
      maxn=sqrt(maxm2*r2);
      mpss=2.*pow(sqrt(M_PI)*maxn,double(ed))/ed/gam;
      break;
    }
  }
}

RamboKK::~RamboKK() 
{
  if (xm2) { delete [] xm2; xm2 = 0; }
  if (p2)  { delete [] p2;  p2  = 0; }
  if (E)   { delete [] E;   E   = 0; }
}

void RamboKK::Set_KKmass()
{
  if (kkp==-1) return;
 
  double * nv = new double[ed];
  double ms2;
  do{
    ms2=0;
    for (short int i=0;i<ed;i++) {
      nv[i]=ran->Get()*maxn;
      ms2+=sqr(nv[i]);
    }
    ms2*=4*sqr(M_PI)/r2;
  }
  while (ms2>maxm2);
  ms[kkp]=ms2;
  delete[] nv;
}


void RamboKK::GenerateWeight(Vec4D * p,Cut_Data * cuts)
{
Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<nin;i++) sump += p[i];
  double ET = sqrt(sump.Abs2());
  weight    = 1.;
  if (massflag) MassiveWeight(p,ET);
  weight   *= exp((2.*nout-4.)*log(ET)+Z_N)/pow(2.*M_PI,nout*3.-4.);
  weight   *= mpss;
}

void RamboKK::GeneratePoint(Vec4D * p,Cut_Data * cuts)
{  

  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<nin;i++) sump += p[i];

  double ET = sqrt(sump.Abs2());

  if (!IsEqual(ET,prevET) && kkp>-1) {
    double mm = prevET = ET;
    for(int j=nin;j<nin+nout;j++)
      if(j!=kkp) mm -= sqrt(ms[j]);
    maxm2=sqr(mm);
    maxn=sqrt(maxm2*r2);
    mpss=2.*pow(sqrt(M_PI)*maxn,double(ed))/ed/gam;
  }

  Set_KKmass();
  
  double Q, S, C, F, G, A, X, RMAS, BQ, e;
  short int i;
  Vec4D R;
  Vec3D B;
  
  for(i=nin;i<nin+nout;i++) {
    C     = 2*ran->Get()-1;
    S     = sqrt(1-C*C);
    F     = 2*M_PI*ran->Get();
    Q     = -log(ran->Get()*ran->Get());
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
  // if (massflag) 
  MassivePoint(p,ET);// The boost is numerically not very precise, MassivePoint is always called for momentum conservation
}

void RamboKK::GeneratePoint(Vec4D * p,Cut_Data * cuts,double * _ran) {
  GeneratePoint(p,cuts);
}

void RamboKK::MassiveWeight(Vec4D* p,double ET)
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

void RamboKK::MassivePoint(Vec4D* p,double ET)
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

namespace PHASIC {

  class RamboKK_Channel_Generator: public Channel_Generator {
  public:
    
    RamboKK_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new RamboKK(p_proc->NIn(),p_proc->NOut(),
			  &p_proc->Flavours().front()));
      return 0;
    }

  };// end of class RamboKK_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(RamboKK_Channel_Generator,"RamboKK",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,RamboKK_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new RamboKK_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    RamboKK_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Rambo integrator for KK states";
}
