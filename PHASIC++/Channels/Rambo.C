#include "PHASIC++/Channels/Rambo.H"
#include "PHASIC++/Channels/Channel_Generator.H"
#include "PHASIC++/Channels/Multi_Channel.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"

using namespace PHASIC;
using namespace ATOOLS;


Rambo::Rambo(size_t _nin,size_t _nout,const Flavour * fl, const Mass_Selector* _ms) :
  p_masssel(_ms)
{
  std::vector<double> masses(_nin+_nout, 0.0);
  for (short int i=0;i<_nin+_nout;i++) {
    masses[i]=0.0;
    for (size_t j(0);j<fl[i].Size();++j) 
      masses[i]+=p_masssel ? p_masssel->Mass(fl[i][j]) : fl[i][j].Mass();
    masses[i]/=fl[i].Size();
  }
  Init(_nin,_nout, masses);
}

Rambo::Rambo(size_t _nin,std::vector<double> masses) :
  p_masssel(NULL)
{
  Init(_nin, masses.size()-_nin, masses);
}

void Rambo::Init(const size_t& _nin,
		 const size_t& _nout,
		 const std::vector<double>& masses)
{
  m_nin    = _nin;
  m_nout   = _nout;
  xm2      = new double[m_nin+m_nout+1];
  p2       = new double[m_nin+m_nout+1];  
  E        = new double[m_nin+m_nout+1];
  p_ms     = new double[m_nin+m_nout+1];
  p_rans   = 0;
  m_rannum = 0;
  massflag = 0;

  for (short int i=0;i<m_nin+m_nout;i++)
    {
      p_ms[i]=sqr(masses[i]);
      if (!ATOOLS::IsZero(p_ms[i])) massflag = 1;
    } 

  double   pi2log = log(M_PI/2.);
  double * Z      = new double[m_nout+1];
  Z[2] = pi2log;
  for (short int k=3;k<=m_nout;k++) Z[k] = Z[k-1]+pi2log-2.*log(double(k-2));
  for (short int k=3;k<=m_nout;k++) Z[k] = Z[k]-log(double(k-1));
  Z_N  = Z[m_nout];
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
  for (short int i=0;i<m_nin;i++) sump += p[i];
  double ET = sqrt(sump.Abs2());
  m_weight    = 1.;
  if (massflag) MassiveWeight(p,ET);
  m_weight   *= exp((2.*m_nout-4.)*log(ET)+Z_N)/pow(2.*M_PI,m_nout*3.-4.);
}

ATOOLS::Vec4D_Vector Rambo::GeneratePoint(const double& E)
{
  ATOOLS::Vec4D_Vector p; p.resize(m_nin+m_nout);
  if (E<p_ms[0]+p_ms[1]) THROW(fatal_error, "sqrt(s) smaller than particle masses");
  double x=1.0/2.0+(p_ms[0]*p_ms[0]-p_ms[1]*p_ms[1])/(2.0*E*E);
  p[0]=ATOOLS::Vec4D(x*E,0.0,0.0,sqrt(ATOOLS::sqr(x*E)-p_ms[0]*p_ms[0]));
  p[1]=ATOOLS::Vec4D((1.0-x)*E,ATOOLS::Vec3D(-p[0]));
  GeneratePoint(&p[0]);
  return p;
}

void Rambo::GeneratePoint(Vec4D * p,Cut_Data * cuts)
{
  Vec4D sump(0.,0.,0.,0.);
  for (short int i=0;i<m_nin;i++) sump += p[i];

  double ET = sqrt(sump.Abs2());
  
  double Q, S, C, F, G, A, X, RMAS, BQ, e;
  short int i;
  Vec4D R;
  Vec3D B;

  for(i=m_nin;i<m_nin+m_nout;i++) {
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
  
  for(i=m_nin;i<m_nin+m_nout;i++) {
    e     = p[i][0];
    BQ    = B*Vec3D(p[i]);
    p[i]  = X*Vec4D((G*e+BQ),Vec3D(p[i])+B*(e+A*BQ));
  }

  m_weight = 1.;
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
  for (short int i=m_nin;i<m_nin+m_nout;i++) {
    xm2[i]   = 0.;
    xmt     += sqrt(p_ms[i]);
    p2[i]    = sqr(Vec3D(p[i]).Abs());
  }
  double x   = 1./sqrt(1.-sqr(xmt/ET));

  // Massive particles : Rescale their momenta by a common factor x

  // Loop to calculate x
  double f0,g0,x2;    
  short int iter = 0; 
  for (;;) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=m_nin;i<m_nin+m_nout;i++) {
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
  for (short int i=m_nin;i<m_nin+m_nout;i++) {
    v    = Vec3D(p[i]).Abs();
    wt2 *= v/p[i][0];
    wt3 += v*v/p[i][0];
  }  
  x      = 1./x;
  m_weight = exp((2.*m_nout-3.)*log(x)+log(wt2/wt3*ET));
}

void Rambo::MassivePoint(Vec4D* p,double ET)
{
  itmax = 6;
  accu  = ET * 1.e-14; //pow(10.,-14.);


  double xmt = 0.;
  double x;
 
  for (short int i=m_nin;i<m_nin+m_nout;i++) {
    xmt   += sqrt(p_ms[i]);
    xm2[i] = p_ms[i];
    p2[i]  = sqr(p[i][0]);
  }

  x = sqrt(1.-sqr(xmt/ET));

  // Massive particles : Rescale their momenta by a common factor x
    
  // Loop to calculate x

  double f0,g0,x2;
  short int iter = 0; 
  for (;;) {
    f0 = -ET;g0 = 0.;x2 = x*x;
    for (short int i=m_nin;i<m_nin+m_nout;i++) {
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
  for (short int i=m_nin;i<m_nin+m_nout;i++) p[i] = Vec4D(E[i],x*Vec3D(p[i]));
}

namespace PHASIC {

  class Rambo_Channel_Generator: public Channel_Generator {
  public:
    
    Rambo_Channel_Generator(const Channel_Generator_Key &key):
    Channel_Generator(key) {}

    int GenerateChannels()
    {
      p_mc->Add(new Rambo(p_proc->NIn(),p_proc->NOut(),
			  &p_proc->Flavours().front(),
			  p_proc->Generator()));
      return 0;
    }

  };// end of class Rambo_Channel_Generator

}// end of namespace PHASIC

DECLARE_GETTER(Rambo_Channel_Generator,"Rambo",
	       Channel_Generator,Channel_Generator_Key);

Channel_Generator *ATOOLS::Getter
<Channel_Generator,Channel_Generator_Key,Rambo_Channel_Generator>::
operator()(const Channel_Generator_Key &args) const
{
  return new Rambo_Channel_Generator(args);
}

void ATOOLS::Getter<Channel_Generator,Channel_Generator_Key,
		    Rambo_Channel_Generator>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"Rambo integrator";
}
