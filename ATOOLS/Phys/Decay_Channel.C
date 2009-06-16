#include "ATOOLS/Phys/Decay_Channel.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Mass_Handler.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"

using namespace ATOOLS;
using namespace std;

Decay_Channel::Decay_Channel(const Flavour & _flin) :
  m_width(0.), m_deltawidth(-1.), m_minmass(0.), m_max(0.), m_flin(_flin)
 { }

Decay_Channel::Decay_Channel(const Decay_Channel & _dec) :
  m_width(_dec.m_width), m_deltawidth(_dec.m_deltawidth), 
  m_minmass(_dec.m_minmass),
  m_flin(_dec.m_flin), m_flouts(_dec.m_flouts) { }

Decay_Channel::~Decay_Channel()
{
  // have to find solution to delete multichannel
}

void Decay_Channel::Output() const
{
  msg_Out()<<(*this);
}

std::ostream &ATOOLS::operator<<(std::ostream &os,const Decay_Channel &dc)
{
  os<<dc.m_flin<<" -> ";
  for (FlSetConstIter fl=dc.m_flouts.begin();fl!=dc.m_flouts.end();++fl) 
    os<<(*fl)<<" ";
  os<<" : "<<dc.m_width;
  if (dc.m_deltawidth>=0.) os<<" (+/-"<<dc.m_deltawidth<<")";
  os<<" GeV";
  os<<"."<<endl;
  return os;
}

string Decay_Channel::Name() const
{
  string name=m_flin.IDName()+string(" --> ");
  for (FlSetConstIter flit=m_flouts.begin();flit!=m_flouts.end();++flit) {
    name+=flit->IDName()+string(" ");
  }
  return name;
}

double DCLambda(double a, double b, double c)
{
  double L = (sqr(a-b-c)-4.*b*c);
  if (L>0.0) return sqrt(L)/2/sqrt(a);
  if (L>-Accu()) return 0.0;
  msg_Error()<<"passed impossible mass combination:"<<std::endl;
  msg_Error()<<"m_a="<<sqrt(a)<<" m_b="<<sqrt(b)<<" m_c="<<sqrt(c)<<endl;
  msg_Error()<<"L="<<L<<endl;
  return 0.;
}

double DCWeight(double s, double sp, double b, double c)
{
  return DCLambda(sp,b,c)/DCLambda(s,b,c)*s/sp;
}

double Decay_Channel::DiceMass(const double& min,const double& max) const
{
  double mass=-1.0;
  double decaymin = MinimalMass();
  DEBUG_VAR(decaymin);
  Mass_Handler masshandler(GetDecaying());
  if(decaymin>max) mass=-1.0;
  else if (decaymin==0.0) mass = masshandler.GetMass(decaymin, max);
  else {
    double s=sqr(GetDecaying().HadMass());
    double mb(0.0), mc(0.0);
    for (int i=0; i<NOut(); ++i) {
      mc+=GetDecayProduct(i).HadMass();
      if(GetDecayProduct(i).HadMass()>mb)
        mb=GetDecayProduct(i).HadMass();
    }
    mc-=mb;
    double b=sqr(mb);
    double c=sqr(mc);
    double spmax=2.0*b+2.0*c+sqrt(sqr(b)+14.0*b*c+sqr(c));
    double wmax=DCWeight(s,spmax,b,c);
    double w=0.0;
    int trials(0);
    do {
      mass = masshandler.GetMass(decaymin, max);
      double sp=sqr(mass);
      w=DCWeight(s,sp,b,c);
      ++trials;
      if (w>wmax+Accu())
        msg_Error()<<METHOD<<" w="<<w<<" > wmax="<<wmax<<std::endl;
    } while (w<ran.Get()*wmax && trials<1000);
  }
  return mass;
}
