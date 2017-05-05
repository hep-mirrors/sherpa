#include "DIRE/Shower/Kernel.H"

#include "DIRE/Tools/Parton.H"
#include "DIRE/Tools/Amplitude.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

#include <typeinfo>

using namespace DIRE;
using namespace ATOOLS;

Kernel::Kernel(Shower *const ps,Kernel_Key key):
  p_ps(ps), p_lf(NULL), p_gf(NULL), m_ef(1.0),
  m_type(key.m_type), m_mode(key.m_mode), m_on(1)
{
  key.p_k=this;
  std::string gauge;
  if (key.p_v) {
    if (key.p_v->order[0]) {
      gauge="QCD";
      gauge+="{"+ToString(key.p_v->in[0].StrongCharge())+"}";
      gauge+="{"+ToString(key.p_v->in[1+key.m_mode].StrongCharge())+"}";
      gauge+="{"+ToString(key.p_v->in[2-key.m_mode].StrongCharge())+"}";
    }
  }
  else {
    gauge="QCD";
    for (size_t i(0);i<key.m_fl.size();++i)
      gauge+="{"+ToString(key.m_fl[i].StrongCharge())+"}";
  }
  p_gf = Gauge_Getter::GetObject(gauge,key);
  if (p_gf==NULL) {
    m_on=-1;
    return;
  }
  std::string type((key.m_type&1)?"I":"F");
  type+=(key.m_type&2)?"I":"F";
  p_lf = Lorentz_Getter::GetObject
    (type+"_"+(key.p_v?key.p_v->Lorentz[0]:key.m_lfid),key);
  if (p_lf==NULL) {
    m_on=-1;
    return;
  }
  std::string ffl(ToString(p_lf->Flav(1)));
  for (size_t i(2);i<p_lf->Flavs().size();++i)
    ffl+=","+ToString(p_lf->Flav(i));
  msg_Debugging()<<"Init("<<m_on<<") "<<p_lf->Flav(0)<<"->"<<ffl
		 <<" => ("<<Demangle(typeid(*p_lf).name()).substr(6)
		 <<","<<Demangle(typeid(*p_gf).name()).substr(6)
		 <<"), mode "<<key.m_mode<<", swap "<<key.m_swap<<"\n";
}

Kernel::~Kernel()
{
  if (p_gf) delete p_gf;
  if (p_lf) delete p_lf;
}

double Kernel::Value(const Splitting &s) const
{
  return p_gf->Value(s)*p_lf->Value(s)*p_lf->Jacobian(s);
}

Weight Kernel::GetWeight(const Splitting &s,const double &o) const
{
  double f(p_gf->Value(s)*p_lf->Value(s)*p_lf->Jacobian(s));
  double h(p_gf->Estimate(s)*p_lf->Estimate(s));
  double g(dabs(f)<h?(f>=0.0?h:-h):o*f);
  return Weight(f,g,m_ef*h);
}

bool Kernel::GeneratePoint(Splitting &s) const
{
  s.p_sk=this;
  if (!p_lf->SetLimits(s)) return false;
  if (!p_lf->GeneratePoint(s)) return false;
  return p_gf->GeneratePoint(s);
}

double Kernel::Integral(Splitting &s) const
{
  if (!p_lf->SetLimits(s)) return false;
  return m_ef*p_lf->Integral(s)*p_gf->Estimate(s);
}

int Kernel::Construct(Splitting &s,const int mode) const
{
  int stat=p_lf->Construct(s,mode);
  if (stat!=1) return stat;
  return p_gf->Construct(s);
}
