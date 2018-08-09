#include "CFPSHOWER++/Shower/Kernel.H"
#define COMPILE__Getter_Function
#define PARAMETER_TYPE CFPSHOWER::Kernel_Info
#define OBJECT_TYPE    CFPSHOWER::Kernel
#define SORT_CRITERION std::less<std::string>
#include "ATOOLS/Org/Getter_Function.C"
#include "ATOOLS/Org/Message.H"
#include <algorithm>

using namespace CFPSHOWER;
using namespace ATOOLS;
using namespace std;

Kernel::Kernel(const Kernel_Info & info) :
  m_type(info.Type()), m_flavs(info.GetFlavs()), m_swapped(info.Swapped()) {}

double Kernel::Integral(Splitting & split,const Mass_Selector * msel) {
  split.SetKernel(this);
  //if (!split.InitKinematics(msel)) return 0.;
  split.InitKinematics(msel);
  double I = (p_gauge->Charge() * p_sfunction->Integral(split) *
	      p_gauge->OverEstimate(split));
  //msg_Out()<<"   * "<<p_sfunction->Name()<<" ("<<m_swapped<<") "
  //	   <<"("<<m_flavs[0]<<" -> "<<m_flavs[1]<<" "<<m_flavs[2]<<"): "
  //	   <<p_gauge->Charge()<<" * "<<p_sfunction->Integral(split)<<" * "
  //	   <<p_gauge->OverEstimate(split)<<" = "<<I<<" / (2 pi).\n";
  return I/(2.*M_PI);
}

bool Kernel::Generate(Splitting & split,const Mass_Selector * msel,
		      const double & overfac) {
  // Sequence of generation of one splitting after kernels is selected:
  // - attach pointer of kernel to Splitting, init kinematics (setting masses,
  //   s_ijk and Q^2, making sure splitting is allowed when splitting is massive).
  // - generate a splitting kinematics, essentially the z and phi
  // - calculate x and y
  // - construct outgoing four-momenta of splitting products and spectator
  // - if this is successful, calculate weights, apply hit-or-miss, attach
  //   acceptance and rejection weights to the splitter-spectator pair
  split.SetKernel(this);
  if (!split.InitKinematics(msel)) {
    //msg_Out()<<"   * rejected (t = "<<sqrt(split.T())<<", no phase space for "
    //	     <<p_sfunction->Name()<<"(F = "<<m_flavs[1]<<", swap = "<<m_swapped<<", "
    //	     <<"Q^2 = "<<split.Q2()<<").\n";
    return false;
  }
  p_sfunction->GeneratePoint(split);
  p_sfunction->InitKinematics(split);
  if (!p_gauge->SetColours(split)) return false;
  if (p_sfunction->Construct(split)==1) {
    if (split.GetWeight()) { delete split.GetWeight(); } 
    split.SetWeight(MakeWeight(split,overfac));
    if ((*split.GetWeight())()>=ran->Get()) {
      //msg_Out()<<"   * add acceptance weight (t = "<<sqrt(split.T())<<", z = "<<split.Z()<<") "
      //       <<"from "<<(*split.GetWeight())<<" = "<<split.GetWeight()->Accept()<<"\n";
      split.GetSplitter()->AddWeight(split,true);
      return true;
    }
    else {
      //msg_Out()<<"   * add rejection weight (t = "<<sqrt(split.T())<<", z = "<<split.Z()<<") "
      //       <<"from "<<(*split.GetWeight())<<" = "<<split.GetWeight()->Reject()<<"\n";
      split.GetSplitter()->AddWeight(split,false);
    }
  }
  //else {
  //msg_Out()<<"   * rejected (t = "<<sqrt(split.T())<<", z = "<<split.Z()<<"), "
  //	   <<"no kinematics found for "<<p_sfunction->Name()<<"(swap = "<<m_swapped<<").\n";
  //}
  return false;
}

Weight * Kernel::MakeWeight(const Splitting & split,const double & overfac) {
  double weight   = (p_gauge->Charge() * (*p_gauge)(split) *
		     (*p_sfunction)(split) *
		     p_sfunction->Jacobean(split));
  double realover = (p_gauge->Charge() * p_gauge->OverEstimate(split) *
		     p_sfunction->OverEstimate(split));
  double over     = (dabs(weight)<realover) ? (weight>0.?1.:-1.)*realover : overfac*weight;
  //msg_Out()<<" *** "<<METHOD<<"("<<weight<<", "<<over<<", "<<realover<<").\n";
  return new Weight(weight,over,realover);
}

bool Kernel::FillOffsprings(Splitting & split) {
  int beam = split.GetSplitter()->Beam();
  for (size_t i=1;i<m_flavs.size();i++) {
    Parton * parton = new Parton(m_flavs[i],split.Mom(i-1));
    if (i==1) parton->SetBeam(beam);
    parton->SetColor(split.Col(i-1));
    split.SetParton(i-1,parton);
  }
  return true;
}



ostream & CFPSHOWER::operator<<(ostream &s,const Kernel & kernel) {
  s<<"Kernel["<<(kernel.GetSF()->Name())<<" / "<<(kernel.GetGauge()->Name())<<"].\n";
  return s;
}

DECLARE_GETTER(Kernel,"Kernel",Kernel,Kernel_Info);

Kernel * Getter<Kernel,Kernel_Info,Kernel>::operator()(const Parameter_Type & info) const
{
  SF_Base    * sfunction = SF_Getter::GetObject(info.SFName(),info);
  Gauge_Base * gaugepart = GP_Getter::GetObject(info.GPName(),info);
  if (!sfunction || !gaugepart) {
    if (sfunction) delete sfunction;
    if (gaugepart) delete gaugepart;
    return NULL;
  }
  Kernel * kernel = new Kernel(info);
  kernel->SetSF(sfunction);
  kernel->SetGauge(gaugepart);
  //msg_Out()<<" *** Init new Kernel("<<info.GetFlavs()[0]<<" -> "
  //	   <<info.GetFlavs()[1]<<" "<<info.GetFlavs()[2]<<"): "
  //	   <<sfunction->Name()<<" for swap = "<<info.Swapped()<<".\n";
  return kernel;
}

void Getter<Kernel,Kernel_Info,Kernel>::PrintInfo(ostream &str,const size_t width) const
{
  str<<"Splitting Kernel";
}
