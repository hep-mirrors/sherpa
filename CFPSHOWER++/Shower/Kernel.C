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
  m_type(info.Type()),
  m_split(info.GetSplit()), m_flavs(info.GetFlavs()),
  p_msel(NULL)
{
  for (size_t beam=0;beam<2;beam++) p_pdf[beam] = NULL;
}

double Kernel::Integral(Splitting & split,const Mass_Selector * msel) {
  split.SetKernel(this);
  split.InitKinematics(msel);
  double SF = p_sf->Integral(split);
  double I = (p_gauge->Charge() * SF *
	      p_gauge->OverEstimate(split));
  //msg_Out()<<METHOD<<"(SF = "<<SF<<" * "
  //	   <<"gauge = "<<p_gauge->OverEstimate(split)<<" * "
  //	   <<"charge = "<<p_gauge->Charge()<<")\n";
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
  split.InitKinematics(msel);
  p_sf->GeneratePoint(split);
  if (!p_sf->InitKinematics(split)) {
    //msg_Out()<<"   * rejected (t = "<<split.t()<<", no phase space for "
    //	     <<p_sf->Name()<<"(F = "<<m_flavs[1]<<", swap = "<<m_swapped<<", "
    //	     <<"Q^2 = "<<split.Q2()<<").\n";
    return false;
  }
  if (!p_gauge->SetColours(split)) {
    //msg_Out()<<"   * no colours for (t = "<<split.t()<<", no phase space for "
    //	     <<p_sf->Name()<<"(F = "<<m_flavs[1]<<", swap = "<<m_swapped<<", "
    //	     <<"Q^2 = "<<split.Q2()<<").\n";
    return false;
  }
  if (p_sf->Construct(split)==1) {
    if (split.GetWeight()) { delete split.GetWeight(); } 
    split.SetWeight(MakeWeight(split,overfac));
    //msg_Out()<<"   * weight = "<<(p_gauge->Charge()*(*p_gauge)(split))<<"(Col)"
    //	     <<" * "<<(*p_sf)(split)<<"(SF) * "<<p_sf->Jacobean(split)<<"(J) "
    //	     <<"= "<<((*split.GetWeight())())<<"\n";
    if ((*split.GetWeight())()>=ran->Get()) {
      //msg_Out()<<"   * add acceptance weight "
      //       <<"(t = "<<split.t()<<", zi = "<<split.zi()<<", "
      //       <<"eta = "<<split.eta()<<") "
      //       <<"from "<<(*split.GetWeight())<<" = "
      //	       <<split.GetWeight()->Accept()<<"\n";
      split.GetSplitter()->AddWeight(split,true);
      return true;
    }
    else {
      //msg_Out()<<"   * add rejection weight "
      //       <<"(t = "<<split.t()<<", zi = "<<split.zi()<<", "
      //       <<"eta = "<<split.eta()<<") "
      //       <<"from "<<(*split.GetWeight())<<" = "
      //       <<split.GetWeight()->Reject()<<"\n";
      split.GetSplitter()->AddWeight(split,false);
      return false;
    }
  }
  //msg_Out()<<"   * "<<METHOD<<" no kinematics for (t = "<<split.t()<<", "
  //	   <<"no phase space for "<<p_sf->Name()
  //	   <<"(F = "<<m_flavs[1]<<", swap = "<<m_swapped<<", "
  //	   <<"Q^2 = "<<split.Q2()<<").\n";
  return false;
}

double Kernel::GetXPDF(const double & x,const double & Q2,
		       const ATOOLS::Flavour & flav,const size_t beam) {
  PDF::PDF_Base * pdf = p_pdf[beam];
  if (pdf==NULL)                          return 0.;
  if (x<pdf->XMin()   || x>pdf->XMax() ||
      Q2<pdf->Q2Min() || Q2>pdf->Q2Max()) return 0.;
  if (Q2<sqr(2.*flav.Mass(true)))         return 0.;
  pdf->Calculate(x,Q2);
  return pdf->GetXPDF(flav.Bar());
}

Weight * Kernel::MakeWeight(const Splitting & split,const double & overfac) {
  double weight   = (p_gauge->Charge() * (*p_gauge)(split) *
		     (*p_sf)(split) * p_sf->Jacobean(split));
  double realover = (p_gauge->Charge() * p_gauge->OverEstimate(split) *
		     p_sf->OverEstimate(split));
  double over     = ( (dabs(weight)<realover) ?
		      (weight>0.?1.:-1.)*realover : overfac*weight );
  return new Weight(weight,over,realover);
}

bool Kernel::FillOffsprings(Splitting & split) {
  int beam = split.GetSplitter()->Beam();
  for (size_t i=0;i<m_flavs.size();i++) {
    Parton * parton = new Parton(m_flavs[i],split.Mom(i));
    if (i==1) parton->SetBeam(beam);
    parton->SetColor(split.Col(i));
    split.SetParton(i,parton);
  }
  return true;
}


ostream & CFPSHOWER::operator<<(ostream &s,const Kernel & kernel) {
  Flavour_Vector flavs = kernel.GetFlavs(); 
  s<<"Kernel["<<(kernel.GetSF()->Name())<<" / "<<(kernel.GetGauge()->Name())<<"]: "
   <<kernel.GetSplit()<<" --> ";
  for (size_t i=0;i<flavs.size();i++) s<<flavs[i]<<" ";
  s<<" tag sequence = { ";
  for (size_t i=0;i<kernel.GetSF()->TagSequence().size();i++)
    s<<kernel.GetSF()->TagSequence()[i]<<" ";
  s<<"}\n";
  return s;
}

DECLARE_GETTER(Kernel,"Kernel",Kernel,Kernel_Info);

Kernel * Getter<Kernel,Kernel_Info,Kernel>::operator()(const Parameter_Type & info) const
{
  SF_Base    * sf = SF_Getter::GetObject(info.SFName(),info);
  Gauge_Base * gp = GP_Getter::GetObject(info.GPName(),info);
  if (sf) {
    msg_Out()<<"   * try to init new Kernel("<<info.GetSplit()<<" -> ";
    for (size_t i=0;i<info.GetFlavs().size();i++) msg_Out()<<info.GetFlavs()[i]<<" ";
    msg_Out()<<" ["<<info.SFName()<<"]\n";
  }
  if (!sf || !gp) {
    if (sf) delete sf;
    if (gp) delete gp;
    msg_Out()<<"   * failed.\n";
    return NULL;
  }
  Kernel * kernel = new Kernel(info);
  kernel->SetSF(sf);
  kernel->SetGauge(gp);
  msg_Out()<<"Found "<<kernel->GetSplit()<<" --> ";
  for (size_t i=0;i<kernel->GetFlavs().size();i++) msg_Out()<<kernel->GetFlavs()[i]<<" ";
  msg_Out()<<" ["<<info.SFName()<<"]\n";
  return kernel;
}

void Getter<Kernel,Kernel_Info,Kernel>::PrintInfo(ostream &str,const size_t width) const
{
  str<<"Splitting Kernel";
}
