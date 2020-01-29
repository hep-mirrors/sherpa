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
  m_type(info.Type()), p_msel(NULL), m_enhance(1.)
{
  for (size_t beam=0;beam<2;beam++) p_pdf[beam] = NULL;
}

double Kernel::Integral(Splitting & split,const Mass_Selector * msel) {
  split.SetKernel(this);
  split.InitSplitting(msel);
  double I = (p_gauge->Charge() * p_sf->Integral(split) *
	      p_gauge->OverEstimate(split));
  // msg_Out()<<"   *** "<<METHOD<<"("<<GetSplit()<<" -->";
  // for (size_t i=0;i<GetFlavs().size();i++)
  //   msg_Out()<<" "<<GetFlavs()[Tags()[i]]<<" ["<<Tags()[i]<<"]";
  // msg_Out()<<"): SF( = "<<p_sf->Integral(split)<<") * gauge( = "
  // 	   <<(p_gauge->Charge() * p_gauge->OverEstimate(split))<<") = "<<I<<" for "
  // 	   <<"Q2/t0 = "<<split.Q2()<<"/"<<split.tcut()<<" (enhance = "<<m_enhance<<").\n";
  return m_enhance * I/(2.*M_PI);
}

bool Kernel::Generate(Splitting & split,const Mass_Selector * msel,
		      const double & overfac) {
  // Sequence of generation of one splitting after kernel is selected:
  // - attach pointer of kernel to Splitting, init invairants for splitting:
  //   setting Q^2, masses, and Q^2_red
  //   (make sure splitting is allowed when splitting is massive).
  // - generate a splitting kinematics, essentially the z and phi
  // - calculate kinematic invariants
  // - construct outgoing four-momenta of splitting products and spectator
  // - if this is successful, calculate weights, apply hit-or-miss, attach
  //   acceptance and rejection weights to the splitter-spectator pair
  split.SetKernel(this);
  split.InitSplitting(msel);
  p_sf->GeneratePoint(split);
  if (!p_gauge->SetColours(split)) {
    //msg_Out()<<"   *** "<<METHOD<<" couldn't set colours.\n";
    return false;
  }
  if (p_sf->Construct(split)) {
    if (split.GetWeight()) { delete split.GetWeight(); split.SetWeight(NULL); } 
    split.SetWeight(MakeWeight(split,overfac));
    if ((*split.GetWeight())()>=ran->Get()) {
      split.GetSplitter()->AddWeight(split,true);
      // msg_Out()<<"   *** "<<METHOD
      // 	       <<" success with weight = "<<(*split.GetWeight())()<<", too small.\n";
      return true;
    }
    else {
      // msg_Out()<<"   *** "<<METHOD
      // 	       <<" rejected with weight = "<<(*split.GetWeight())()<<", too small.\n";
      split.GetSplitter()->AddWeight(split,false);
      return false;
    }
  }
  //msg_Out()<<"   *** "<<METHOD<<" couldn't construct splitting.\n";
  return false;
}

bool Kernel::UpdateKinematics(Splitting & split) {
  return p_sf->Construct(split,1);
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
  double SF       = (*p_sf)(split);
  double weight   = (p_gauge->Charge() * (*p_gauge)(split) *
		     SF * p_sf->Jacobean(split));
  double realover = (p_gauge->Charge() * p_gauge->OverEstimate(split) *
		     p_sf->OverEstimate(split));
  double over     = overfac * weight;
  if (dabs(weight)<realover) {
    if (weight>=0.) over = realover;
    else            over = -realover;
  }
  // msg_Out()<<METHOD<<": "<<(p_gauge->Charge() * (*p_gauge)(split))<<"(Col) * "
  // 	   <<SF<<"(SF, z = "<<split.z()<<") * "
  // 	   <<p_sf->Jacobean(split)<<"(J) / "
  // 	   <<"over = "<<p_sf->OverEstimate(split)
  // 	   <<" -> "<<weight<<"/"<<over<<".\n";
  return new Weight(weight,over,m_enhance*realover);
}

bool Kernel::FillOffsprings(Splitting & split) {
  int beam = split.GetSplitter()->Beam();
  for (size_t i=0;i<GetFlavs().size();i++) {
    Parton * parton = new Parton(GetFlavs()[i],p_sf->GetMoms()[i]);
    if (i==0) parton->SetBeam(beam);
    parton->SetColor(p_gauge->GetColor(i));
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
  for (size_t i=0;i<kernel.GetSF()->Tags().size();i++) s<<kernel.GetSF()->Tags()[i]<<" ";
  s<<"}\n";
  return s;
}

DECLARE_GETTER(Kernel,"Kernel",Kernel,Kernel_Info);

Kernel * Getter<Kernel,Kernel_Info,Kernel>::operator()(const Parameter_Type & info) const
{
  SF_Base    * sf = SF_Getter::GetObject(info.SFName(),info);
  Gauge_Base * gp = GP_Getter::GetObject(info.GPName(),info);
  if (!sf || !gp) {
    if (sf) delete sf;
    if (gp) delete gp;
    return NULL;
  }
  Kernel * kernel = new Kernel(info);
  kernel->SetSF(sf);
  kernel->SetGauge(gp);
  if (msg_LevelIsDebugging()) {
    msg_Out()<<"***** Found "<<kernel->GetSplit()<<" --> ";
    for (size_t i=0;i<kernel->GetFlavs().size();i++)
      msg_Out()<<kernel->GetFlavs()[i]<<" ";
    msg_Out()<<" ["<<info.SFName()<<"]\n";
  }
  return kernel;
}

void Getter<Kernel,Kernel_Info,Kernel>::PrintInfo(ostream &str,const size_t width) const
{
  str<<"Splitting Kernel";
}
