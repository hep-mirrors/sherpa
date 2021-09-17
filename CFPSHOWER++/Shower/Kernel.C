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
  p_msel(NULL), m_enhance(1.)
{
  for (size_t beam=0;beam<2;beam++) p_pdf[beam] = NULL;
}

double Kernel::Integral(Splitting & split,Configuration & config,
			const Mass_Selector * msel) {
  split.Init(this,msel);
  if (p_kin->Init(split,config,msel)) {
    double SF = p_sf->Integral(split), GF = p_gauge->OverEstimate(split); 
    double I = (p_gauge->Charge() * SF * GF);
    //msg_Out()<<"   * "<<GetSplit()<<" -->";
    //for (size_t i=0;i<GetFlavs().size();i++)
    //  msg_Out()<<" "<<GetFlavs()[i]<<" ["<<Tags()[i]<<"]";
    //msg_Out()<<" ("<<p_kin->Scheme()<<"): "
    //	     <<"SF(["<<split.Zmin()<<", "<<split.Zmax()<<"] = "
    //	     <<SF<<") * "
    //	     <<"gauge( = "<<p_gauge->Charge()<<" * "
    //	     <<GF<<" = "<<(p_gauge->Charge()*GF)<<") --> "
    //	     <<"I/(2pi) = "<<(I/(2.*M_PI))<<" for "
    //	     <<"Q2/t0 = "<<split.Q2()<<"/"<<split.Tcut()<<"\n"
    //	     <<"     (enhance = "<<m_enhance<<", alphaS max = "
    //	     <<p_gauge->AlphaSMax(split)<<").\n";
    return m_enhance * I/(2.*M_PI);
  }
  return 0.;
}

bool Kernel::Generate(Splitting & split,Configuration & config,
		      const Mass_Selector * msel,const double & overfac) {
  // Sequence of generation of one splitting after kernel is selected:
  // - attach pointer of kernel to Splitting, init invariants for splitting:
  //   setting Q^2, masses, and Q^2_red
  //   (make sure splitting is allowed when splitting is massive).
  // - generate a splitting kinematics, essentially the z and phi
  // - calculate kinematic invariants
  // - construct outgoing four-momenta of splitting products and spectator
  // - if this is successful, calculate weights, apply hit-or-miss, attach
  //   acceptance and rejection weights to the splitter-spectator pair
  split.Init(this,msel);
  if (p_kin->Init(split,config,msel)) {
    p_sf->GeneratePoint(split);
    if (!(*p_kin)(split,config)) return false; 
    //msg_Out()<<" --> kt2 = "<<split.T()<<", z = "<<split.Z()<<"\n"; 
    //msg_Out()<<METHOD<<"("<<this<<", "<<p_sf<<") "
    //     <<"with the following momenta: \n"
    //	     <<"   pi = "<<split.Mom(0)<<"\n"
    //	     <<"   pj = "<<split.Mom(1)<<"\n"
    //	     <<"   pk = "<<split.Mom(2)<<"\n";
    p_kin->CalculateJacobean(split,config);
    split.SetWeight(MakeWeight(split,overfac));
    if ((*split.GetWeight())()>=ran->Get()) {
      split.GetSplitter()->AddWeight(split,true);
      if (!p_gauge->SetColours(split)) {
	/*msg_Out()<<" ---> colour setting failed.\n";*/
	return false;
      }
      return true;
    }
    split.GetSplitter()->AddWeight(split,false);
  }
  return false;
}

bool Kernel::UpdateSystem(Splitting & split,Configuration & config) {
  return p_kin->UpdateSystem(split,config);
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
  //msg_Out()<<"------------------------------------------------------\n";
  double weight   = (p_gauge->Charge() * (*p_gauge)(split) *
		     (*p_sf)(split) * p_kin->Weight());
  //msg_Out()<<METHOD<<" for colours: "<<p_gauge->Charge()<<" * "
  //	   <<"as("<<split.T()<<") = "
  //	   <<(*p_gauge)(split)<<" = "<<(p_gauge->Charge()*(*p_gauge)(split))
  //	   <<"/"<<p_gauge->OverEstimate(split)<<"\n";
  double realover = (p_gauge->Charge() * p_gauge->OverEstimate(split) *
		     p_sf->OverEstimate(split));
  double over     = overfac * weight;
  if (dabs(weight)<realover) {
    if (weight>=0.) over = realover;
    else            over = -realover;
  }
  //msg_Out()<<METHOD<<":\n"
  //	   <<"   [kin = "<<p_kin->Weight()<<"] * [alpha = "
  //	   <<(*p_gauge)(split)<<" / "<<p_gauge->OverEstimate(split)<<" = "
  //	   <<((*p_gauge)(split)/p_gauge->OverEstimate(split))<<"]\n"
  //	   <<"   [SF = "<<(p_gauge->Charge() * (*p_sf)(split))<<"] / "
  //	   <<(p_gauge->Charge() * p_sf->OverEstimate(split))
  //	   <<" for kt2 = "<<split.T()<<")\n"
  //	   <<"   ==> weight = "<<(weight/over)<<" "
  //	   <<"(realover = "<<realover<<").\n";
  return new Weight(weight,over,m_enhance*realover);
}

bool Kernel::FillOffsprings(Splitting & split) {
  int beam = split.GetSplitter()->Beam();
  split.Clear();
  //msg_Out()<<METHOD<<" for "<<GetFlavs().size()<<" outgoing partons.\n";
  for (size_t i=0;i<GetFlavs().size();i++) {
    Parton * parton = new Parton(GetFlavs()[i],split.Mom(i));
    if (i==0) parton->SetBeam(beam);
    parton->SetColor(p_gauge->GetColor(i));
    split.AddParton(parton);
    //msg_Out()<<" * "<<(*parton)<<"\n";
  }
  return true;
}


ostream & CFPSHOWER::operator<<(ostream &s,const Kernel & kernel) {
  Flavour_Vector flavs = kernel.GetFlavs(); 
  s<<"Kernel["<<kernel.GetSF()->Name()<<" / "<<kernel.GetGauge()->Name()<<"]: "
   <<kernel.GetSplit()<<" --> ";
  for (size_t i=0;i<flavs.size();i++) s<<flavs[i]<<" ";
  s<<" tag sequence = { ";
  for (size_t i=0;i<kernel.GetSF()->Tags().size();i++)
    s<<kernel.GetSF()->Tags()[i]<<" ";
  s<<"}\n";
  return s;
}

DECLARE_GETTER(Kernel,"Kernel",Kernel,Kernel_Info);

Kernel * Getter<Kernel,Kernel_Info,Kernel>::
operator()(const Parameter_Type & info) const
{
  //msg_Out()<<"\n\n***** looking for "<<info
  //	   <<"      "
  //	   <<info.SFName()<<", "<<info.GPName()<<", "<<info.KinName()<<"\n";
  SF_Base         * sf    = SF_Getter::GetObject(info.SFName(),info);
  Gauge_Base      * gauge = GP_Getter::GetObject(info.GPName(),info);
  Kinematics_Base * kin   = Kinematics_Getter::GetObject(info.KinName(),info);
  if (!sf || !gauge || !kin) {
    //msg_Out()<<"--> failed: "
    //	     <<"SF: "<<sf<<", Gauge: "<<gauge<<", Kin: "<<kin<<"\n";
    if (sf)    delete sf;
    if (gauge) delete gauge;
    if (kin)   delete kin;
    return NULL;
  }
  Kernel * kernel = new Kernel(info);
  kernel->SetSF(sf);
  kernel->SetGauge(gauge);
  kernel->SetKinematics(kin);
  //if (msg_LevelIsDebugging()) {
    msg_Out()<<"***** Found "<<kernel->GetSplit()<<" --> ";
    for (size_t i=0;i<kernel->GetFlavs().size();i++)
      msg_Out()<<kernel->GetFlavs()[i]<<" ("<<kernel->GetSF()->Tags()[i]<<") ";
    msg_Out()<<" ["
	     <<kernel->GetSF()->Name()<<" + "
	     <<kernel->GetGauge()->Name()<<" + "
	     <<kernel->GetKinematics()->Name()<<"]\n";
    //}
  return kernel;
}

void Getter<Kernel,Kernel_Info,Kernel>::PrintInfo(ostream &str,const size_t width) const
{
  str<<"Parton Shower Kernel\n";
}
