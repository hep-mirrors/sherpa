#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Single_Events/Event_Handler.H"

#include "ATOOLS/Phys/Weight_Info.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "PHASIC++/Process/MCatNLO_Process.H"
#include "PHASIC++/Process/Process_Base.H"

#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace ATOOLS;
using namespace PHASIC;
using namespace SHERPA;

class DY_Weight_Hook : public Userhook_Base, public Tag_Replacer {
private:

  Sherpa* p_sherpa;
  int m_jetmode;
  double m_rt, m_lnrt, m_dphi, m_lndphi, m_gating[2];
  std::vector<double> m_lss;
  std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;
  std::vector<std::string> m_names, m_vtags;

public:

  DY_Weight_Hook(const Userhook_Arguments args) :
    Userhook_Base("DY_Weight"),
    p_sherpa(args.p_sherpa)
  {
    DEBUG_FUNC("");
    Settings& s = Settings::GetMainSettings();
    m_jetmode=s["DY_JET_MODE"].SetDefault(0).Get<int>();
    std::string fname=s["DY_WEIGHT_FILE"].
      SetDefault("lambda_export_variations.json").Get<std::string>();
    msg_Debugging()<<"DY_Weight user hook reading from '"
		   <<fname<<"'."<<std::endl;
    std::ifstream f(fname);
    json data=json::parse(f);
    auto tags=data["moments"];
    msg_Debugging()<<"tags = "<<tags<<"\n";
    auto gating=data["gating"]["window_GeV"];
    msg_Debugging()<<"gating = "<<gating<<"\n";
    m_gating[0]=gating[0];
    m_gating[1]=gating[1];
    m_names=data["scheme_names"];
    for (size_t i(0);i<m_names.size();++i) {
      std::string var(m_names[i]);
      auto vals = data["schemes"][var]["lambda_physical"];
      m_lss.push_back(data["schemes"][var]["log_norm_shift"]);
      msg_Debugging()<<"vals[\""<<var<<"\"] = "<<vals<<"\n";
      msg_Debugging()<<"lss[\""<<var<<"\"] = "<<m_lss.back()<<"\n";
      m_calcs.push_back(new Algebra_Interpreter());
      m_calcs.back()->SetTagReplacer(this);
      m_calcs.back()->AddTag("const","1.0");
      m_calcs.back()->AddTag("rt","1.0");
      m_calcs.back()->AddTag("lnrt","1.0");
      m_calcs.back()->AddTag("dphi","1.0");
      m_calcs.back()->AddTag("lndphi","1.0");
      std::string expr, repl("×");
      for (size_t pos, j(0);j<vals.size();++j) {
	std::string tag(tags[j]);
	while ((pos=tag.find(repl))!=std::string::npos)
	  tag.replace(pos,repl.length(),"*");
	expr+=(j&&vals[j]>0?"+":"")+ToString(vals[j])
	  +"*"+std::string(tag);
      }
      msg_Debugging()<<"expr = "<<expr<<"\n";
      m_calcs.back()->Interprete(expr);
      if (msg_LevelIsIODebugging()) m_calcs.back()->PrintEquation();
      m_vtags.push_back("Nominal");
      if (m_names[i]=="0p5MuR") m_vtags.back()="MUR=0.5__MUF=1__";
      if (m_names[i]=="0p5MuF") m_vtags.back()="MUR=1__MUF=0.5__";
      if (m_names[i]=="0p5MuRF") m_vtags.back()="MUR=0.5__MUF=0.5__";
      if (m_names[i]=="2MuR") m_vtags.back()="MUR=2__MUF=1__";
      if (m_names[i]=="2MuF") m_vtags.back()="MUR=1__MUF=2__";
      if (m_names[i]=="2MuRF") m_vtags.back()="MUR=2__MUF=2__";
    }
  }

  ~DY_Weight_Hook()
  {
    for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
  }

  std::string ReplaceTags(std::string &expr) const
  {
    return m_calcs.front()->ReplaceTags(expr);
  }

  Term *ReplaceTags(Term *term) const
  {
    switch (term->Id()) {
    case 1: term->Set(1.); return term;
    case 2: term->Set(m_rt); return term;
    case 3: term->Set(m_lnrt); return term;
    case 4: term->Set(m_dphi); return term;
    case 5: term->Set(m_lndphi); return term;
    }
    return term;
  }

  void AssignId(Term *term)
  {
    if (term->Tag()=="const") term->SetId(1);
    if (term->Tag()=="rt") term->SetId(2);
    if (term->Tag()=="lnrt") term->SetId(3);
    if (term->Tag()=="dphi") term->SetId(4);
    if (term->Tag()=="lndphi") term->SetId(5);
  }

  double Beta(const double &qt)
  {
    double t((qt-m_gating[0])/(m_gating[1]-m_gating[0]));
    if (t<0.) t=0.;
    if (t>1.) t=1.;
    return 1.-(6.*pow(t,5)-15.*pow(t,4)+10*pow(t,3));
  }

  ATOOLS::Return_Value::code Run(ATOOLS::Blob_List* blobs)
  {
    DEBUG_FUNC(p_sherpa->GetInitHandler()->
	       GetMatrixElementHandler()->
	       Process()->Parent()->Name());
    Vec4D l1, l2;
    Blob *psb(blobs->FindFirst(btp::Shower));
    for (size_t i(0);i<psb->NOutP();++i)
      if (psb->OutParticle(i)->Flav().IsLepton()) {
	if (l1==Vec4D()) l1=psb->OutParticle(i)->Momentum();
	else l2=psb->OutParticle(i)->Momentum();
      }
    m_lnrt=log(m_rt=(l1+l2).PPerp()/(l1+l2).Mass());
    m_lndphi=log(m_dphi=M_PI-l1.DPhi(l2));
    msg_Debugging()<<"q_T = "<<(l1+l2).PPerp()<<", r_T = "
		   <<m_rt<<", \\Delta\\phi = "<<m_dphi<<"\n";
    Weights_Map &wmap = (*blobs->FindFirst(btp::Signal_Process))
      ["WeightsMap"]->Get<Weights_Map>();
    double wnom(1.0), wew(1.0);
    Weights_Map::const_iterator wit(wmap.find("ASSOCIATED_CONTRIBUTIONS"));
    if (wit!=wmap.end()) {
      for (size_t l(0);l<wit->second.Size();++l)
	if (wit->second.Name(l)=="EW") {
	  wew=wit->second[l]/wit->second.Nominal();
	  break;
	}
    }
    msg_Debugging()<<"EW weight: "<<wew<<"\n";
    wit=wmap.find("All");
    if (wit==wmap.end()) THROW(fatal_error,"Variation set 'All' not found");
    for (size_t k(0);k<wit->second.Size();++k)
      if (wit->second.Name(k).find("Nominal")==0) {
	wnom=wit->second[k];
	break;
      }
    msg_Debugging()<<"nominal weight: "<<wnom<<"\n";
    for (size_t i(0);i<m_calcs.size();++i) {
      double w=m_calcs[i]->Calculate()->Get<double>();
      msg_Debugging()<<m_names[i]<<": ln(w) = "<<w<<", shift = "<<m_lss[i]<<"\n";
      double svweight(1.0);
      std::string svname("None");
      if (m_vtags[i]!="") {
	for (size_t k(0);k<wit->second.Size();++k)
	  if (wit->second.Name(k).find(m_vtags[i])==0) {
	    svweight=wit->second[k]/wnom;
	    svname=wit->second.Name(k);
	    break;
	  }
      }
      if (m_jetmode==1) {
	double beta(Beta((l1+l2).PPerp()));
	w=beta*exp(w-m_lss[i])+(1.-beta)*svweight;
	msg_Debugging()<<m_names[i]<<": w = "<<w<<" (\\beta = "<<beta
		       <<") <-> "<<svweight<<" ("<<svname<<")\n";
      }
      else {
	Process_Base *proc=p_sherpa->GetInitHandler()->
	  GetMatrixElementHandler()->Process()->Parent();
	size_t nout=proc->NOut();
	if (proc->Get<MCatNLO_Process>()!=nullptr) --nout;
	if (nout>3) w=svweight;
	else w=exp(w-m_lss[i]);
	msg_Debugging()<<m_names[i]<<": w = "<<w<<" (n_{jet} = "<<nout-2
		       <<") <-> "<<svweight<<" ("<<svname<<")\n";
      }
      wmap["MaxEnt_QCD"][m_names[i]]=w;
      wmap["MaxEnt_EW"][m_names[i]]=w*wew;
    }
    return Return_Value::Nothing;
  }

  void Finish() {}

};

DECLARE_GETTER(DY_Weight_Hook,"DY_MaxEnt_Weight",
               Userhook_Base,Userhook_Arguments);

Userhook_Base *ATOOLS::Getter<Userhook_Base,Userhook_Arguments,DY_Weight_Hook>::
operator()(const Userhook_Arguments &args) const
{
  return new DY_Weight_Hook(args);
}

void ATOOLS::Getter<Userhook_Base,Userhook_Arguments,DY_Weight_Hook>::
PrintInfo(std::ostream &str,const size_t width) const
{
  str<<"DY maxent weight userhook";
}
