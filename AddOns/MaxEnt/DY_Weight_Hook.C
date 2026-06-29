#include "SHERPA/Tools/Userhook_Base.H"
#include "ATOOLS/Org/Message.H"
#include "SHERPA/Single_Events/Event_Handler.H"

#include "ATOOLS/Phys/Weight_Info.H"
#include <algorithm>
#include "MODEL/Main/Running_AlphaS.H"
#include "SHERPA/Main/Sherpa.H"
#include "SHERPA/Initialization/Initialization_Handler.H"
#include "CSSHOWER++/Main/CS_Shower.H"
#include "PDF/Main/Shower_Base.H"
#include "PHASIC++/Process/MCatNLO_Process.H"

#include <nlohmann/json.hpp>

using json = nlohmann::json;
using namespace ATOOLS;
using namespace SHERPA;

class DY_Weight_Hook : public Userhook_Base, public Tag_Replacer {
private:

  Sherpa* p_sherpa;
  double m_rt, m_lnrt, m_dphi, m_lndphi, m_gating[2];
  std::vector<double> m_lss;
  std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;
  std::vector<std::string> m_names;

public:

  DY_Weight_Hook(const Userhook_Arguments args) :
    Userhook_Base("DY_Weight"),
    p_sherpa(args.p_sherpa)
  {
    DEBUG_FUNC("");
    Settings& s = Settings::GetMainSettings();
    std::string fname = s["DY_WEIGHT_FILE"].
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
    m_lndphi=log(m_dphi=l1.DPhi(l2));
    double beta(Beta((l1+l2).PPerp())), wnom;
    msg_Debugging()<<"q_T = "<<(l1+l2).PPerp()<<", r_T = "
		   <<m_rt<<", \\Delta\\phi = "<<m_dphi<<"\n";
    auto me_w_info = (*blobs->FindFirst(btp::Signal_Process))
      ["MEWeightInfo"]->Get<ME_Weight_Info*>();
    Weights_Map &wmap = (*blobs->FindFirst(btp::Signal_Process))
      ["WeightsMap"]->Get<Weights_Map>();
    for (size_t i(0);i<m_calcs.size();++i)
      if (m_names[i]=="central") {
	wnom=m_calcs[i]->Calculate()->Get<double>();
	msg_Debugging()<<"nominal: ln(w) = "<<wnom<<", shift = "<<m_lss[i]<<"\n";
	wnom=beta*exp(wnom-m_lss[i])+(1.-beta);
	msg_Debugging()<<"nominal: w = "<<wnom<<" (\\beta = "<<beta<<")\n";
	wmap*=wnom;
	*me_w_info*=wnom;
      }
    for (size_t i(0);i<m_calcs.size();++i)
      if (m_names[i]!="central") {
	double w=m_calcs[i]->Calculate()->Get<double>();
	msg_Debugging()<<m_names[i]<<": ln(w) = "<<w<<", shift = "<<m_lss[i]<<"\n";
	w=beta*exp(w-m_lss[i])+(1.-beta);
	msg_Debugging()<<m_names[i]<<": w = "<<w<<" (\\beta = "<<beta<<")\n";
	wmap["MaxEnt"][m_names[i]]=w/wnom;
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
