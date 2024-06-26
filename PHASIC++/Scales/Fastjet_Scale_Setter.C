#include "ATOOLS/Org/CXXFLAGS_PACKAGES.H"

#include "PHASIC++/Scales/Scale_Setter_Base.H"

#include "PHASIC++/Scales/Tag_Setter.H"
#include "PHASIC++/Process/Process_Base.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Model_Base.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Fastjet_Helpers.H"

namespace PHASIC {

  class Fastjet_Scale_Setter: public Scale_Setter_Base {
  protected:

    std::vector<ATOOLS::Algebra_Interpreter*> m_calcs;

    Tag_Setter m_tagset;

    fjcore::JetDefinition *p_jdef;

    double m_ptmin, m_etmin, m_eta, m_y;

    ATOOLS::Flavour_Vector m_f;

    int m_mode, m_bmode;

    double ASMeanScale(const std::vector<double> &mu,
		       const size_t &offset) const;

  public:

    Fastjet_Scale_Setter(const Scale_Setter_Arguments &args);

    ~Fastjet_Scale_Setter();

    double Calculate(const std::vector<ATOOLS::Vec4D> &p,
		     const size_t &mode);

    void SetScale(const std::string &mu2tag,
		  ATOOLS::Algebra_Interpreter &mu2calc);

    const ATOOLS::Vec4D_Vector &Momenta() const;

  };// end of class Fastjet_Scale_Setter

}// end of namespace PHASIC

using namespace PHASIC;
using namespace ATOOLS;

DECLARE_GETTER(Fastjet_Scale_Setter,"FASTJET",
	       Scale_Setter_Base,Scale_Setter_Arguments);

Scale_Setter_Base *ATOOLS::Getter
<Scale_Setter_Base,Scale_Setter_Arguments,Fastjet_Scale_Setter>::
operator()(const Scale_Setter_Arguments &args) const
{
  return new Fastjet_Scale_Setter(args);
}

void ATOOLS::Getter<Scale_Setter_Base,Scale_Setter_Arguments,
		    Fastjet_Scale_Setter>::
PrintInfo(std::ostream &str,const size_t width) const
{ 
  str<<"variable scale scheme using fast jets";
}

Fastjet_Scale_Setter::Fastjet_Scale_Setter
(const Scale_Setter_Arguments &args):
  Scale_Setter_Base(args), m_tagset(this),
  p_jdef(NULL)
{
  rpa->gen.AddCitation(1,"FastJet is published under \\cite{Cacciari:2011ma}.");
  std::string jtag(args.m_scale);
  size_t pos(jtag.find("FASTJET["));
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
  jtag=jtag.substr(pos+8);
  pos=jtag.find(']');
  if (pos==std::string::npos)
    THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
  jtag=jtag.substr(0,pos);
  Data_Reader read(" ",",","#","=");
  read.AddIgnore(":");
  read.SetString(jtag);
  m_mode=read.StringValue<int>("M",1);
  m_bmode=read.StringValue<int>("B",0);
  m_ptmin=read.StringValue<double>("PT",0.0);
  m_etmin=read.StringValue<double>("ET",0.0);
  m_eta=read.StringValue<double>("Eta",100.0);
  m_y=read.StringValue<double>("Y",100.0);
  double R(read.StringValue<double>("R",0.4));
  double f(read.StringValue<double>("f",0.75));
  std::string algo(read.StringValue<std::string>("A","antikt"));
  fjcore::JetAlgorithm ja(fjcore::kt_algorithm);
  if (algo=="cambridge") ja=fjcore::cambridge_algorithm;
  if (algo=="antikt") ja=fjcore::antikt_algorithm;
  std::string reco(read.StringValue<std::string>("C","E"));
  fjcore::RecombinationScheme recom(fjcore::E_scheme);
  if (reco=="pt") recom=fjcore::pt_scheme;
  if (reco=="pt2") recom=fjcore::pt2_scheme;
  if (reco=="Et") recom=fjcore::Et_scheme;
  if (reco=="Et2") recom=fjcore::Et2_scheme;
  if (reco=="BIpt") recom=fjcore::BIpt_scheme;
  if (reco=="BIpt2") recom=fjcore::BIpt2_scheme;
  p_jdef=new fjcore::JetDefinition(ja,R,recom);
  m_f=p_proc->Flavours();
  m_p.resize(p_proc->NIn()+p_proc->NOut());
  std::string tag(args.m_scale);
  std::vector<std::string> ctags;
  while (true) {
    size_t pos(tag.find('{'));
    if (pos==std::string::npos) {
      if (!m_calcs.empty()) break;
      else { THROW(fatal_error,"Invalid scale '"+args.m_scale+"'"); }
    }
    tag=tag.substr(pos+1);
    pos=tag.find('}');
    if (pos==std::string::npos) 
      THROW(fatal_error,"Invalid scale '"+args.m_scale+"'");
    std::string ctag(tag.substr(0,pos));
    tag=tag.substr(pos+1);
    m_calcs.push_back(new Algebra_Interpreter());
    m_calcs.back()->AddFunction(MODEL::as->GetAIGMeanFunction());
    m_calcs.back()->SetTagReplacer(&m_tagset);
    if (m_calcs.size()==1) m_tagset.SetCalculator(m_calcs.back());
    ctags.push_back(ctag);
  }
  for (size_t i(p_proc->NIn());i<m_f.size();++i)
    if (!ToBeClustered(m_f[i], m_bmode)) m_scale.push_back(0.0);
  m_scale.resize(Max(m_scale.size(),m_calcs.size()+stp::size));
  for (size_t i(0);i<m_calcs.size();++i)
    SetScale(ctags[i],*m_calcs[i]);
  m_scale.resize(stp::size);
  SetCouplings();
}

Fastjet_Scale_Setter::~Fastjet_Scale_Setter()
{
  for (size_t i(0);i<m_calcs.size();++i) delete m_calcs[i];
  delete p_jdef;
}

const Vec4D_Vector &Fastjet_Scale_Setter::Momenta() const
{
  return m_p;
}

double Fastjet_Scale_Setter::Calculate
(const std::vector<ATOOLS::Vec4D> &_momenta,const size_t &mode) 
{
  std::vector<ATOOLS::Vec4D> momenta(_momenta);
  if (p_proc->Flavours()[0].IsLepton()&&rpa->gen.Beam2().IsHadron()) {
    msg_Debugging()<<METHOD<<"(): Boost to Breit frame {\n";
    Vec4D pp(rpa->gen.PBeam(1)), qq(momenta[0]-momenta[2]);
    Poincare cms(pp+qq);
    double Q2(-qq.Abs2()), x(Q2/(2.0*pp*qq)), E(sqrt(Q2)/(2.0*x));
    Vec4D p(sqrt(E*E+pp.Abs2()),0.0,0.0,-E);
    Vec4D q(0.0,0.0,0.0,2.0*x*E);
    cms.Boost(pp);
    cms.Boost(qq);
    Poincare zrot(pp,-Vec4D::ZVEC);
    zrot.Rotate(pp);
    zrot.Rotate(qq);
    Poincare breit(p+q);
    breit.BoostBack(pp);
    breit.BoostBack(qq);
    if (!IsEqual(pp,p,1.0e-3) || !IsEqual(qq,q,1.0e-3))
      msg_Error()<<METHOD<<"(): Boost error."<<std::endl;
    for (int i(0);i<momenta.size();++i) {
      msg_Debugging()<<"  "<<i<<": "<<momenta[i];
      cms.Boost(momenta[i]);
      zrot.Rotate(momenta[i]);
      breit.BoostBack(momenta[i]);
      msg_Debugging()<<" -> "<<momenta[i]<<"\n";
    }
    msg_Debugging()<<"}\n";
  }
  m_p.resize(2);
  m_p[0]=-momenta[0];
  m_p[1]=-momenta[1];
  std::vector<fjcore::PseudoJet> input;
  for (size_t i(p_proc->NIn());i<momenta.size();++i)
    if (!ToBeClustered(m_f[i], m_bmode)) m_p.push_back(momenta[i]);
    else {
      input.push_back(MakePseudoJet(m_f[i], momenta[i]));
    }
  fjcore::ClusterSequence cs(input,*p_jdef);
  std::vector<fjcore::PseudoJet>
    jets(fjcore::sorted_by_pt(cs.inclusive_jets()));
  for (size_t i(0);i<jets.size();++i) {
    if (m_bmode==0 || BTag(jets[i], m_bmode)) {
      Vec4D pj(jets[i].E(),jets[i].px(),jets[i].py(),jets[i].pz());
      if (pj.PPerp()>m_ptmin && pj.EPerp()>m_etmin &&
          (m_eta==100.0 || dabs(pj.Eta())<m_eta) &&
          (m_y==100.0 || dabs(pj.Y())<m_y)) m_p.push_back(pj);
    }
  }
  m_scale.resize(Max(m_scale.size(),m_calcs.size()+stp::size));
  for (size_t idx(stp::size), i(0);i<input.size();++i)
    m_scale[idx++]=cs.exclusive_dmerge_max(i);
  for (size_t i(0);i<m_calcs.size();++i)
    m_scale[m_mode==1?i+stp::size:i]=m_calcs[i]->Calculate()->Get<double>();
  if (m_mode==1) m_scale[stp::res]=m_scale[stp::fac]=
		   m_scale[stp::ren]=ASMeanScale(m_scale,stp::size);
  else for (size_t i(m_calcs.size());i<stp::size;++i) m_scale[i]=m_scale[0];
  msg_Debugging()<<METHOD<<"(): Set {\n"
		 <<"  Q     = "<<sqrt(m_scale[stp::res])<<"\n"
		 <<"  \\mu_f = "<<sqrt(m_scale[stp::fac])<<"\n"
		 <<"  \\mu_r = "<<sqrt(m_scale[stp::ren])<<"\n";
  for (size_t i(stp::size);i<m_scale.size();++i)
    msg_Debugging()<<"  \\mu_"<<i<<" = "<<sqrt(m_scale[i])<<"\n";
  msg_Debugging()<<"} <- "<<p_proc->Name()<<"\n";
  m_scale.resize(stp::size);
  return m_scale[stp::fac];
}

double Fastjet_Scale_Setter::ASMeanScale
(const std::vector<double> &mu,const size_t &offset) const
{
  msg_Debugging()<<"Setting scales {\n";
  double mur2(1.0), as(1.0), oqcd(0.0);
  for (size_t i(offset);i<offset+m_calcs.size();++i) {
    double cas(MODEL::as->BoundedAlphaS(mu[i]));
    msg_Debugging()<<"  \\mu_{"<<i<<"} = "
		   <<sqrt(mu[i])<<", as = "<<cas<<"\n";
    mur2*=mu[i];
    as*=cas;
    ++oqcd;
  }
  if (oqcd==0.0) THROW(fatal_error,"No jets!");
  as=pow(as,1.0/oqcd);
  mur2=MODEL::as->WDBSolve(as,MODEL::as->CutQ2(),sqr(rpa->gen.Ecms()));
  if (!IsEqual((*MODEL::as)(mur2),as))
    msg_Error()<<METHOD<<"(): Failed to determine \\mu."<<std::endl; 
  msg_Debugging()<<"} -> as = "<<as<<" -> \\mu = "<<sqrt(mur2)<<"\n";
  return mur2;
}

void Fastjet_Scale_Setter::SetScale
(const std::string &mu2tag,Algebra_Interpreter &mu2calc)
{ 
  if (mu2tag=="" || mu2tag=="0") THROW(fatal_error,"No scale specified");
  msg_Debugging()<<METHOD<<"(): scale '"<<mu2tag
		 <<"' in '"<<p_proc->Name()<<"' {\n";
  msg_Indent();
  m_tagset.SetTags(&mu2calc);
  mu2calc.Interprete(mu2tag);
  if (msg_LevelIsDebugging()) mu2calc.PrintEquation();
  msg_Debugging()<<"}\n";
}

