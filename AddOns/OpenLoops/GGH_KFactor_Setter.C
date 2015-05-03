#include "ATOOLS/Phys/NLO_Subevt.H"
#include "ATOOLS/Phys/Flavour.H"

#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"

#include "AMEGIC++/Main/Single_Process_External.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"

#include "COMIX/Main/Single_Process.H"

#include "AddOns/OpenLoops/GGH_KFactor_Setter.H"
#include "AddOns/OpenLoops/GGH_Process_Manager.H"

using namespace PHASIC;
using namespace ATOOLS;

GGH_Process_Manager PHASIC::s_procmanager = GGH_Process_Manager();

GGH_KFactor_Setter::GGH_KFactor_Setter
(const KFactor_Setter_Arguments &args):
  KFactor_Setter_Base(args) 
{
  p_ampl = GetAmpl();
  m_set_vcorrection = false;
  if (p_proc->Name().find("__QCD(R)")==std::string::npos)
    m_real_corr=false;
  // This constructor will be called during runs in which only
  // libraries are written out. Switch off in that case
  if (p_proc->Generator()->Generators()->NewLibraries())
    {
      msg_Out() << METHOD <<
	": Libraries created, no KFactor will be applied in this run" << std::endl;
      SetOn(false);
    }
  else
    {
      s_procmanager.SetGenerators(p_proc->Generator()->Generators());
      p_default_tree = static_cast<AMEGIC::Single_Process*>
	(s_procmanager.GetProcess(*p_ampl, false));
      p_default_loop = static_cast<AMEGIC::Single_Process_External*>
	(s_procmanager.GetProcess(*p_ampl, true ));
    }
}
  
GGH_KFactor_Setter::~GGH_KFactor_Setter()
{
  if(p_ampl) p_ampl->Delete();
}

double GGH_KFactor_Setter::KFactor() 
{
  if(!m_on) return 1.;
  const Vec4D_Vector& p = GetMomenta();
  if(!m_real_corr || (p_ampl->Legs().size())<5)
    return MassCorrectionFactor(p);
  if(IsCollinear(p))
    return ClusterMassCorrectionFactor();
  else 
    return MassCorrectionFactor(p);
}

double GGH_KFactor_Setter::KFactor(const NLO_subevt& evt) 
{
  if(!m_on) return 1.;
  Vec4D_Vector p(evt.p_mom, &(evt.p_mom[evt.m_n]));
  for (size_t i(0);i<p_proc->NIn();++i) p[i]=-p[i];
  if(!m_real_corr || (p_ampl->Legs().size())<5)
    return MassCorrectionFactor(evt.m_pname, p);
  if(IsCollinear(Vec4D_Vector(evt.p_mom, &(evt.p_mom[evt.m_n]))))
    return ClusterMassCorrectionFactor();
  else 
    return MassCorrectionFactor(evt.m_pname, p);
}

Vec4D_Vector GGH_KFactor_Setter::GetMomenta() const 
{
  return p_proc->Integrator()->Momenta();
}

bool GGH_KFactor_Setter::IsCollinear(const Vec4D_Vector& p) const
{
  for (size_t i(3); i<p.size(); i++){
    if ( p[i].PPerp2() < IR_CO) return true;
    for (size_t j(i+1); j<p.size(); j++){
      if ( (fabs(p[j].PPerp(p[i])) < IR_CO)
	   || (fabs(p[i].PPerp(p[j])) < IR_CO))
	return true;
    }
  }
  return false;
}

double GGH_KFactor_Setter::OSVertexCorrection(){
  if(!m_set_vcorrection) SetOSVertexCorrection();
  return vertex_correction;
}

void GGH_KFactor_Setter::SetOSVertexCorrection(){
  double m_h = Flavour(kf_h0).Mass();
  Vec4D_Vector fake_p;
  Vec4D gmom0(m_h/2.,0,0,+m_h/2.);
  Vec4D gmom1(m_h/2.,0,0,-m_h/2.);
  Vec4D hmom(m_h,0,0,0);
  Cluster_Amplitude *amp = Cluster_Amplitude::New();
  amp->SetNIn(2);
  amp->CreateLeg(gmom0,Flavour(kf_gluon));
  amp->CreateLeg(gmom1,Flavour(kf_gluon));
  amp->CreateLeg(hmom, Flavour(kf_h0));
  fake_p.push_back(gmom0); fake_p.push_back(gmom1); fake_p.push_back(hmom);
  AMEGIC::Single_Process*          tree = static_cast<AMEGIC::Single_Process*>(s_procmanager.GetProcess(*amp, false));
  AMEGIC::Single_Process_External* loop = static_cast<AMEGIC::Single_Process_External*>(s_procmanager.GetProcess(*amp, true));
  vertex_correction = loop->DSigma(fake_p,false)/tree->DSigma(fake_p,false);
  m_set_vcorrection = true;
}

double GGH_KFactor_Setter::MassCorrectionFactor(const Vec4D_Vector& p)
{
  return p_default_loop->DSigma(p, false)/p_default_tree->DSigma(p, false);
}

double GGH_KFactor_Setter::MassCorrectionFactor(const Cluster_Amplitude& ampl)
{
  Vec4D_Vector p;
  p.push_back(-p_next_ampl->Leg(0)->Mom());
  p.push_back(-p_next_ampl->Leg(1)->Mom());
  for(size_t i(2); i<p_next_ampl->Legs().size(); i++){
    p.push_back(p_next_ampl->Leg(i)->Mom());
  }
  AMEGIC::Single_Process*          tree = static_cast<AMEGIC::Single_Process*>
    (s_procmanager.GetProcess(ampl, false));
  AMEGIC::Single_Process_External* loop = static_cast<AMEGIC::Single_Process_External*>
    (s_procmanager.GetProcess(ampl, true ));
  return loop->DSigma(p,false)/tree->DSigma(p,false);
}

double GGH_KFactor_Setter::MassCorrectionFactor(const std::string& name, const Vec4D_Vector& p)
{
  AMEGIC::Single_Process*          tree = static_cast<AMEGIC::Single_Process*>
    (s_procmanager.GetProcess(name, false));
  AMEGIC::Single_Process_External* loop = static_cast<AMEGIC::Single_Process_External*>
    (s_procmanager.GetProcess(name, true ));
  return loop->DSigma(p,false)/tree->DSigma(p,false);
}

double GGH_KFactor_Setter::MassCorrectionFactor(const NLO_subevt& evt){
  return MassCorrectionFactor(evt.m_pname,
			      Vec4D_Vector(evt.p_mom, &(evt.p_mom[evt.m_n])));
}

double GGH_KFactor_Setter::ClusterMassCorrectionFactor()
{
  SetNextAmplitude();
  if(!p_next_ampl){
      msg_Out() << METHOD << ": Warning, no cluster amplitude found for reweighting" << std::endl;
      msg_Out() << METHOD << ": Falling back to vertex correction" << std::endl; 
      return OSVertexCorrection();
    }
  if(p_next_ampl->Leg(2)->Mom().PPerp()< IR_CO){
    msg_Out() << METHOD << ": Falling back to vertex correction" << std::endl; 
    return OSVertexCorrection();
  }
  return MassCorrectionFactor(*p_next_ampl);
}

double GGH_KFactor_Setter::ClusterMassCorrectionFactor(NLO_subevtlist* subevts)
{
  if(!subevts->size()>1)
    THROW(fatal_error, "Internal error");
  NLO_subevtlist::const_iterator it=subevts->begin();
  NLO_subevt* select_proc= *subevts->begin();
  double minkt2 = (**it).m_kt2;
  for(; it!=subevts->end(); ++it)
    {
      if(dynamic_cast<AMEGIC::Single_Real_Correction*>(static_cast<Process_Base*>((*(it))->p_proc)))
	continue;
      if ((**it).m_kt2<minkt2)
	select_proc = *it;
    }
  if((select_proc->p_mom)[2].PPerp() < IR_CO){
    msg_Out() << METHOD << ": Falling back to vertex correction" << std::endl;
    return OSVertexCorrection();
  }
  else
    return MassCorrectionFactor(*select_proc);
}

Cluster_Amplitude* GGH_KFactor_Setter::GetAmpl() const {
  Cluster_Amplitude* ret(Cluster_Amplitude::New());
  ret->SetNIn(2);
  ret->CreateLeg(Vec4D(), p_proc->Flavours()[0].Bar());
  ret->CreateLeg(Vec4D(), p_proc->Flavours()[1].Bar());
  for(Flavour_Vector::const_iterator it=++++p_proc->Flavours().begin(); it!=p_proc->Flavours().end(); ++it)
    ret->CreateLeg(Vec4D(), *it);
  Process_Base::SortFlavours(ret);
  return ret;
}

void GGH_KFactor_Setter::SetNextAmplitude() {
  p_next_ampl = p_proc->ScaleSetter(1)->Amplitudes().front()->Next();
  Process_Base::SortFlavours(p_next_ampl);
}

bool GGH_KFactor_Setter::ContainsDecays(const Process_Base& proc) const {
  for(SubprocInfo_Vector::const_iterator it = p_proc->Info().m_fi.m_ps.begin();
      it != p_proc->Info().m_fi.m_ps.end(); ++it){
    if (it->GetExternal().size()>1)
      return true;
  }
  return false;
}

DECLARE_GETTER(GGH_KFactor_Setter,"GGH",
	       KFactor_Setter_Base,
	       KFactor_Setter_Arguments);

KFactor_Setter_Base *Getter <KFactor_Setter_Base,
					     KFactor_Setter_Arguments,
					     GGH_KFactor_Setter>::
operator()(const KFactor_Setter_Arguments &args) const
{
  return new GGH_KFactor_Setter(args);
}

void Getter<KFactor_Setter_Base,
		    KFactor_Setter_Arguments,
		    GGH_KFactor_Setter>::
PrintInfo(std::ostream &str,const size_t width) const 
{ str<<"GGH K-Factor\n"; }
