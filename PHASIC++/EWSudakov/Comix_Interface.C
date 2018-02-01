#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "COMIX/Main/Single_Process.H"
#include "ATOOLS/Phys/Color.H"

#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Shell_Tools.H"
#include "ATOOLS/Org/Message.H"

using namespace PHASIC;
using namespace COMIX;
using namespace ATOOLS;

Comix_Interface::Comix_Interface(Process_Base * const proc,
                 Cluster_Amplitude * ampl) :
  p_proc(proc),p_ampl(ampl)
{
  InitializeProcesses();
}

Comix_Interface::~Comix_Interface()
{ }

/*
** Initialize Comix processes
*/

PHASIC::Process_Base * Comix_Interface::
GetComixProc(Cluster_Amplitude * const ampl)
{
  PHASIC::Process_Base::SortFlavours(ampl);
  std::string name(PHASIC::Process_Base::GenerateName(ampl)+"__Sudakov");
  StringProcess_Map *pm(m_apmap[nlo_type::lo]);
  StringProcess_Map::const_iterator pit(pm->find(name));
  if (pit==pm->end()) THROW(fatal_error,"Initialization failed");
  return pit->second;
}

void Comix_Interface::InitializeProcesses()
{
  if(m_apmap.find(nlo_type::lo)==m_apmap.end())
    m_apmap[nlo_type::lo] = new StringProcess_Map();
  StringProcess_Map *pm(m_apmap[nlo_type::lo]);
  PHASIC::Process_Base::SortFlavours(p_ampl);
  std::string name(PHASIC::Process_Base::GenerateName(p_ampl)+"__Sudakov");
  StringProcess_Map::const_iterator pit(pm->find(name));
  PHASIC::Process_Base *xs=NULL;
  if (pit!=pm->end() && pit->second->
      Integrator()->ColorIntegrator()!=NULL) xs=pit->second;
  if (xs) return;
  //My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  //My_In_File::OpenDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
  Process_Info pi;
  pi.m_addname="__Sudakov";
  pi.m_megenerator="Comix";
  for (size_t i(0);i<p_ampl->NIn();++i) {
    Flavour fl(p_ampl->Leg(i)->Flav().Bar());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  for (size_t i(p_ampl->NIn());i<p_ampl->Legs().size();++i) {
    Flavour fl(p_ampl->Leg(i)->Flav());
    if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
    pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
  }
  pi.m_maxcpl=p_proc->Info().m_maxcpl;
  pi.m_mincpl=p_proc->Info().m_mincpl;
  DEBUG_VAR(pi.m_maxcpl<<" " <<pi.m_mincpl);
  PHASIC::Process_Base *proc=
    p_ampl->Proc<PHASIC::Process_Base>()->
    Generator()->Generators()->InitializeProcess(pi,false);
  if (proc==NULL) THROW(fatal_error,"Invalid process");
  Selector_Key skey(NULL,NULL,true);
  proc->SetSelector(skey);
  proc->SetScale
    (Scale_Setter_Arguments
     (MODEL::s_model,"VAR{"+ToString(sqr(rpa->gen.Ecms()))+"}","Alpha_QCD 1"));
  proc->SetKFactor(KFactor_Setter_Arguments("None"));
  proc->Get<COMIX::Process_Base>()->Tests();
  proc->FillProcessMap(&m_apmap);
  //My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Comix/");
  //My_In_File::CloseDB(rpa->gen.Variable("SHERPA_CPP_PATH")+"/Process/Sherpa/");
  if ((pit=pm->find(name))==pm->end()) THROW(fatal_error,"Internal error");
  xs=pit->second;
  if (xs==NULL) THROW(fatal_error,"Invalid process");
}

void Comix_Interface::Construct_Ampli(Complex ampl1, Complex ampl2)
{
  Cluster_Amplitude *ca_ampl1, *ca_ampl2;
  /* construct them in some way starting from p_ampl ... */
  ampl1 = amplitude(ca_ampl1);
  ampl2 = amplitude(ca_ampl2);
}

Complex Comix_Interface::amplitude(ATOOLS::Cluster_Amplitude *const ampl)
{
  Complex res;
  COMIX::Single_Process *xs(GetComixProc(ampl)->Get<COMIX::Single_Process>());
  DEBUG_FUNC(xs->Name());
  msg_Debugging()<<*ampl<<"\n";
  Vec4D_Vector p(ampl->Legs().size());
  for (size_t i(0);i<p.size();++i)
    p[i]=i<ampl->NIn()?-ampl->Leg(i)->Mom():ampl->Leg(i)->Mom();
  xs->GetAmplitude()->SetMomenta(p);
  std::vector<double> s(xs->ScaleSetter(1)->Scales().size(),0.0);
  s[stp::fac]=p_proc->ScaleSetter()->Scale(stp::fac);
  s[stp::ren]=p_proc->ScaleSetter()->Scale(stp::ren);
  s[stp::res]=p_proc->ScaleSetter()->Scale(stp::res);
  xs->SetFixedScale(s);
  xs->ScaleSetter(1)->CalculateScale(p);
  std::vector<Spin_Amplitudes> amps;
  std::vector<std::vector<Complex> > cols;
  xs->FillAmplitudes(amps,cols);
  res = amps.front().front();
  return res;
}

double Comix_Interface::Differential(ATOOLS::Cluster_Amplitude *const ampl)
{
  /* 1st create the two ampitudes based on the emission of a W, Z or a photon.
   * then we have ampl1 * conj(ampl2). Later we multiply by the relevant
   * Sudakov log (in Sudakov.C), and the coupling factor in KFactor.C
   */
  Complex ampl1, ampl2;
  Construct_Ampli(ampl1,ampl2);
  return (ampl1*std::conj(ampl2)).real();
}

double Comix_Interface::Born(){
  return (amplitude(p_ampl)*std::conj(amplitude(p_ampl))).real();
}
