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

Comix_Interface::Comix_Interface(Process_Base& proc,
                                 ATOOLS::Cluster_Amplitude* ampl) :
  m_proc{ proc },
  p_ampl{ ampl }
{
  // m_proc.FillProcessMap(&m_apmap);
  InitializeSU2RotatedProcesses();
  InitializeSU2RotatedProcesses(1);
}

void Comix_Interface::FillSpinAmplitudes(
    std::vector<Spin_Amplitudes>& spinampls,
    ATOOLS::Cluster_Amplitude* ampl) const
{
  const auto loprocmapit{ m_apmap.find(nlo_type::lo) };
  if (loprocmapit == m_apmap.end())
    THROW(fatal_error, "LO entry in process map not found");
  Cluster_Amplitude* campl(ampl->Copy());
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  Process_Base::SortFlavours(campl); 
  std::string pname(Process_Base::GenerateName(campl));
  StringProcess_Map::const_iterator
    pit(loprocmapit->second->find(pname));
  if (pit->second==NULL) 
    THROW(fatal_error, "Process not found");
    
  pit->second->Differential(*campl,2|4|128);
  campl->Delete();
  std::vector<std::vector<Complex> > cols;
  pit->second->FillAmplitudes(spinampls,cols);
}

void Comix_Interface::InitializeSU2RotatedProcesses(const size_t mode)
{
  for (size_t i{ 0 }; i < p_ampl->Legs().size(); ++i) {
    auto* ampl = p_ampl->Copy();
    if(mode == 1 ) { // rotate 
      auto* leg = ampl->Leg(i);
      auto flav = leg->Flav();
      Flavour newflav;
      // TODO: generalise to other flavours (preferredly, add an WeakIsoPartner
      // function to the Flavour class
      if (flav.IsPhoton()) {
	newflav = Flavour{kf_Z};
      } else if (flav.Kfcode() == kf_Z) {
	newflav = Flavour{kf_photon};
      } else {
	ampl->Delete();
	continue;
      }
      leg->SetFlav(newflav);
    }
    PHASIC::Process_Base::SortFlavours(ampl);
    std::string name(PHASIC::Process_Base::GenerateName(ampl)+"__Sudakov");
    Process_Info pi;
    pi.m_addname="__Sudakov";
    pi.m_megenerator="Comix";
    for (size_t i(0);i<ampl->NIn();++i) {
      Flavour fl(ampl->Leg(i)->Flav().Bar());
      if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
      pi.m_ii.m_ps.push_back(Subprocess_Info(fl,"",""));
    }
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      Flavour fl(ampl->Leg(i)->Flav());
      if (Flavour(kf_jet).Includes(fl)) fl=Flavour(kf_jet);
      pi.m_fi.m_ps.push_back(Subprocess_Info(fl,"",""));
    }
    pi.m_maxcpl=m_proc.Info().m_maxcpl;
    pi.m_mincpl=m_proc.Info().m_mincpl;
    PHASIC::Process_Base *proc=
      m_proc.Generator()->Generators()->InitializeProcess(pi,false);
    if (proc==NULL) THROW(fatal_error,"Invalid process");
    Selector_Key skey(NULL,NULL,true);
    proc->SetSelector(skey);
    proc->SetScale
      (Scale_Setter_Arguments
       (MODEL::s_model,"VAR{"+ToString(sqr(rpa->gen.Ecms()))+"}","Alpha_QCD 1"));
    proc->SetKFactor(KFactor_Setter_Arguments("None"));
    proc->Get<COMIX::Process_Base>()->Tests();
    proc->FillProcessMap(&m_apmap);
    ampl->Delete();
    ampl = NULL;
  }
}
