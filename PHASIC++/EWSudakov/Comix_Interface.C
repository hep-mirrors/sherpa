#include "PHASIC++/EWSudakov/Comix_Interface.H"

#include "PHASIC++/EWSudakov/EWSudakov_Amplitudes.H"
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

Comix_Interface::Comix_Interface(Process_Base* proc,
                                 EWSudakov_Amplitudes& ampls):
  p_proc{ proc }
{
  InitializeProcesses(ampls);
}

void Comix_Interface::FillSpinAmplitudes(
    std::vector<Spin_Amplitudes>& spinampls,
    ATOOLS::Cluster_Amplitude& ampl) const
{
  const auto loprocmapit = m_apmap.find(nlo_type::lo);
  if (loprocmapit == m_apmap.end())
    THROW(fatal_error, "LO entry in process map not found");
  Cluster_Amplitude* campl(ampl.Copy());
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  std::string pname(Process_Base::GenerateName(campl));
  auto pit = loprocmapit->second->find(pname);
  if (pit->second == NULL)
    THROW(fatal_error, "Process not found");
  pit->second->Differential(*campl, 2|4|128);
  campl->Delete();
  std::vector<std::vector<Complex>> cols;
  pit->second->FillAmplitudes(spinampls, cols);

  // TODO: understand the following sign corrections that we had to hard-code
  // in order to get amplitude ratios with the expected sign from COMIX

  // flip signs of all ee->ZP amplitudes
  if (ampl.Legs()[0]->Flav().Kfcode() == kf_e
      && ampl.Legs()[1]->Flav().Kfcode() == kf_e
      && ampl.Legs()[2]->Flav().Kfcode() == kf_Z
      && ampl.Legs()[3]->Flav().Kfcode() == kf_photon) {
    for (size_t i {0}; i < spinampls[0].size(); ++i) {
      const auto& spins = spinampls[0].GetSpinCombination(i);
      spinampls[0][i] = -spinampls[0][i];
    }
  }

  // flip signs of all e- ve -> mu- vmu amplitudes
  else if (ampl.Legs()[0]->Flav().Kfcode() == kf_e
             && !ampl.Legs()[0]->Flav().IsAnti()
             && ampl.Legs()[1]->Flav().Kfcode() == kf_nue
             && ampl.Legs()[2]->Flav().Kfcode() == kf_mu
             && ampl.Legs()[3]->Flav().Kfcode() == kf_numu) {
    for (size_t i {0}; i < spinampls[0].size(); ++i) {
      const auto& spins = spinampls[0].GetSpinCombination(i);
      spinampls[0][i] = -spinampls[0][i];
    }
  }

  // flip signs of the e- ve -> Z/h0 W- amplitudes when the bosons are
  // longitudinally polarised
  else if (ampl.Legs()[0]->Flav().Kfcode() == kf_e
      && !ampl.Legs()[0]->Flav().IsAnti()
      && ampl.Legs()[1]->Flav().Kfcode() == kf_nue
      && (ampl.Legs()[2]->Flav().Kfcode() == kf_Z || ampl.Legs()[2]->Flav().Kfcode() == kf_h0)
      && ampl.Legs()[3]->Flav().Kfcode() == kf_Wplus) {
    for (size_t i {0}; i < spinampls[0].size(); ++i) {
      const auto& spins = spinampls[0].GetSpinCombination(i);
      if ((spins[2] == 2 || spins[2] == 0) && spins[3] == 2)
        spinampls[0][i] = -spinampls[0][i];
    }
  }
}

void Comix_Interface::InitializeProcesses(EWSudakov_Amplitudes& ampls)
{
  DEBUG_FUNC("");
  const auto gpath =
    Default_Reader{}.Get<std::string>("PRINT_EWSUDAKOV_GRAPHS", "");
  for (auto& kv : ampls) {
    auto& ampl = kv.second;
    msg_Debugging() << "Initialize process for ampl=" << *ampl << std::endl;
    std::string name(PHASIC::Process_Base::GenerateName(ampl.get())
                     + "__Sudakov");
    Process_Info pi;
    pi.m_addname="__Sudakov";
    pi.m_megenerator="Comix";
    if (gpath != "")
      pi.m_gpath=gpath;
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
    pi.m_maxcpl=p_proc->Info().m_maxcpl;
    pi.m_mincpl=p_proc->Info().m_mincpl;
    PHASIC::Process_Base *proc=
      p_proc->Generator()->Generators()->InitializeProcess(pi,false);
    if (proc==NULL) THROW(fatal_error,"Invalid process");
    Selector_Key skey(NULL,NULL,true);
    proc->SetSelector(skey);
    proc->SetScale
      (Scale_Setter_Arguments
       (MODEL::s_model,"VAR{"+ToString(sqr(rpa->gen.Ecms()))+"}","Alpha_QCD 1"));
    proc->SetKFactor(KFactor_Setter_Arguments("None"));
    proc->Get<COMIX::Process_Base>()->Tests();
    proc->FillProcessMap(&m_apmap);
  }
}
