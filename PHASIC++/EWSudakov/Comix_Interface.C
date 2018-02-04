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
                                 Cluster_Amplitude& ampl) :
  m_proc{ proc },
  m_ampl{ ampl }
{ }

Comix_Interface::~Comix_Interface()
{ }

/*
** Initialize Comix processes
*/

void Comix_Interface::FillSpinAmplitudes(
    std::vector<Spin_Amplitudes>& spinampls)
{
  DEBUG_VAR(m_ampl);
  NLOTypeStringProcessMap_Map *procs
    (m_ampl.Procs<NLOTypeStringProcessMap_Map>());
  if (procs==NULL)
    THROW(fatal_error, "Process map not found");
  nlo_type::code type=nlo_type::lo;
  if (procs->find(type)==procs->end())
    THROW(fatal_error, "LO entry in process map not found");
  Cluster_Amplitude *campl(m_ampl.Copy());
  campl->SetMuR2(sqr(rpa->gen.Ecms()));
  campl->SetMuF2(sqr(rpa->gen.Ecms()));
  campl->SetMuQ2(sqr(rpa->gen.Ecms()));
  Process_Base::SortFlavours(campl);
  std::string pname(Process_Base::GenerateName(campl));
  StringProcess_Map::const_iterator pit((*(*procs)[type]).find(pname));
  if (pit==(*(*procs)[type]).end()) {
    (*(*procs)[type])[pname]=NULL;
    pit=(*procs)[type]->find(pname);
  }
  if (pit->second==NULL) {
    THROW(fatal_error, "Process not found");
  }
  // TODO: Somehow, this puts event generation into an infinite loop
  // do we need it for diagonal terms at all?
  //int kfon(pit->second->KFactorSetter(true)->On());
  //pit->second->KFactorSetter(true)->SetOn(false);
  //pit->second->Differential(*campl,2|4|128);
  //pit->second->KFactorSetter(true)->SetOn(kfon);
  campl->Delete();
  std::vector<std::vector<Complex> > cols;
  pit->second->FillAmplitudes(spinampls,cols);
}
