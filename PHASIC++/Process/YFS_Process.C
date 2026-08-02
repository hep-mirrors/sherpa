#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/NLO_Types.H"
#include "PHASIC++/Process/Process_Info.H"
#include "PHASIC++/Process/YFS_Process.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "MODEL/Main/Single_Vertex.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Hard_Process_Variation_Generator.H"
#include "ATOOLS/Phys/Weight_Info.H"
#include "SHERPA/SoftPhysics/Resonance_Finder.H"


#include <cassert>
#include <cstddef>

using namespace ATOOLS;
using namespace PHASIC;
using namespace PDF;

using namespace SHERPA;
using namespace MODEL;
using namespace std;

YFS_Process::YFS_Process
(ME_Generators& gens, NLOTypeStringProcessMap_Map *pmap):
  m_gens(gens), p_bornproc(NULL), p_realproc(NULL), p_realrealproc(NULL),
  p_realvirtproc(NULL),
  p_yfsvirt(NULL), p_yfsreal(NULL), p_yfsrealvirt(NULL), p_yfsrealreal(NULL),
  p_yfsvv(NULL)
{
  m_lastxs= 0.0;
  RegisterDefaults();
  static bool ref(false);
  // m_gens.LoadGenerator(m_name);
  p_apmap = pmap;
}

YFS_Process::~YFS_Process() {
  if (p_bornproc) delete p_bornproc;
  if (p_realproc) delete p_realproc;
  if (p_virtproc) delete p_virtproc;
  if (p_realrealproc) delete p_realrealproc;
  if (p_realvirtproc) delete p_realvirtproc;
  if (p_vv) delete p_vv;
  if (p_int) delete p_int;
  // ME providers for this process's YFS corrections; the shared NLO_Base only
  // borrows them, so they are released here.
  if (p_yfsvirt) delete p_yfsvirt;
  if (p_yfsreal) delete p_yfsreal;
  if (p_yfsrealvirt) delete p_yfsrealvirt;
  if (p_yfsrealreal) delete p_yfsrealreal;
  if (p_yfsvv) delete p_yfsvv;
  // p_yfs is NOT deleted: every process is handed the one handler owned by
  // Initialization_Handler, so deleting it here was a double free as soon as
  // there was more than one process.
  if (p_apmap) delete p_apmap;
  // if (p_rampl) delete p_rampl;
  if (p_gen) delete p_gen;
}

void YFS_Process::Init(const Process_Info &pi,
                       BEAM::Beam_Spectra_Handler *const beam,
                       PDF::ISR_Handler *const isr,
                       YFS::YFS_Handler *const yfs, const int mode)
{
  for (size_t beam=0;beam<2;++beam) {
      yfs->SetInFlav(pi.ExtractFlavours()[beam]);
      // yfs->SetBeam(p_beamspectra);
  }
  m_gens.InitializeProcess(pi, true);
  p_yfs = yfs;
  PRINT_VAR(pi.ExtractFlavours());
  p_yfs->resetparticles();
  p_yfs->SetFlavours(pi.ExtractFlavours());
  Process_Info ypi(pi), vpi(pi);
  if (pi.m_fi.m_nlocpl[0] != 0) THROW(not_implemented, "YFS cannot do NLO QCD.");
  if (pi.Has(nlo_type::rsub) || pi.Has(nlo_type::vsub)) {
    THROW(fatal_error, "YFS subtraction terms cannot be seperated. Only use BVR in NLO_Part");
  }
  m_name = GenerateName(ypi.m_ii, ypi.m_fi);
  Process_Base::Init(ypi, beam, isr, yfs);
  p_bornproc = InitProcess(ypi, nlo_type::born, false);
  p_yfs->SetNLOType(nlo_type::born);
  if (pi.Has(nlo_type::real)) {
    Process_Info rpi(pi); // real process info
    for (size_t i = 0; i < pi.m_fi.m_nlocpl.size(); ++i)
    {
      rpi.m_maxcpl[i] += rpi.m_fi.m_nlocpl[i];
    }
    rpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon, "", ""));
    p_realproc = InitProcess(rpi, nlo_type::real, false);
    // p_realproc->SetLookUp(false);
    p_realproc->SetParent(this);
    // p_realproc->FillProcessMap(p_apmap);
    // p_realproc->SetLookUp(true);
    if (p_yfs->NLO()->NeedsRealProvider()) {
      p_yfsreal = new YFS::Real(rpi);
      p_yfsreal->SetProc(p_realproc);
    }
    p_yfs->SetNLOType(nlo_type::real);
  }
  if (pi.Has(nlo_type::loop)) {
    vpi.m_fi.SetNLOType(nlo_type::loop);
    // Process_Base::Init(vpi, beam, isr, yfs);
    p_virtproc = InitProcess(vpi, nlo_type::loop, false);
    p_virtproc->FillProcessMap(p_apmap);
    // p_virtproc->p_mapproc=NULL
    p_virtproc->SetParent(this);
    p_virtproc->SetLookUp(false);
    if (p_yfs->NLO()->NeedsVirtualProvider()) {
      p_yfsvirt = new YFS::Virtual(vpi);
      p_yfsvirt->SetProc(p_bornproc);
    }
    p_yfs->SetNLOType(nlo_type::loop);
  }
  if (pi.Has(nlo_type::rvirt)) {
    Process_Info rvpi(pi);
    for (size_t i = 0; i < pi.m_fi.m_nlocpl.size(); ++i)
    {
      rvpi.m_maxcpl[i] += rvpi.m_fi.m_nlocpl[i];
      rvpi.m_mincpl[i] += rvpi.m_fi.m_nlocpl[i];
    }
    rvpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon, "", ""));
    p_realvirtproc = InitProcess(rvpi, nlo_type::rvirt, false);
    p_yfsrealvirt = new YFS::RealVirtual(rvpi);
    p_yfsrealvirt->SetProc(p_realvirtproc);
    p_yfs->SetNLOType(nlo_type::loop);
  }
  if (pi.Has(nlo_type::realreal)) {
    Process_Info rrpi(pi);
    for (size_t i = 0; i < pi.m_fi.m_nlocpl.size(); ++i)
    {
      rrpi.m_maxcpl[i] += rrpi.m_fi.m_nlocpl[i]+rrpi.m_fi.m_nlocpl[i];
      rrpi.m_mincpl[i] += rrpi.m_fi.m_nlocpl[i]+rrpi.m_fi.m_nlocpl[i];
    }
    rrpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon, "", ""));
    rrpi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon, "", ""));
    p_realrealproc = InitProcess(rrpi, nlo_type::real, false);
    if (!p_yfsreal) {
      msg_Error() << "No real corrections matrix element provider found.\n"
                  << "    Double-real (RR) subtraction terms for NNLO cannot be\n"
                  << "    constructed without them. Check that the real correction\n"
                  << "    process has been set up and linked correctly before\n"
                  << "    initialising the NNLO subtraction.\n";
      THROW(fatal_error, "Missing real corrections provider for NNLO RR subtraction");
    }
    p_yfsrealreal = new YFS::RealReal(rrpi);
    p_realrealproc->SetParent(this);
    p_yfsrealreal->SetProc(p_realrealproc);
    p_yfs->SetNLOType(nlo_type::realreal);
  }
  if(pi.Has(nlo_type::vv)) {
    Process_Info vvpi(pi);
    for(size_t i=0;  i < pi.m_fi.m_nlocpl.size(); ++i)
    {
      vvpi.m_maxacpl[i] += vvpi.m_fi.m_nlocpl[i]+vvpi.m_fi.m_nlocpl[i];
      vvpi.m_mincpl[i] += vvpi.m_fi.m_nlocpl[i]+vvpi.m_fi.m_nlocpl[i];
    }
    p_vv = InitProcess(vvpi, nlo_type::vv, false);
    p_yfsvv = new YFS::VirtualVirtual(vvpi);
  }
  p_bornproc->SetLookUp(false);
  // p_bornproc->SetParent(p_bornproc);
  p_bornproc->SetSelected(this);
  // The handler is shared with every other process, so hand it this process's
  // providers now and again whenever this process becomes the selected one.
  MakeActive();
  FindResonances();
}

void YFS_Process::MakeActive()
{
  // Everything the shared YFS_Handler holds that is process specific is
  // (re)pointed here: the flavours/particles and the NLO ME providers. Keep
  // these two in step -- a handler carrying one process's flavours and
  // another's loop ME is how the multi-process breakage showed up.
  p_yfs->resetparticles();
  p_yfs->SetFlavours(Flavours());
  p_yfs->NLO()->SetProviders(p_yfsvirt, p_yfsreal, p_yfsrealvirt,
                             p_yfsrealreal, p_yfsvv);
}


Process_Base* YFS_Process::InitProcess
(const Process_Info &pi, nlo_type::code nlotype, const bool real)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlotype);
  Process_Base* proc = m_gens.InitializeProcess(pi, false);
  if (!proc)
  {
    std::stringstream msg;
    msg_Error()<< METHOD <<": Unable to initialize process:\n" << cpi;
    THROW(fatal_error,  msg.str());
  }
  return proc;
}

bool YFS_Process::CalculateTotalXSec(const std::string &resultpath,
                                     const bool create)
{
  p_int = p_bornproc->Integrator();
  p_int->Reset();
  auto psh = p_int->PSHandler();
  MakeActive();
  psh->InitCuts();
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  double var(p_int->TotalVar());
  msg_Info() << METHOD << "(): Calculate xs for '"
             // << m_name << "' (" << (p_gen ? p_gen->Name() : "") << ")" << std::endl;
             << m_name <<  std::endl;
  double totalxsborn(psh->Integrate() / rpa->Picobarn());
  if (!IsEqual(totalxsborn, p_int->TotalResult())) {
    msg_Error() << "Result of PS-Integrator and summation do not coincide!\n"
                << "  '" << m_name << "': " << totalxsborn
                << " vs. " << p_int->TotalResult() << std::endl;
    }
  if (p_int->Points()) {
    p_int->SetTotal();
    if (var == p_int->TotalVar()) {
      exh->RemoveTerminatorObject(p_int);
      return 1;
    }
    p_int->StoreResults();
    exh->RemoveTerminatorObject(p_int);
    return 1;
    }
  exh->RemoveTerminatorObject(p_int);
  return 0;
}


 ATOOLS::Weights_Map YFS_Process::Differential(const Vec4D_Vector& p,
                                         Variations_Mode varmode) {
  THROW(fatal_error,"Invalid function call");
}

void YFS_Process::InitPSHandler(const double &maxerror,
                                const std::string eobs,
                                const std::string efunc) {
  if (!p_bornproc) return;
  p_bornproc->InitPSHandler(maxerror, eobs, efunc);
  p_bornproc->Integrator()->SetPSHandler
  (p_bornproc->Integrator()->PSHandler());
}


Cluster_Amplitude *YFS_Process::CreateAmplitude(const NLO_subevt *sub) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(m_nin);
  PRINT_INFO(p_realproc->Generator());
  ampl->SetMS(p_realproc->Generator());
  ampl->SetMuF2(sub->m_mu2[stp::fac]);
  ampl->SetMuR2(sub->m_mu2[stp::ren]);
  Int_Vector ci(sub->m_n, 0), cj(sub->m_n, 0);
  for (size_t i = 0; i < sub->m_n; ++i) {
    ampl->CreateLeg(i < m_nin ? -sub->p_mom[i] : sub->p_mom[i],
                    i < m_nin ? sub->p_fl[i].Bar() : sub->p_fl[i],
                    ColorID(ci[i], cj[i]), sub->p_id[i]);
    if (!sub->IsReal() && sub->p_id[i] & (1 << sub->m_i)) {
      if ((sub->p_id[i] & (1 << sub->m_j)) == 0)
        THROW(fatal_error, "Internal error");
      ampl->Legs().back()->SetK(1 << sub->m_k);
    }
  }
  ampl->Decays() = *sub->p_dec;
  return ampl;
}




Weight_Info *YFS_Process::OneEvent(const int wmode,ATOOLS::Variations_Mode varmode,
                                   const int mode)
{
  auto psh = p_int->PSHandler();
  psh->InitCuts();
  MakeActive();
  p_selected = p_bornproc;
  Weight_Info *winfo(NULL);
  winfo = p_int->PSHandler()->OneEvent(this, varmode, mode);
  // PRINT_VAR(METHOD);
  // if(p_realproc) OneRealEvent();
  return winfo;
}


void YFS_Process::FindResonances() {
  std::map<std::string, MODEL::Vertex_List> restab_map;
  Vertex_List vlist;
  FindProcessPossibleResonances(m_flavs, vlist);
  msg_Debugging() << "Process: " << this->Name() << " -> "
            << vlist.size() << " non-QCD resonances.\n";
  for (size_t k = 0; k < vlist.size(); ++k) msg_Out() << vlist[k] << endl;
  restab_map[this->Name()] = vlist;
  p_yfs->p_dipoles->SetProcResMap(restab_map);
  if(p_yfs->NLO()) p_yfs->NLO()->p_nlodipoles->SetProcResMap(restab_map);
}

void YFS_Process::OneRealEvent(){
  // Vec4D_Vector &p(p_yfs);
  // Cluster_Amplitude *rampl;
  // Cluster_Amplitude(rampl);
  Weight_Info *winfo(NULL);
}



void YFS_Process::FindProcessPossibleResonances
(const Flavour_Vector& fv, MODEL::Vertex_List& vlist)
{
  const Vertex_Table *vtab(s_model->VertexTable());
  Flavour_Vector fslep;
  for (size_t i(2); i < fv.size(); ++i) {
    if (!fv[i].Strong()) fslep.push_back(fv[i]);
  }
  for (Vertex_Table::const_iterator it(vtab->begin()); it != vtab->end(); ++it) {
    if (it->first.IsOn()      && !it->first.Strong() &&
        it->first.IsMassive() && !it->first.IsDummy()) {
      for (size_t i(0); i < it->second.size(); ++i) {
        bool on(true);
        double m(it->first.Mass());
        if(it->first.Width()==0) {on = false; break;}
        Single_Vertex * v(it->second[i]);
        for (size_t j(1); j < v->in.size(); ++j) {
          if (v->dec)        { on = false; break; }
          if (v->in[j] == v->in[0].Bar()) { on = false; break; }
          if (v->in[j].IsDummy())      { on = false; break; }
          if ((m -= v->in[j].Mass()) < 0.) { on = false; break; }
          bool flavfound(false);
          for (size_t k(0); k < fslep.size(); ++k)
            if (v->in[j] == fslep[k])    { flavfound = true; break; }
          if (!flavfound)              { on = false; break; }
        }
        if (on) vlist.push_back(v);
      }
    }
  }
}


void YFS_Process::SetLookUp(const bool lookup) {
  return;
}

bool YFS_Process::InitScale() {
  return true;
}
void YFS_Process::SetScale(const Scale_Setter_Arguments &scale) {
  if (p_bornproc) p_bornproc->SetScale(scale);
  if (p_realproc) p_realproc->SetScale(scale);
  if (p_realrealproc) p_realrealproc->SetScale(scale);
}


void YFS_Process::SetKFactor(const KFactor_Setter_Arguments &args) {
  if (p_bornproc) p_bornproc->SetKFactor(args);
}


void YFS_Process::RegisterDefaults() {
  Scoped_Settings s{ Settings::GetMainSettings()["ISR_YFS"] };
  s["QEDMODE"].SetDefault(0);
}


void YFS_Process::SetFixedScale(const std::vector<double> &s)
{
  if (p_bornproc) p_bornproc->SetFixedScale(s);
}

bool YFS_Process::IsGroup() const
{
  return false;
}

size_t YFS_Process::Size() const
{

  return 1;
}

Process_Base *YFS_Process::operator[](const size_t &i)
{
  if (i == 1) return p_realproc;
  return p_bornproc;
}

void YFS_Process::SetSelector(const Selector_Key &key)
{
  if (p_bornproc) p_bornproc->SetSelector(key);
  if (p_realproc) p_realproc->SetSelector(key);
  if (p_realrealproc) p_realrealproc->SetSelector(key);
  if (p_realvirtproc) p_realvirtproc->SetSelector(key);
}

void YFS_Process::SetGenerator(ME_Generator_Base *const gen)
{
  PRINT_INFO("HERE");
  p_gen = gen;
}



void YFS_Process::SetShower(PDF::Shower_Base *const ps)
{
  p_shower = ps;
  if (p_bornproc) p_bornproc->SetShower(ps);
  if (p_realproc) p_realproc->SetShower(ps);
  if (p_realrealproc) p_realrealproc->SetShower(ps);
}

void YFS_Process::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc = mc;
  if (p_bornproc) p_bornproc->SetNLOMC(mc);
  if (p_realproc) p_realproc->SetNLOMC(mc);
  if (p_realrealproc) p_realrealproc->SetNLOMC(mc);
}

void YFS_Process::InitializeTheReweighting(ATOOLS::Variations_Mode mode)
{
  if (mode == Variations_Mode::nominal_only)
    return;
  // Parse settings for hard process variation generators; note that this can
  // not be done in the ctor, because the matrix element handler has not yet
  // fully configured the process at this point, and some variation generators
  // might require that (an example would be EWSud)
  Settings& s = Settings::GetMainSettings();
  for (auto varitem : s["VARIATIONS"].GetItems()) {
    if (varitem.IsScalar()) {
      const auto name = varitem.Get<std::string>();
      if (name == "None") {
        return;
      } else {
        m_hard_process_variation_generators.push_back(
            Hard_Process_Variation_Generator_Getter_Function::
            GetObject(name, Hard_Process_Variation_Generator_Arguments{this})
            );
      }
    }
  }
}
