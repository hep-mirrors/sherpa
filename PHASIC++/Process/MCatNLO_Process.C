#include "PHASIC++/Process/MCatNLO_Process.H"

#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "PHASIC++/Selectors/Jet_Finder.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PHASIC++/Process/ME_Generators.H"
#include "PHASIC++/Scales/Scale_Setter_Base.H"
#include "PHASIC++/Selectors/Combined_Selector.H"
#include "PHASIC++/Process/Single_Process.H"
#include "PHASIC++/Channels/BBar_Multi_Channel.H"
#include "PDF/Main/NLOMC_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/Jet_Criterion.H"
#include "PDF/Main/Cluster_Definitions_Base.H"
#include "PDF/Main/ISR_Handler.H"
#include "MODEL/Main/Running_AlphaS.H"
#include "ATOOLS/Math/MathTools.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/MyStrStream.H"
#include "ATOOLS/Org/Scoped_Settings.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Weight_Info.H"

#include <cassert>

using namespace ATOOLS;
using namespace PHASIC;
using namespace PDF;

MCatNLO_Process::MCatNLO_Process
(ME_Generators& gens,NLOTypeStringProcessMap_Map *pmap):
  m_gens(gens), p_bviproc(NULL), p_rsproc(NULL),
  p_bproc(NULL), p_rproc(NULL), p_ddproc(NULL),
  p_ampl(NULL)
{
  static bool ref(false);
  p_apmap=pmap;
  if (!ref) {
    ref=true;
    rpa->gen.AddCitation
      (1,"The automation of MCatNLO is published in \\cite{Hoeche:2011fd}.");
  }
  m_wassevent=false;
}

MCatNLO_Process::~MCatNLO_Process()
{
  if (p_ampl) p_ampl->Delete();
  if (p_rproc) delete p_rproc;
  if (p_bproc) delete p_bproc;
  if (p_ddproc) delete p_ddproc;
  if (p_rsproc) delete p_rsproc;
  if (p_bviproc) delete p_bviproc;
}

void MCatNLO_Process::Init(const Process_Info &pi,
			  BEAM::Beam_Spectra_Handler *const beam,
			   PDF::ISR_Handler *const isr,
			   YFS::YFS_Handler *const yfs,
			   const int mode)
{
  RegisterDefaults();
  Scoped_Settings s{ Settings::GetMainSettings()["MC@NLO"] };
  Process_Info cpi(pi);
  // book-keep user settings,
  // if not set, need Amegic for BVI and DADS, Comix for RS, B, R
  // (for historic reasons, still default to Amegic for B)
  std::string usermegen(pi.m_megenerator), userrsmegen(pi.m_rsmegenerator);
  std::string defbvimegen("Amegic"), defddmegen("Amegic");
  std::string defrsmegen("Comix"), defbmegen("Amegic"), defrmegen("Comix");
  cpi.m_fi.SetNLOType(nlo_type::born|nlo_type::loop|
		      nlo_type::vsub|nlo_type::real|nlo_type::rsub);
  if (pi.m_fi.m_nlocpl[0]==0. && pi.m_fi.m_nlocpl[1]==1.)
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
  else if (pi.m_fi.m_nlocpl[0]==1. && pi.m_fi.m_nlocpl[1]==0.)
    cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
  else THROW(not_implemented,"Cannot do NLO QCD+EW yet.");
  Process_Base::Init(cpi,beam,isr,yfs);
  m_pinfo.m_fi.SetNLOType(pi.m_fi.NLOType());
  Process_Info npi(m_pinfo);
  npi.m_fi.m_ps.pop_back();
  m_name=GenerateName(m_pinfo.m_ii,npi.m_fi);
  if (pi.Has(nlo_type::real)!=pi.Has(nlo_type::rsub))
    THROW(fatal_error, "R/S can't be initialised separately.");

  // init B proc
  Process_Info spi(pi);
  ++spi.m_fi.m_nmax;
  spi.m_fi.SetNLOType(cpi.m_fi.NLOType());
  Process_Info bpi(spi);
  if (bpi.m_megenerator=="") bpi.m_megenerator=defbmegen;
  p_bproc=InitProcess(bpi,nlo_type::lo,false);
  if (p_bproc==NULL) return;

  for (size_t i(0);i<m_pinfo.m_fi.m_nlocpl.size();++i) {
    spi.m_maxcpl[i]+=spi.m_fi.m_nlocpl[i];
    spi.m_mincpl[i]+=spi.m_fi.m_nlocpl[i];
  }

  // init R proc
  Process_Info rpi(spi);
  rpi.m_megenerator=rpi.m_rsmegenerator;
  if (rpi.m_megenerator=="") rpi.m_megenerator=defrmegen;
  p_rproc=InitProcess(rpi,nlo_type::lo,true);

  // init BVI proc
  nlo_type::code bvicode = (nlo_type::code) s["PP_BVI_MODE"].Get<int>();
  Process_Info bvipi(spi);
  if (bvipi.m_megenerator=="") bvipi.m_megenerator=defbvimegen;
  p_bviproc=InitProcess(bvipi, bvicode, false);

  // init DADS proc
  Process_Info ddpi(spi);
  if (ddpi.m_megenerator=="") ddpi.m_megenerator=defddmegen;
  p_ddproc=InitProcess(ddpi,nlo_type::rsub,1);

  // init RS proc
  Process_Info rspi(spi);
  rspi.m_integrator=rspi.m_rsintegrator;
  rspi.m_itmin=rspi.m_rsitmin;
  rspi.m_itmax=rspi.m_rsitmax;
  rspi.m_megenerator=rspi.m_rsmegenerator;
  if (rspi.m_megenerator=="") rspi.m_megenerator=defrsmegen;
  p_rsproc=InitProcess(rspi,nlo_type::real|nlo_type::rsub,1|2);

  p_rsproc->FillProcessMap(p_apmap);
  p_bviproc->FillProcessMap(p_apmap);
  p_ddproc->FillProcessMap(p_apmap);
  p_bproc->SetLookUp(false);
  p_rproc->SetLookUp(false);
  p_bproc->SetParent(this);
  p_rproc->SetParent(this);
  p_bproc->FillProcessMap(p_apmap);
  p_rproc->FillProcessMap(p_apmap);
  m_psmode   = s["PSMODE"].Get<int>();
  m_hpsmode  = s["HPSMODE"].Get<int>();
  if (m_hpsmode == -1) {
    if (m_pinfo.m_ckkw&1) m_hpsmode=0;
    else m_hpsmode=4;
  }
  m_kfacmode = s["KFACTOR_MODE"].Get<int>();
  m_fomode   = s["FOMODE"].Get<int>();
  m_rsscale  = s["RS_SCALE"].Get<std::string>();
  if (!m_fomode) {
    p_bviproc->SetSProc(p_ddproc);
    p_bviproc->SetMCMode(1);
    p_ddproc->SetMCMode(2);
    p_rsproc->SetMCMode(2);
  }
  if (p_rsproc->Size()!=p_rproc->Size()) {
    for (size_t i(0);i<p_rproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]R : "<<(*p_rproc)[i]->Name()<<std::endl;
    for (size_t i(0);i<p_rsproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]RS:  "<<(*p_rsproc)[i]->Name()<<std::endl;
    THROW(fatal_error,"R and RS have different size");
  }
  if (p_bproc->Size()!=p_bviproc->Size()) {
    for (size_t i(0);i<p_bproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]B:   "<<(*p_bproc)[i]->Name()<<std::endl;
    for (size_t i(0);i<p_bviproc->Size();++i)
      msg_Debugging()<<"["<<i<<"]BVI: "<<(*p_bviproc)[i]->Name()<<std::endl;
    THROW(fatal_error,"B and BVI have different size");
  }
  for (size_t i(0);i<p_rsproc->Size();++i)
    if ((*p_rsproc)[i]->Flavours()!=(*p_rproc)[i]->Flavours())
      THROW(fatal_error,"Ordering differs in R and RS");
  for (size_t i(0);i<p_bviproc->Size();++i)
    if ((*p_bviproc)[i]->Flavours()!=(*p_bproc)[i]->Flavours())
      THROW(fatal_error,"Ordering differs in B and BVI");

  for (size_t i(0);i<p_rsproc->Size();++i)
    (*p_rsproc)[i]->GetMEwgtinfo()->m_type=mewgttype::H;
}

void MCatNLO_Process::RegisterDefaults() const
{
  Scoped_Settings s{ Settings::GetMainSettings()["MC@NLO"] };
  s["HPSMODE"].SetDefault(-1);  // H event shower mode
  s["KFACTOR_MODE"].SetDefault(14);  // K-factor mode
  s["FOMODE"].SetDefault(0);  // fixed order mode
  s["RS_SCALE"].SetDefault("");  // RS scale
  s["PP_BVI_MODE"].SetDefault(7);  // BVI mode
}

Process_Base* MCatNLO_Process::InitProcess
(const Process_Info &pi,nlo_type::code nlotype,const int real)
{
  Process_Info cpi(pi);
  cpi.m_fi.SetNLOType(nlotype);
  if (real) {
    if (pi.m_fi.m_nlocpl[1]==0. && pi.m_fi.m_nlocpl[0]==1.)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
    else if (pi.m_fi.m_nlocpl[0]==0. && pi.m_fi.m_nlocpl[1]==1.)
      cpi.m_fi.m_ps.push_back(Subprocess_Info(kf_ewjet,"",""));
    else THROW(fatal_error, "Internal error.");
  }
  Process_Base* proc = m_gens.InitializeProcess(cpi,false);
  return proc;
}

bool MCatNLO_Process::InitSubtermInfo()
{
  for (size_t i(0);i<p_ddproc->Size();++i) {
    NLO_subevtlist *subs((*p_ddproc)[i]->GetSubevtList());
    for (size_t j(0);j<subs->size()-1;++j) {
      NLO_subevt *sub((*subs)[j]);
      for (size_t ij(0);ij<sub->m_n;++ij)
	for (size_t k(0);k<sub->m_n;++k)
          if (k != ij && sub->p_fl[k] == sub->p_fl[sub->m_kt] &&
              sub->p_fl[ij]==sub->p_fl[sub->m_ijt]) {
	    m_iinfo[sub->m_pname].insert(IDip_ID(ij,k));
	  }
      m_dinfo[subs->back()->m_pname].insert(*sub);
    }
  }
  return true;
}

Process_Base *MCatNLO_Process::FindProcess
(const Cluster_Amplitude *ampl,const nlo_type::code type,
 const bool error) const
{
  std::string name(Process_Base::GenerateName(ampl));
  StringProcess_Map::const_iterator pit(p_apmap->find(type)->second->find(name));
  if (pit!=p_apmap->find(type)->second->end()) {
    if (type==nlo_type::lo) {
      for (size_t i(0);i<p_bproc->Size();++i)
	if ((*p_bproc)[i]==pit->second) return pit->second;
      for (size_t i(0);i<p_rproc->Size();++i)
	if ((*p_rproc)[i]==pit->second) return pit->second;
    }
    if (type&nlo_type::vsub) {
      for (size_t i(0);i<p_bviproc->Size();++i)
	if ((*p_bviproc)[i]==pit->second) return pit->second;
    }
    if (type&nlo_type::rsub) {
      for (size_t i(0);i<p_rsproc->Size();++i)
	if ((*p_rsproc)[i]==pit->second) return pit->second;
    }
  }
  if (error)
    THROW(fatal_error,"Process '"+name+"'("+ToString(type)+") not found");
  return NULL;
}

Cluster_Amplitude *MCatNLO_Process::CreateAmplitude(const NLO_subevt *sub) const
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetNIn(m_nin);
  ampl->SetMS(p_rsproc->Generator());
  ampl->SetMuF2(sub->m_mu2[stp::fac]);
  ampl->SetMuR2(sub->m_mu2[stp::ren]);
  Int_Vector ci(sub->m_n,0), cj(sub->m_n,0);
  for (size_t i=0;i<sub->m_n;++i) {
    ampl->CreateLeg(i<m_nin?-sub->p_mom[i]:sub->p_mom[i],
		    i<m_nin?sub->p_fl[i].Bar():sub->p_fl[i],
		    ColorID(ci[i],cj[i]),sub->p_id[i]);
    if (!sub->IsReal() && sub->p_id[i]&(1<<sub->m_i)) {
      if ((sub->p_id[i]&(1<<sub->m_j))==0)
	THROW(fatal_error,"Internal error");
      ampl->Legs().back()->SetK(1<<sub->m_k);
    }
  }
  ampl->Decays()=*sub->p_dec;
  return ampl;
}

Weights_Map MCatNLO_Process::Differential(const Vec4D_Vector &p,
					  Variations_Mode varmode)
{
  THROW(fatal_error,"Invalid function call");
}

Weights_Map MCatNLO_Process::LocalKFactor(Cluster_Amplitude& ampl)
{
  DEBUG_FUNC(Name());

  // validate ampl
  Cluster_Amplitude *rampl(ampl.Prev());
  if (rampl->Legs().size()!=m_nin+m_nout) return 0.0;

  // evaluate RS process
  msg_Debugging()<<*rampl<<"\n";
  Process_Base *rsproc(FindProcess(rampl,nlo_type::rsub,false));
  if (rsproc==NULL) return 0.0;
  msg_Debugging()<<"Found '"<<rsproc->Name()<<"'\n";
  Cluster_Amplitude *crampl(rampl->Copy());
  int mm(p_bproc->Generator()->SetMassMode(0));
  p_bproc->Generator()->ShiftMasses(crampl);
  p_bproc->Generator()->SetMassMode(mm);
  int rmode = rampl->ColorMap().empty() ? 128 : 0;
  Weights_Map rs = rsproc->Differential(*crampl, Variations_Mode::all, rmode);
  Weights_Map r = rsproc->GetSubevtList()->back()->m_results;
  if (rmode && r.Nominal() == 0.0 &&
      rsproc->Differential(*crampl, Variations_Mode::nominal_only, 64).Nominal() != 0.0) {
    for (int i(0); i < 100 && r.Nominal() == 0.0; ++i) {
      rs = rsproc->Differential(*crampl, Variations_Mode::all, rmode);
      r = rsproc->GetSubevtList()->back()->m_results;
    }
  }
  crampl->Delete();
  msg_Debugging() << "H = " << rs.Nominal() << ", R = " << r.Nominal() << " -> "
                  << rs.Nominal() / r.Nominal() << "\n";
  if (r.Nominal() == 0.0) {
    return 0.0;
  }

  // evaluate BVI process
  msg_Debugging()<<ampl<<"\n";
  Process_Base *bviproc(FindProcess(&ampl,nlo_type::vsub,false));
  if (bviproc==NULL) return 0.0;
  msg_Debugging()<<"Found '"<<bviproc->Name()<<"'\n";
  Process_Base *bproc(FindProcess(&ampl));
  Cluster_Amplitude *campl(ampl.Copy());
  p_bproc->Generator()->SetMassMode(0);
  p_bproc->Generator()->ShiftMasses(campl);
  p_bproc->Generator()->SetMassMode(mm);
  bproc->GetMEwgtinfo()->m_type = mewgttype::none;
  int bmode = ampl.ColorMap().empty() ? 128 : 0;
  Weights_Map b = bproc->Differential(*campl, Variations_Mode::all, bmode);
  if (bmode && b.Nominal() == 0.0 &&
      bproc->Differential(*campl, Variations_Mode::nominal_only, 64).Nominal() != 0.0) {
    for (int i(0); i < 100 && b.Nominal() == 0.0; ++i) {
      b = bproc->Differential(*campl, Variations_Mode::all, bmode);
    }
  }
  if (b.Nominal() == 0.0) {
    campl->Delete();
    return 0.0;
  }
  bviproc->BBarMC()->GenerateEmissionPoint(*campl);
  Weights_Map bvi = bviproc->Differential(*campl, Variations_Mode::all, bmode);
  campl->Delete();

  // eventually calculate local K factor
  msg_Debugging() << "Calc'ing LocalKFactor weights.\n";
  const double random(ran->Get());

  // Below, we'll Iterate through all weight tuples, calculating the
  // corresponding K factors. For this, we use the absolute version of the
  // variations, which makes it more straightforward to deal with each
  // variation individually in the loop below.
  bvi.MakeAbsolute();
  b.MakeAbsolute();
  rs.MakeAbsolute();
  r.MakeAbsolute();

  // Copy BVI variations to get the same variations structure; the actual K
  // factor data is set in the loop below; note that via bvi.MakeAbsolute(), we
  // can be certain that global factors such as base_weight and
  // nominal_prefactor are set to 1.0
  Weights_Map kfacs = bvi;

  for (const auto& kv : bvi) {
    const std::string& key = kv.first;
    const Weights& weights = kv.second;
    const auto num_weights = weights.Size();
    for (int i {0}; i < num_weights; i++) {
      msg_Debugging() << weights.Name(i) << '\n';

      const LocalKFactorInfo info = CalculateLocalKFactorInfo(
          bvi.Get(key, i), b.Get(key, i), rs.Get(key, i), r.Get(key, i));
      if (info.s == 0.0 && info.h == 0.0) {
        kfacs[key][i] = 0.0;
        continue;
      }
      // select S or H and set corresponding K factor; update cluster
      // amplitudes if this is the nominal calculation
      const double selectionwgt {dabs(info.s) / (dabs(info.s) + dabs(info.h))};
      if (selectionwgt > random) {
        const double kfac {info.s / selectionwgt};
        msg_Debugging() << "S selected ( w = " << kfac << " )\n";
        if (key == "Main" && i == 0 && m_kfacmode / 10) {
          for (Cluster_Amplitude* campl(ampl.Next()); campl;
              campl = campl->Next()) {
            campl->SetLKF(bvi[key][i] / b[key][i]);
            campl->SetNLO(2);
          }
        }
        kfacs[key][i] = kfac;
      } else {
        const double kfac {info.h / (1.0 - selectionwgt)};
        msg_Debugging() << "H selected ( w = " << kfac << " )\n";
        if (key == "Main" && i == 0 && m_kfacmode / 10)
          ampl.SetNLO(m_hpsmode);
        kfacs[key][i] = kfac;
      }
      for (Cluster_Amplitude *campl(&ampl);campl;
	   campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
    }
  }
  kfacs.MakeRelative();
  return kfacs;
}

MCatNLO_Process::LocalKFactorInfo MCatNLO_Process::CalculateLocalKFactorInfo(
    double bvi, double b, double rs, double r)
{
  LocalKFactorInfo info;
  const double bvib(b ? bvi / b : 0.0), rsr(r ? rs / r : 0.0);
  if (m_kfacmode % 10 == 0) {
    info.s = bvib * (1.0 - rsr);
    info.h = rsr;
  } else if (m_kfacmode % 10 == 1) {
    info.s = bvib * (1.0 - rsr);
    info.h = 0;
  } else if (m_kfacmode % 10 == 2) {
    info.s = 0;
    info.h = rsr;
  } else if (m_kfacmode % 10 == 3) {
    info.s = bvib;
    info.h = 0.;
  } else if (m_kfacmode % 10 == 4) {
    info.s = bvib + (b?rs/b:0.);
    info.h = 0.;
  } else
    THROW(fatal_error, "Unknown Kfactor mode.");
  msg_Debugging() << "BVI = " << bvi << ", B = " << b << " -> S = " << info.s
                  << ", H = " << info.h << "\n";
  return info;
}

Cluster_Amplitude *MCatNLO_Process::GetAmplitude()
{
  if (p_ampl==NULL) return NULL;
  Cluster_Amplitude *ampl(p_ampl->CopyAll());
  ME_Generator_Base *gen((ME_Generator_Base*)ampl->MS());
  int mm(gen->SetMassMode(1));
  int stat(gen->ShiftMasses(ampl));
  if (stat<0) {
    msg_Tracking() << METHOD << "(): Mass shift failed." << std::endl;
    gen->SetMassMode(mm);
    return NULL;
  }
  if (stat==1) {
    stat=p_shower->GetClusterDefinitions()->ReCluster(ampl);
    if (stat!=1) {
      msg_Debugging()<<METHOD<<"(): Reclustering failed."<<std::endl;
    }
  }
  gen->SetMassMode(mm);
  return ampl;
}

Weights_Map MCatNLO_Process::OneHEvent(const int wmode)
{
  msg_Debugging()<<"H event\n";
  m_wassevent=false;
  if (p_ampl) p_ampl->Delete();
  p_selected=p_rproc;
  const size_t selectedindex(p_rproc->SynchronizeSelectedIndex(*p_rsproc));
  assert(selectedindex != std::numeric_limits<size_t>::max());
  Process_Base *rproc((*p_rproc)[selectedindex]);
  rproc->Integrator()->SetMax
    (p_rsproc->Selected()->Integrator()->Max());
  Vec4D_Vector &p(p_rsproc->Selected()->Integrator()->Momenta());
  rproc->SetFixedScale(p_rsproc->Selected()->ScaleSetter(1)->Scales());
  rproc->ScaleSetter(1)->CalculateScale(Vec4D_Vector());
  rproc->SetFixedScale(std::vector<double>());
  rproc->GetMEwgtinfo()->m_mur2=
    p_rsproc->Selected()->GetMEwgtinfo()->m_mur2;
  rproc->GetMEwgtinfo()->m_muf2=
    p_rsproc->Selected()->GetMEwgtinfo()->m_muf2;
  rproc->Integrator()->SetMomenta(p);
  Color_Integrator *ci(rproc->Integrator()->ColorIntegrator().get()),
    *rci(p_rsproc->Selected()->Integrator()->ColorIntegrator().get());
  if (ci && rci) {
    ci->SetI(rci->I());
    ci->SetJ(rci->J());
  }
  p_ampl=p_rsproc->Selected()->GetSubevtList()->back()->p_ampl;
  if (p_ampl) p_ampl=p_ampl->CopyAll();
  else p_ampl=((Single_Process*)(rproc))->Cluster(p,1);
  p_selected->Selected()->SetMEwgtinfo(*p_rsproc->Selected()->GetMEwgtinfo());
  if (p_ampl==NULL) {
    msg_Error()<<METHOD<<"(): No valid clustering. Skip event."<<std::endl;
    return {0};
  }
  Scale_Setter_Base *scs(p_rsproc->Selected()->ScaleSetter(1));
  for (Cluster_Amplitude *ampl(p_ampl);ampl;ampl=ampl->Next()) {
    ampl->SetMuR2(scs->Scale(stp::ren));
    ampl->SetMuF2(scs->Scale(stp::fac));
    ampl->SetMuQ2(scs->Scale(stp::res));
  }
  if (p_ampl->Next()) {
    p_ampl->Next()->SetNLO
      (p_ampl->Next()->NLO()|m_hpsmode);
  }
  Selector_Base *jf=p_rsproc->Selected()->
    Selector()->GetSelector("Jetfinder");
  Weights_Map wgtmap;
  wgtmap.insert({"QCUT", Weights {Variations_Type::qcut}});
  if (jf && m_nout-1<m_pinfo.m_fi.NMaxExternal()) {
    for (Cluster_Amplitude *ampl(p_ampl);
	 ampl;ampl=ampl->Next()) ampl->SetJF(jf);
    Jet_Finder *cjf(static_cast<Jet_Finder*>(jf));
    HEventVeto_Args hva(cjf,cjf->JC()->Value(p_ampl));
    for (size_t i {0}; i < s_variations->Size(Variations_Type::qcut) + 1;
         ++i) {
      const double fac {
          i != 0 ? s_variations->Qcut_Parameters(i - 1).m_scale_factor : 1.0};
      int stat = hva.m_jcv < sqr(hva.p_jf->Qcut() * fac);
      wgtmap["QCUT"][i] = (stat ? 1.0 : 0.0);
    }
  }
  msg_Debugging()<<"R selected via NLO "<<*p_ampl<<"\n";
  for (Cluster_Amplitude *campl(p_ampl->Next());campl;
       campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
  return wgtmap;
}

Weights_Map MCatNLO_Process::OneSEvent(const int wmode)
{
  msg_Debugging()<<"S event\n";
  m_wassevent=true;
  if (p_ampl) p_ampl->Delete();
  p_ampl=NULL;
  const size_t selectedindex(p_bproc->SynchronizeSelectedIndex(*p_bviproc));
  assert(selectedindex != std::numeric_limits<size_t>::max());
  Process_Base *bproc((*p_bproc)[selectedindex]);
  Vec4D_Vector &p(p_bviproc->Selected()->Integrator()->Momenta());
  p_ampl = dynamic_cast<Single_Process*>
    (p_bviproc->Selected())->Cluster(p);
  SortFlavours(p_ampl);
  p_ampl->SetProcs(p_apmap);
  p_ampl->SetIInfo(&m_iinfo);
  p_ampl->SetDInfo(&m_dinfo);
  p_ampl->Decays()=m_decins;
  if (p_ampl->JF<void>())
    for (Cluster_Amplitude *sampl(p_ampl);
	 sampl;sampl=sampl->Next()) sampl->SetNLO(2);
  double lkf(p_bviproc->Selected()->Last()/p_bviproc->Selected()->LastB());
  for (Cluster_Amplitude *ampl(p_ampl);
       ampl;ampl=ampl->Next()) ampl->SetLKF(lkf);
  int stat(1);
  if (!(m_psmode&2)) {
    p_nlomc->SetShower(p_shower);
    stat=p_nlomc->GeneratePoint(p_ampl);
  }
  Cluster_Amplitude *next(p_ampl), *ampl(p_ampl->Prev());
  if (ampl) {
    p_ampl=NULL;
    Process_Base *rproc(NULL);
    Flavour_Vector afl(ampl->Legs().size());
    for (size_t i(0);i<afl.size();++i)
      afl[i]=i<m_nin?ampl->Leg(i)->Flav().Bar():ampl->Leg(i)->Flav();
    for (size_t i(0);i<p_rproc->Size();++i)
      if (afl==(*p_rproc)[i]->Flavours()) {
	rproc=(*p_rproc)[i];
	break;
      }
    if (rproc==NULL) THROW(fatal_error,"Invalid splitting");
    p_selected=p_rproc;
    p_rproc->SetSelected(rproc);
    rproc->Integrator()->PSHandler()->SetEnhanceWeight(bproc->Integrator()->PSHandler()->EnhanceWeight());
    rproc->Integrator()->SetMax(bproc->Integrator()->Max());
    rproc->Integrator()->SetMomenta(*ampl);
    rproc->GetMEwgtinfo()->m_mur2=bproc->GetMEwgtinfo()->m_mur2;
    Cluster_Leg *lij(NULL);
    for (size_t i(0);i<next->Legs().size();++i) {
      next->Leg(i)->SetStat(1);
      if (next->Leg(i)->K()) {
	lij=next->Leg(i);
      }
    }
    if (!lij) THROW(fatal_error,"Internal error");
    std::vector<int> ids(ID(lij->Id()));
    size_t iid(1<<ids.front()), jid(1<<ids.back()), kid(lij->K());
    ids.push_back(ID(lij->K()).front());
    if (ids.size()!=3) THROW(fatal_error,"Internal error");
    Cluster_Amplitude *campl(ampl->Copy());
    campl->IdSort();
    Cluster_Definitions_Base *clus(p_shower->GetClusterDefinitions());
    Cluster_Param kt2=clus->Cluster
      (Cluster_Config(campl,ids[0],ids[1],ids[2],lij->Flav(),
		      p_bviproc->Generator(),NULL,1));
    if (kt2.m_kt2<=0.0) {
      kt2=clus->Cluster
	(Cluster_Config(campl,ids[1],ids[0],ids[2],lij->Flav(),
			p_bviproc->Generator(),NULL,1));
    }
    campl->Delete();
    ampl->SetKT2(kt2.m_kt2);
    ampl->SetMu2(kt2.m_kt2);
    p_ampl=ampl;
    ampl->SetMuF2(next->MuF2());
    ampl->SetMuR2(next->MuR2());
    ampl->SetOrderQCD(next->OrderQCD()+1);
    ampl->Next()->SetNLO(4);
    ampl->SetJF(ampl->Next()->JF<void>());
    ampl->Next()->SetCA(clus);
    next->SetKin(kt2.m_kin);
    while (ampl->Next()) {
      ampl=ampl->Next();
      if (!(ampl->Leg(0)->Id()&p_ampl->Leg(0)->Id()))
	std::swap<Cluster_Leg*>(ampl->Legs()[0],ampl->Legs()[1]);
      if (ampl->IdNew()&iid) ampl->SetIdNew(ampl->IdNew()|jid);
      for (size_t i(0);i<ampl->Legs().size();++i) {
	ampl->Leg(i)->SetNMax(p_ampl->Leg(i)->NMax());
	size_t cid(ampl->Leg(i)->Id());
	if (cid&iid) {
	  for (size_t j(0);j<ampl->Legs().size();++j)
	    if (ampl->Leg(j)->K()==cid)
	      ampl->Leg(j)->SetK(cid|jid);
	  ampl->Leg(i)->SetId(cid|jid);
	  if (ampl->Prev()->Prev()==NULL) {
	    ampl->Leg(i)->SetK(kid);
	  }
	}
      }
    }
    msg_Debugging()<<"R selected via Sudakov "<<*p_ampl
		   <<" ( w = "<<p_nlomc->WeightsMap().Nominal()<<" )\n";
    p_selected->Selected()->SetMEwgtinfo(*p_bviproc->Selected()->GetMEwgtinfo());
    return p_nlomc->WeightsMap();
  }
  p_selected=p_bproc;
  ampl=p_ampl;
  if (!(m_psmode&2)) ampl->SetNLO(4);
  else ampl->SetNLO(0);
  bproc->Integrator()->SetMomenta(*p_ampl);
  msg_Debugging()<<"B selected "<<*p_ampl
		 <<" ( w = "<<p_nlomc->WeightsMap().Nominal()<<" )\n";
  for (Cluster_Amplitude *campl(p_ampl->Next());campl;
       campl=campl->Next()) msg_Debugging()<<*campl<<"\n";
  p_selected->Selected()->SetMEwgtinfo(*p_bviproc->Selected()->GetMEwgtinfo());
  if (m_psmode&2) return Weights_Map{1.0};
  return stat ? p_nlomc->WeightsMap() : Weights_Map{0.0};
}

Weight_Info *MCatNLO_Process::OneEvent(const int wmode,
                                       Variations_Mode varmode,
                                       const int mode)
{
  DEBUG_FUNC("");
  const double S(p_bviproc->Integrator()->SelectionWeight(wmode));
  const double H(p_rsproc->Integrator()->SelectionWeight(wmode));
  Weight_Info *winfo(NULL);
  if (S > ran->Get() * (S + H)) {
    p_selected = p_bviproc;
    winfo = p_bviproc->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode == 0) {
      // calculate and apply weight factor
      const Weights_Map Swgts {OneSEvent(wmode)};
      assert(p_ampl);
      const double Bsel(p_bproc->Selected()->Integrator()->SelectionWeight(wmode));
      const double Ssel(p_bviproc->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Bsel / Ssel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Swgts;
      *(p_selected->Selected()->GetMEwgtinfo()) *= Swgts.Nominal()*wgtfac;
    }
  } else {
    p_selected = p_rsproc;
    winfo = p_rsproc->OneEvent(wmode, varmode, mode);
    if (winfo && m_fomode == 0) {
      // calculate and apply weight factor
      const Weights_Map Hwgts {OneHEvent(wmode)};
      const double Rsel(p_rproc->Selected()->Integrator()->SelectionWeight(wmode));
      const double RSsel(p_rsproc->Selected()->Integrator()->SelectionWeight(wmode));
      const double selwgtratio(Rsel / RSsel);
      const double wgtfac(selwgtratio);
      winfo->m_weightsmap *= wgtfac;
      winfo->m_weightsmap *= Hwgts;
      *p_selected->Selected()->GetMEwgtinfo() *= Hwgts.Nominal()*wgtfac;
    }
  }
  Mass_Selector *ms(Selected()->Generator());
  for (Cluster_Amplitude *ampl(p_ampl);
       ampl;ampl=ampl->Next()) {
    ampl->SetNLO(1|ampl->NLO());
    ampl->SetMS(ms);
  }
  if (winfo && winfo->m_weightsmap.IsZero()) {
    delete winfo;
    winfo=NULL;
  }
  if (winfo==NULL) return winfo;

  if (rpa->gen.HardSC() || (rpa->gen.SoftSC() && !Flavour(kf_tau).IsStable())) {
    msg_Debugging()<<"Calcing Differential for spin correlations using "
                   <<Selected()->Generator()->Name()<<":"<<std::endl;
    ME_Weight_Info mwi(*Selected()->GetMEwgtinfo());
    if (Selected()->Integrator()->ColorIntegrator()!=NULL)
      while (Selected()->Differential(*p_ampl,Variations_Mode::nominal_only,1|2|4|128).Nominal()==0.0);
    else
      Selected()->Differential(*p_ampl,Variations_Mode::nominal_only,1|2|4|128);
    Selected()->SetMEwgtinfo(mwi);
  }
  return winfo;
}

void MCatNLO_Process::InitPSHandler
(const double &maxerror,const std::string eobs,const std::string efunc)
{
  p_bviproc->InitPSHandler(maxerror,eobs,efunc);
  p_ddproc->InitPSHandler(maxerror,eobs,efunc);
  p_rsproc->InitPSHandler(maxerror,eobs,efunc);
  p_rsproc->Integrator()->SetEnhanceFactor
    (p_rsproc->Integrator()->EnhanceFactor()*p_int->RSEnhanceFactor());
  p_rproc->Integrator()->SetPSHandler
    (p_rsproc->Integrator()->PSHandler());
  p_bproc->Integrator()->SetPSHandler
    (p_bviproc->Integrator()->PSHandler());
}

bool MCatNLO_Process::CalculateTotalXSec(const std::string &resultpath,
                                         const bool         create)
{
  Vec4D_Vector p(p_rsproc->NIn()+p_rsproc->NOut());
  Cluster_Amplitude *ampl(Cluster_Amplitude::New());
  for (int i(0);i<p.size();++i)
    ampl->CreateLeg(Vec4D(),Flavour(kf_jet));
  auto psh = p_ddproc->Integrator()->PSHandler();
  do {
    psh->TestPoint(&p.front(),&p_ddproc->Info(),p_ddproc->Generator());
    for (size_t i(0);i<p.size();++i) ampl->Leg(i)->SetMom(p[i]);
    p_ddproc->Differential(*ampl,Variations_Mode::nominal_only,4);
  } while (!InitSubtermInfo());
  ampl->Delete();
  bool res(p_bviproc->CalculateTotalXSec(resultpath,create));
  psh=p_rsproc->Integrator()->PSHandler();
  if (psh->AbsError()==0.0)
    psh->SetAbsError(psh->Error()*rpa->Picobarn()*
		     dabs(p_bviproc->Integrator()->TotalResult()));
#ifndef USING__Threading
  if (!p_rsproc->CalculateTotalXSec(resultpath,create)) res=false;
#endif
  for (size_t i(0);i<p_bviproc->Size();++i)
    (*p_bproc)[i]->Integrator()->SetMax
      ((*p_bviproc)[i]->Integrator()->Max());
  return res;
}

void MCatNLO_Process::SetLookUp(const bool lookup)
{
  m_lookup=lookup;
  p_bviproc->SetLookUp(lookup);
  p_ddproc->SetLookUp(lookup);
  p_rsproc->SetLookUp(lookup);
}

bool MCatNLO_Process::InitScale()
{
  bool res(true);
  if (!p_bviproc->InitScale()) res=false;
  if (!p_ddproc->InitScale()) res=false;
  if (!p_rsproc->InitScale()) res=false;
  if (!p_rproc->InitScale()) res=false;
  if (!p_bproc->InitScale()) res=false;
  return res;
}

void MCatNLO_Process::SetScale(const Scale_Setter_Arguments &scale)
{
  p_bviproc->SetScale(scale);
  p_ddproc->SetScale(scale);
  if (m_rsscale!="") {
    Scale_Setter_Arguments rsscale(scale);
    rsscale.m_scale=m_rsscale;
    p_rsproc->SetScale(rsscale);
    p_rproc->SetScale(rsscale);
  }
  else {
    p_rsproc->SetScale(scale);
    p_rproc->SetScale(scale);
  }
  p_bproc->SetScale(scale);
}

void MCatNLO_Process::SetKFactor(const KFactor_Setter_Arguments &args)
{
  p_bviproc->SetKFactor(args);
  p_ddproc->SetKFactor(args);
  p_rsproc->SetKFactor(args);
  p_rproc->SetKFactor(args);
  p_bproc->SetKFactor(args);
}

void MCatNLO_Process::InitializeTheReweighting(ATOOLS::Variations_Mode mode)
{
  p_bviproc->InitializeTheReweighting(mode);
  p_ddproc->InitializeTheReweighting(mode);
  p_rsproc->InitializeTheReweighting(mode);
  p_rproc->InitializeTheReweighting(mode);
  p_bproc->InitializeTheReweighting(mode);
}

void MCatNLO_Process::SetFixedScale(const std::vector<double> &s)
{
  p_bviproc->SetFixedScale(s);
  p_ddproc->SetFixedScale(s);
  p_rsproc->SetFixedScale(s);
  p_rproc->SetFixedScale(s);
  p_bproc->SetFixedScale(s);
}

bool MCatNLO_Process::IsGroup() const
{
  return true;
}

size_t MCatNLO_Process::Size() const
{
  return 2;
}

Process_Base *MCatNLO_Process::operator[](const size_t &i)
{
  if (i==1) return p_rsproc;
  return p_bviproc;
}

void MCatNLO_Process::SetSelector(const Selector_Key &key)
{
  p_bviproc->SetSelector(key);
  p_ddproc->SetSelector(key);
  p_rsproc->SetSelector(key);
  p_rproc->SetSelector(key);
  p_bproc->SetSelector(key);
}

void MCatNLO_Process::SetShower(PDF::Shower_Base *const ps)
{
  p_shower=ps;
  p_bviproc->SetShower(ps);
  p_ddproc->SetShower(ps);
  p_rsproc->SetShower(ps);
  p_rproc->SetShower(ps);
  p_bproc->SetShower(ps);
}

void MCatNLO_Process::SetNLOMC(PDF::NLOMC_Base *const mc)
{
  p_nlomc=mc;
  p_bviproc->SetNLOMC(mc);
  p_ddproc->SetNLOMC(mc);
  p_rsproc->SetNLOMC(mc);
  p_rproc->SetNLOMC(mc);
  p_bproc->SetNLOMC(mc);
}
