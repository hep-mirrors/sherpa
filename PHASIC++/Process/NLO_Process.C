#include "PHASIC++/Process/NLO_Process.H"

#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Exception.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "AMEGIC++/Main/Single_Process.H"
#include "AMEGIC++/DipoleSubtraction/Single_Real_Correction.H"
#include "AMEGIC++/DipoleSubtraction/Single_DipoleTerm.H"
#include "PHASIC++/Selectors/Selector.H"
#include "PHASIC++/Main/Process_Integrator.H"
#include "PHASIC++/Main/Phase_Space_Handler.H"
#include "PHASIC++/Channels/Extra_Emission_Generator.H"
#include "PHASIC++/Process/ME_Generator_Base.H"
#include "PDF/Main/Shower_Base.H"
#include "PDF/Main/ISR_Handler.H"

using namespace std;
using namespace ATOOLS;
using namespace PHASIC;

NLO_Process::NLO_Process(ME_Generators& gens) :
  m_gens(gens),
  p_bproc(NULL), p_vproc(NULL), p_iproc(NULL), p_rproc(NULL), p_sproc(NULL)
{
  THROW(not_implemented, "Not implemented.");
}

void NLO_Process::Init(const Process_Info &pi,
                       BEAM::Beam_Spectra_Handler *const beam,
                       PDF::ISR_Handler *const isr)
{
  Process_Base::Init(pi, beam, isr);

  if (! (m_pinfo.m_fi.NLOType()&nlo_type::born &&
         m_pinfo.m_fi.NLOType()&nlo_type::real)) {
    THROW(fatal_error, "Born and real part must be initialised.");
  }

  if (m_pinfo.m_fi.NLOType()&nlo_type::born) p_bproc=InitBornLike(nlo_type::lo);
  if (m_pinfo.m_fi.NLOType()&nlo_type::loop) p_vproc=InitBornLike(nlo_type::loop);
  if (m_pinfo.m_fi.NLOType()&nlo_type::vsub) p_iproc=InitBornLike(nlo_type::vsub);
  if (m_pinfo.m_fi.NLOType()&nlo_type::real) p_rproc=InitRealLike(nlo_type::lo);
  if (m_pinfo.m_fi.NLOType()&nlo_type::rsub) p_sproc=InitRealLike(nlo_type::rsub);

  m_bxsecs.resize(p_bproc->Size());

  msg_Info()<<endl<<METHOD<<": "<<"Initialised NLO process for "<<Name()<<endl;
  if (p_bproc) {
    msg_Info()<<"  "<<p_bproc->Name()<<" "<<p_bproc->Generator()->Name()
              <<" with "<<p_bproc->Size()<<" subprocess(es):"<<endl;
    for (size_t i=0;i<p_bproc->Size();++i) {
      msg_Info()<<"    "<<(*p_bproc)[i]->Name();
      if ((*p_bproc)[i]->IsMapped())
        msg_Info()<<" -> "<<(*p_bproc)[i]->Get<AMEGIC::Process_Base>()->Partner()->Name();
      msg_Info()<<endl;
    }
  }
  if (p_vproc) {
    msg_Info()<<"  "<<p_vproc->Name()<<" "<<p_vproc->Generator()->Name()
              <<" with "<<p_vproc->Size()<<" subprocess(es),"<<endl;
    for (size_t i=0;i<p_vproc->Size();++i) {
      msg_Info()<<"    "<<(*p_vproc)[i]->Name();
      if ((*p_vproc)[i]->IsMapped())
        msg_Info()<<" -> "<<(*p_vproc)[i]->Get<AMEGIC::Process_Base>()->Partner()->Name();
      msg_Info()<<endl;
    }
  }
  if (p_iproc) {
    msg_Info()<<"  "<<p_iproc->Name()<<" "<<p_iproc->Generator()->Name()
              <<" with "<<p_iproc->Size()<<" subprocess(es),"<<endl;
    for (size_t i=0;i<p_iproc->Size();++i) {
      msg_Info()<<"    "<<(*p_iproc)[i]->Name();
      if ((*p_iproc)[i]->IsMapped())
        msg_Info()<<" -> "<<(*p_iproc)[i]->Get<AMEGIC::Process_Base>()->Partner()->Name();
      msg_Info()<<endl;
    }
  }
  if (p_rproc) {
    msg_Info()<<"  "<<p_rproc->Name()<<" "<<p_rproc->Generator()->Name()
              <<" with "<<p_rproc->Size()<<" subprocess(es),"<<endl;
    for (size_t i=0;i<p_rproc->Size();++i) {
      msg_Info()<<"    "<<(*p_rproc)[i]->Name();
      if ((*p_rproc)[i]->IsMapped())
        msg_Info()<<" -> "<<(*p_rproc)[i]->Get<AMEGIC::Process_Base>()->Partner()->Name();
      msg_Info()<<endl;
    }
  }
  if (p_sproc) {
    msg_Info()<<"  "<<p_sproc->Name()<<" "<<p_sproc->Generator()->Name()
              <<" with "<<p_sproc->Size()<<" subprocess(es)."<<endl;
    for (size_t i=0;i<p_sproc->Size();++i) {
      msg_Info()<<"    "<<(*p_sproc)[i]->Name();
      if ((*p_sproc)[i]->IsMapped())
        msg_Info()<<" -> "<<(*p_sproc)[i]->Get<AMEGIC::Process_Base>()->Partner()->Name();
      msg_Info()<<endl;
    }
  }

  // fill maps needed for POWHEG
  for (size_t i=0; i<p_bproc->Size(); ++i) {
    m_pmap_born[(*p_bproc)[i]->Name()]=(*p_bproc)[i];
  }
  for (size_t i=0; i<p_rproc->Size(); ++i) {
    m_pmap_real[(*p_rproc)[i]->Name()]=(*p_rproc)[i];
  }
  FillRBMap();
  msg_Info()<<"rbmap:  for "<<m_rbmap.size()<<" real processes"<<endl;
  for (size_t i=0;i<m_rbmap.size();++i) {
    msg_Info()<<" "<<m_rbmap[i].size()<<":  "<<(*p_rproc)[i]->Name();
    if ((*p_rproc)[i]->IsMapped())
      msg_Info()<<" -> "<<(*p_rproc)[i]->Get<AMEGIC::Single_Process>()->Partner()->Name();
    msg_Info()<<"\n      ";
    for (size_t j=0;j<m_rbmap[i].size();++j) {
      msg_Info()<<m_rbmap[i][j]<<" ,  ";
    }
    msg_Info()<<endl;
  }
}

Process_Base* NLO_Process::InitBornLike(nlo_type::code nlotype)
{
  Process_Info pi(m_pinfo);
  pi.m_fi.SetNLOType(nlotype);
  Process_Base *proc(m_gens.InitializeProcess(pi,false));
  if (proc==NULL) {
    msg_Error()<<pi<<std::endl;
    THROW(not_implemented, "Process not found.");
  }
  else {
    proc->SetParent(this);
    return proc;
  }
}

Process_Base* NLO_Process::InitRealLike(nlo_type::code nlotype)
{
  Process_Info pi(m_pinfo);
  pi.m_fi.SetNLOType(nlotype);
  if (pi.m_fi.m_nloqcdtype==nlotype) {
    pi.m_fi.m_ps.push_back(Subprocess_Info(kf_jet,"",""));
  }
  else if (pi.m_fi.m_nloewtype==nlotype) {
    pi.m_fi.m_ps.push_back(Subprocess_Info(kf_photon,"",""));
  }
  else {
    THROW(fatal_error, "Internal error.");
  }
  Process_Base *proc(m_gens.InitializeProcess(pi,false));
  if (proc==NULL) {
    msg_Error()<<pi<<std::endl;
    THROW(not_implemented, "Process not found.");
  }
  else {
    proc->SetParent(this);
    return proc;
  }
}


double NLO_Process::Differential(const Vec4D_Vector &p)
{
  m_last=0.0;

  m_lastsingle=p_bproc->Differential(p);
  if (m_lastsingle==0.0) return 0.0;
  // fill B xsecs
  for (size_t i=0; i<p_bproc->Size(); ++i) {
    m_bxsecs[i]=(*p_bproc)[i]->Last();
  }

  if (p_vproc) m_last+=p_vproc->Differential(p);
  if (p_iproc) m_last+=p_iproc->Differential(p);

  if (p_rproc && p_sproc) {
    double rs_sum(0.0), asw(0.0);
    while (rs_sum==0.0) {
      // create extra emission phase space point
      Cluster_Amplitude *ampl(NULL), *next(NULL);
      while (true) {
        ampl=CreateAmplitude(this,p);
        if (!p_int->EEG()->GeneratePoint(ampl)) {
          ampl->Delete();
          msg_Error()<<METHOD<<"(): EEG failure. Discard point."<<std::endl;
          return m_last=0.0;
        }
        next=ampl->Next();
        p_int->EEG()->GenerateWeight(next);
        if (p_int->EEG()->Weight()>0.0) break;
        ampl->Delete();
      }
      // calculate dsigma(R-S) for N+1 phase space
      if (p_rproc->IsGroup() && p_sproc->IsGroup()) {
        Process_Base* winner(NULL);
        double wmax(0.0);
        // only (possibly) add a single flavour R-S ?
        for (size_t i(0);i<p_rproc->Size();++i) {
          Cluster_Amplitude * bam(NULL), * ram(NULL);
          // calculate subtraction terms
          double valrsub=(*p_sproc)[i]->Differential(*next);
          // if some subtraction terms accepted event, such that the
          // underlying born is defined, add real contribution
          if (valrsub!=0.0) {
            double valreal=(*p_rproc)[i]->Differential(*next);
            // further calculate the (RB_ME)/(RB_PS) ratio for every
            // born momentum configuration
            // only if shower is set
            if (p_shower && false) { // disable R/B calculation until it works
              // extract momenta of N+1 phase space
              ATOOLS::Vec4D_Vector pNLO(next->Legs().size(),Vec4D(0.,0.,0.,0.));
              for (size_t j=0;j<next->Legs().size();++j)
                pNLO[j] = next->Leg(j)->Mom();
              // create clus. ampl. for real
              ram = CreateAmplitude((*p_rproc)[i],pNLO);
              // calculate RB for all possible born processes
              for (size_t j=0;j<m_rbmap[i].size();++j) {
                // set LO momenta with which to calculate the born again
                AMEGIC::Single_DipoleTerm *
                    sdt((*p_sproc)[i]->Get<AMEGIC::Single_Real_Correction>()
                                     ->Partner()
                                     ->Get<AMEGIC::Single_Real_Correction>()
                                     ->GetSubTerm(j));
                ATOOLS::Vec4D_Vector pLO(sdt->ActiveMom());
                // create clus. ampl. for born and set real as next
                bam = CreateAmplitude((*p_bproc)[m_rbmap[i][j]],pLO);
                bam->SetNext(ram);
                SetIds(bam,ram,sdt);
                // calculate RB_ME for every possible clustering
                (*p_bproc)[m_rbmap[i][j]]->SetLookUp(false);
                double bme((*p_bproc)[m_rbmap[i][j]]->Differential(pLO));
                  (*p_bproc)[m_rbmap[i][j]]->SetLookUp(true);
                double rbme(valreal*p_int->EEG()->Weight()/bme);
                // calculate RB_PS for every possible clustering
                // take active moms, ask shower for its weight for this splitting
                // take amplitude with its next
                // set which one is the new leg (-> ask Frank)
                // ask shower for trial weight
                double rbps(p_shower->TrialWeight(bam));
                // fill into the integrator of the born
                (*p_bproc)[m_rbmap[i][j]]->Integrator()->AddRBPoint(rbme/rbps);
                // output of (RB_ME/RB_PS)
                msg_Info()<<"me : "<<rbme<<"   :    "
                          <<(*p_rproc)[i]->Name()<<" / "
                          <<(*p_bproc)[m_rbmap[i][j]]->Name()<<std::endl;
                msg_Info()<<"ps : "<<rbps<<"   :    "
                          <<(*p_rproc)[i]->Name()<<" / "
                          <<(*p_bproc)[m_rbmap[i][j]]->Name()<<std::endl;
                msg_Info()<<"(R/B)_ME / (R/B)_PS : "<<rbme/rbps<<std::endl;
              }
            }
            if (valreal>wmax) {
              wmax=valreal;
              winner=(*p_rproc)[i];
            }
            valrsub+=valreal;
          }
          rs_sum+=valrsub;
        }
        // reweight alpha_s with winning process' kT^2
        if (winner && p_shower) {
          ME_Generator_Base *gen(winner->Generator());
          gen->SetClusterDefinitions(p_shower->GetClusterDefinitions());
          Cluster_Amplitude *ampl(gen->ClusterConfiguration(winner,1));
          asw=p_shower->CouplingWeight(ampl);
        }
      }
      else {
        // todo
        rs_sum=p_sproc->Differential(*next);
        if (rs_sum!=0.0) rs_sum+=p_rproc->Differential(*next);
      }
//       THROW(normal_exit,"manual exit");
      ampl->Delete();
      if (rs_sum==0.0) AddPoint(0.0);
    }
    m_last+=rs_sum*p_int->EEG()->Weight();
    m_last*=asw; // this probably shouldn't be done
  }
  m_last+=m_lastsingle;
  return m_last;
}

double NLO_Process::Differential2()
{
  m_last=0.0;
  
  if (p_bproc) m_last+=p_bproc->Differential2();
  if (p_vproc) m_last+=p_vproc->Differential2();
  if (p_iproc) m_last+=p_iproc->Differential2();
  if (p_rproc) m_last+=p_rproc->Differential2();
  if (p_sproc) m_last+=p_sproc->Differential2();
  
  return m_last;
}

double NLO_Process::RealEmissionWeight(const Cluster_Amplitude * realampl)
{
  Process_Base::SortFlavours((Cluster_Amplitude*)realampl);
  const Cluster_Amplitude * bornampl(realampl->Prev());
  Process_Base * born = FindBorn(bornampl);
  Process_Base * real = FindReal(realampl);
  born->SetKFactorOn(false);
  real->SetKFactorOn(false);
  double wb(born->Differential(*bornampl));
  double wr(real->Differential(*realampl));
  born->SetKFactorOn(true);
  real->SetKFactorOn(true);
  return wb/wr;
}

void NLO_Process::FillRBMap()
{
  DEBUG_FUNC("");
  for (size_t i=0; i<p_sproc->Size(); ++i) {
    m_rbmap.push_back(std::vector<size_t>());
    AMEGIC::Single_Real_Correction *
             src((*p_sproc)[i]->Get<AMEGIC::Single_Real_Correction>());
    AMEGIC::Single_Real_Correction *
            srcp((*p_sproc)[i]->Get<AMEGIC::Single_Real_Correction>()->Partner()
                              ->Get<AMEGIC::Single_Real_Correction>());
    std::vector<std::pair<Flavour,Flavour> > mappedflavs;
    if (src->IsMapped()) {
      DEBUG_INFO("ismapped");
      ATOOLS::Flavour_Vector procflavs(src->Flavours());
      ATOOLS::Flavour_Vector partflavs(srcp->Flavours());
      for (size_t j=0;j<procflavs.size();++j) {
        if (procflavs[j] != partflavs[j]) {
          mappedflavs.push_back(make_pair(partflavs[j],procflavs[j]));
        }
      }
    }
    DEBUG_INFO("flavs for "<<src->Name()<<" to be mapped are:");
    DEBUG_INFO("          "<<srcp->Name());
    for (size_t j=0;j<mappedflavs.size();++j)
      msg_Info()<<mappedflavs[j].first<<" "<<mappedflavs[j].second<<endl;
    for (size_t j=0; j<srcp->GetSubTermNumber(); ++j) {
      AMEGIC::Single_DipoleTerm * sdtp(srcp->GetSubTerm(j));
      if (sdtp->IsValid()) {
        // rename the un.-born to find the right born
        Process_Info unbornpi(sdtp->GetLOProcess()->Info());
        unbornpi.m_fi.SetNLOType(nlo_type::lo);
        // replace the flavours as they are mapped
        if (src->IsMapped()) {
          ATOOLS::Flavour_Vector inflavs;
          unbornpi.m_ii.GetExternal(inflavs);
          ATOOLS::Flavour_Vector outflavs;
          unbornpi.m_fi.GetExternal(outflavs);
          for (size_t k=0;k<mappedflavs.size();++k) {
            for (size_t l=0;l<inflavs.size();++l) {
              if (inflavs[l]==mappedflavs[k].first) {
                unbornpi.m_ii.SetExternal(mappedflavs[k].second,l);
                msg_Debugging()<<"replaced "<<mappedflavs[k].first
                               <<" with "<<mappedflavs[k].second
                               <<" in initial state at position "<<l<<std::endl;
              }
            }
            for (size_t l=0;l<outflavs.size();++l) {
              if (outflavs[l]==mappedflavs[k].first) {
                unbornpi.m_fi.SetExternal(mappedflavs[k].second,l);
                msg_Debugging()<<"replaced "<<mappedflavs[k].first
                               <<" with "<<mappedflavs[k].second
                               <<" in final state at position "<<l<<std::endl;
              }
            }
          }
        }
        SortFlavours(unbornpi);
        DEBUG_INFO(unbornpi);
        std::string name(GenerateName(unbornpi.m_ii,unbornpi.m_fi));
        DEBUG_VAR(name);
        for (unsigned int k=0;k<p_bproc->Size();++k) {
          if ((*p_bproc)[k]->Name() == name) {
            m_rbmap.back().push_back(k);
            DEBUG_INFO("found born of name '"<<(*p_bproc)[k]->Name()<<"'");
            break;
          }
          if (k==p_bproc->Size()-1)
            THROW(fatal_error,"born of name '"+name+"' not recovered");
        }
      }
    }
  }
}

Process_Base * NLO_Process::FindBorn(const Cluster_Amplitude * born)
{
  std::string name(Process_Base::GenerateName(born));
  Process_Map::iterator pit(m_pmap_born.find(name));
  if (pit != m_pmap_born.end()) return pit->second;
  THROW(fatal_error,"Born process '"+name+"' not found");
  return NULL;
}

Process_Base * NLO_Process::FindReal(const Cluster_Amplitude * real)
{
  std::string name(Process_Base::GenerateName(real));
  Process_Map::iterator pit(m_pmap_real.find(name));
  if (pit != m_pmap_real.end()) return pit->second;
  THROW(fatal_error,"Real process '"+name+"' not found");
  return NULL;
}

Cluster_Amplitude *NLO_Process::CreateAmplitude(const Process_Base * proc,
                                                const Vec4D_Vector &p)
{
  Cluster_Amplitude *ampl = Cluster_Amplitude::New();
  ampl->SetMS(proc->Generator());
  ampl->SetNIn(proc->NIn());
  Int_Vector ci(p.size(),0), cj(p.size(),0);
  SP(Color_Integrator) cint(proc->Integrator()->ColorIntegrator());
  if (cint!=NULL) {
    ci=cint->I();
    cj=cint->J();
    for (size_t i(0);i<ci.size();++i)
      ampl->ColorMap()[ci[i]]=ci[i];
  }
  for (size_t i=0;i<proc->NIn();++i)
    ampl->CreateLeg(-p[i],proc->Info().m_ii.GetExternal(i).Bar(),
		    ColorID(ci[i],cj[i]));
  for (size_t i=proc->NIn();i<proc->NIn()+proc->NOut();++i)
    ampl->CreateLeg(p[i],proc->Info().m_fi.GetExternal(i-proc->NIn()),
		    ColorID(ci[i],cj[i]));
  return ampl;
}

void NLO_Process::SetIds(ATOOLS::Cluster_Amplitude *bam,
                         ATOOLS::Cluster_Amplitude *ram,
                         AMEGIC::Single_DipoleTerm *sdt)
{
  DEBUG_VAR("");
  // set IDs for em/spec
  bam->Leg(sdt->Lijt())->SetId((1<<sdt->Li())|(1<<sdt->Lj()));
  bam->Leg(sdt->Lkt())->SetId(1<<sdt->Lk());
  ram->SetIdNew((1<<sdt->Lj()));
  // reset IDs for uninvolved partons
  // -> ordering unchanged, but possibly moved backwards
  // careful: exclusion for uninvolved partons of same
  //          name needs to be implemented
  for (size_t k=0;k<bam->Legs().size();++k) {
    if ((k != sdt->Lijt()) && (k != sdt->Lkt())) {
      for (size_t l=0;l<ram->Legs().size();++l) {
        if ((l != sdt->Li()) && (l != sdt->Lj()) &&
            (l != sdt->Lk()) &&
            (bam->Leg(k)->Flav() == ram->Leg(l)->Flav())) {
          bam->Leg(k)->SetId(ram->Leg(l)->Id());
          break;
        }
        if (l == ram->Legs().size()-1) {
          msg_Debugging()<<*bam<<"\n"<<*ram<<"\n"<<k<<"\n";
          THROW(fatal_error,"unmatched leg found : "+ToString(k));
        }
      }
    }
  }
  DEBUG_VAR(*bam);
  DEBUG_VAR(*ram);
  DEBUG_VAR(sdt->Lijt()<<" -> "<<sdt->Li()<<" "<<sdt->Lj()<<" , "<<
            sdt->Lkt()<<" -> "<<sdt->Lk());
}


size_t NLO_Process::Size() const
{
  return p_bproc->Size();
}

Process_Base *NLO_Process::operator[](const size_t &i)
{
  return (*p_bproc)[i];
}

void NLO_Process::DeSelect()
{
  p_selected=NULL;
}

bool NLO_Process::SelectOne()
{
  double wborn_max=0.0;
  for (size_t i=0; i<m_bxsecs.size(); ++i) {
    double wborn=m_bxsecs[i];
    if (wborn>wborn_max) {
      wborn_max=wborn;
      p_selected=(*p_bproc)[i];
    }
  }
  p_selected->Integrator()->SetSelectionWeight
    (p_int->SelectionWeight()/m_bxsecs.size());
  p_selected->Integrator()->SetMomenta(p_int->Momenta());
  return true;
}

Weight_Info *NLO_Process::OneEvent() 
{
  p_selected=this;
  Weight_Info* nloweight=p_int->PSHandler()->OneEvent(this);
  SelectOne();
  return nloweight;
}

Weight_Info *NLO_Process::WeightedEvent(const int mode) 
{
  p_selected=this;
  Weight_Info* nloweight=p_int->PSHandler()->WeightedEvent(this,mode);
  SelectOne();
  return nloweight;
}

bool NLO_Process::CalculateTotalXSec(const std::string &resultpath,
					const bool create) 
{ 
  p_int->Reset();
  SP(Phase_Space_Handler) psh(p_int->PSHandler());
  if (p_int->ISR()) {
    if (m_nin==2) {
      if (m_flavs[0].Mass()!=p_int->ISR()->Flav(0).Mass() ||
          m_flavs[1].Mass()!=p_int->ISR()->Flav(1).Mass()) {
        p_int->ISR()->SetPartonMasses(&m_flavs.front());
      }
    }
  }
  psh->InitCuts();
  if (p_int->ISR())
    p_int->ISR()->SetSprimeMin(psh->Cuts()->Smin());
  psh->CreateIntegrators();
  p_int->SetResultPath(resultpath);
  p_int->ReadResults();
  exh->AddTerminatorObject(p_int);
  psh->InitIncoming();
  double var(p_int->TotalVar());
  msg_Info()<<METHOD<<"(): Calculate xs for '"<<m_name<<"'"<<std::endl;
  double totalxs(psh->Integrate()/rpa.Picobarn());
  if (!IsEqual(totalxs,p_int->TotalResult())) {
    msg_Error()<<"Result of PS-Integrator and summation do not coincide!\n"
	       <<"  '"<<m_name<<"': "<<totalxs
	       <<" vs. "<<p_int->TotalResult()<<std::endl;
  }
  if (p_int->TotalXS()>0.0) {
    p_int->SetTotal();
    if (var==p_int->TotalVar()) {
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

void NLO_Process::SetLookUp(const bool lookup)
{
  if (p_bproc) p_bproc->SetLookUp(lookup);
  if (p_vproc) p_vproc->SetLookUp(lookup);
  if (p_iproc) p_iproc->SetLookUp(lookup);
  if (p_rproc) p_rproc->SetLookUp(lookup);
  if (p_sproc) p_sproc->SetLookUp(lookup);
}

void NLO_Process::SetScale(const std::string &scale)
{
  if (p_bproc) p_bproc->SetScale(scale);
  if (p_vproc) p_vproc->SetScale(scale);
  if (p_iproc) p_iproc->SetScale(scale);
  if (p_rproc) p_rproc->SetScale(scale);
  if (p_sproc) p_sproc->SetScale(scale);
}

void NLO_Process::SetKFactor(const std::string &kfactor,
			     const size_t &oqcdlo,const size_t &oewlo)
{
  if (p_bproc) p_bproc->SetKFactor(kfactor,oqcdlo,oewlo);
  if (p_vproc) p_vproc->SetKFactor(kfactor,oqcdlo,oewlo);
  if (p_iproc) p_iproc->SetKFactor(kfactor,oqcdlo,oewlo);
  if (p_rproc) p_rproc->SetKFactor(kfactor,oqcdlo,oewlo);
  if (p_sproc) p_sproc->SetKFactor(kfactor,oqcdlo,oewlo);
}

void NLO_Process::SetKFactorOn(const bool on)
{
  if (p_bproc) p_bproc->SetKFactorOn(on);
  if (p_vproc) p_vproc->SetKFactorOn(on);
  if (p_iproc) p_iproc->SetKFactorOn(on);
  if (p_rproc) p_rproc->SetKFactorOn(on);
  if (p_sproc) p_sproc->SetKFactorOn(on);
}

void NLO_Process::SetClusterDefinitions
(Cluster_Definitions_Base *const cluster)
{
  if (p_rproc) 
    for (size_t i(0);i<p_rproc->Size();++i)
      (*p_rproc)[i]->Generator()->SetClusterDefinitions(cluster);
}

void NLO_Process::SetSelector(const Selector_Key &key)
{
  if (p_bproc) p_bproc->SetSelector(key);
  if (p_vproc) p_vproc->SetSelector(key);
  if (p_iproc) p_iproc->SetSelector(key);
  if (p_sproc) p_sproc->SetSelector(key);
}

bool NLO_Process::Trigger(const ATOOLS::Vec4D_Vector &p)
{
  return p_bproc->Trigger(p);  
}

bool NLO_Process::FillIntegrator(Phase_Space_Handler *const psh)
{
  return p_bproc->FillIntegrator(psh);
}

void NLO_Process::UpdateIntegrator(Phase_Space_Handler *const psh)
{
  p_bproc->UpdateIntegrator(psh);
}

void NLO_Process::BuildCuts(Cut_Data *const cuts)
{
  return p_bproc->BuildCuts(cuts);
}

void NLO_Process::UpdateCuts(const double &sp,const double &y,
			      Cut_Data *const cuts)
{
  p_bproc->UpdateCuts(sp,y,cuts);
}

void NLO_Process::AddPoint(const double &value)
{
  p_int->EEG()->AddPoint(m_last-m_lastsingle);
}
