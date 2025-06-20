#include "SHERPA/Single_Events/Decay_Handler_Base.H"

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Run_Parameter.H"

#include "PHASIC++/Decays/Decay_Channel.H"
#include "PHASIC++/Decays/Decay_Table.H"
#include "PHASIC++/Decays/Decay_Map.H"

#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Decay_Handler_Base::Decay_Handler_Base() :
  p_softphotons(NULL), p_polarization_handler(NULL), p_decaymap(NULL), p_bloblist(NULL), p_ampl(NULL),
  m_decaymatrices(std::vector<METOOLS::Decay_Matrix>()),
  m_stretcher(Momenta_Stretcher("Decay_Handler")),
  m_qedmode(0), m_spincorr(false), m_polcrosssec(false),
  m_decaychainend(false), m_cluster(true), m_mass_smearing(1), m_oserrors(0)
{
  auto& s = Settings::GetMainSettings();
  m_specialtauspincorr = s["SPECIAL_TAU_SPIN_CORRELATIONS"].SetDefault(0).Get<size_t>();
}

Decay_Handler_Base::~Decay_Handler_Base()
{
  if (m_oserrors>0)
    msg_Error()<<METHOD<<" with "<<m_oserrors<<" particles not on their mass shell.\n";
  if (p_decaymap) delete p_decaymap; p_decaymap=NULL;
  if (p_polarization_handler) { delete p_polarization_handler;  p_polarization_handler = NULL; }
}

class Decay_Width_Sorter {
  PHASIC::Decay_Map* p_decaymap;
public:
  Decay_Width_Sorter(PHASIC::Decay_Map* decaymap) : p_decaymap(decaymap) {}

  bool operator()(Particle* p1, Particle* p2) {
    Decay_Table* table1=p_decaymap->FindDecay(p1->Flav());
    Decay_Table* table2=p_decaymap->FindDecay(p2->Flav());
    double width1(0.0), width2(0.0);
    if (table1) width1=table1->TotalWidth();
    if (table2) width2=table2->TotalWidth();
    return width1 < width2;
  }
};

void Decay_Handler_Base::SetMasses(ATOOLS::Blob* blob, bool usefinalmass)
{
  DEBUG_FUNC(blob->GetOutParticles().size());
  Particle_Vector daughters = blob->GetOutParticles();
  if (m_mass_smearing==0 || daughters.size()==1) return;
  sort(daughters.begin(), daughters.end(), Decay_Width_Sorter(p_decaymap));

  Vec4D total(0.0,0.0,0.0,0.0);
  size_t nr_daughters(0);
  for(size_t i=0; i<daughters.size(); i++) {
    if (!daughters[i]->DecayBlob() || daughters[i]->DecayBlob()->NOutP()==0) {
      total+=daughters[i]->Momentum();
      ++nr_daughters;
    }
  }
  double max_mass;
  if (usefinalmass) max_mass=blob->InParticle(0)->FinalMass();
  else max_mass=total.Mass();
  if (nr_daughters<2) return;

  bool success=true; size_t cnt=0;
  do {
    success=true;
    double max = max_mass;
    for(PVIt it=daughters.begin();it!=daughters.end();++it) {
      if(m_mass_smearing==2 && !Decays((*it)->Flav())) continue;
      if ((*it)->DecayBlob() && (*it)->DecayBlob()->NOutP()>0) continue;
      else if ((*it)->DecayBlob()) {
        if (!DiceMass(*it,max)) {
          success=false; ++cnt;
	  if (cnt>22) {
            msg_Error()<<METHOD<<" failed to set masses, retrying event."<<endl;
            throw Return_Value::Retry_Event;
          }
          break;
        }
      }
      else {
        double mass = (*it)->RefFlav().RelBWMass(0.0, max,
						 this->Mass((*it)->RefFlav()));
        (*it)->SetFinalMass(mass);
        DEBUG_INFO(max<<" > "<<"m["<<(*it)->RefFlav()<<"]="<<mass);
      }
      max-=(*it)->FinalMass();
    }
  } while(success==false);
}

bool Decay_Handler_Base::DiceMass(ATOOLS::Particle* p, double max)
{
  Blob* decayblob=p->DecayBlob();
  Blob_Data_Base* data = (*decayblob)["dc"];
  if (data) {
    Decay_Channel* dc = data->Get<Decay_Channel*>();
    if (!dc) THROW(fatal_error,"Missing decay channel for "
                               +decayblob->ShortProcessName()+".");
    double width = p_decaymap->FindDecay(p->Flav())->TotalWidth();
    double mass=dc->GenerateMass(max, width);
    if (mass>0.0) p->SetFinalMass(mass);
    else return false;
  }
  return true;
}

void Decay_Handler_Base::BoostAndStretch(Blob* blob, const Vec4D& labmom)
{
  DEBUG_FUNC("");
  DEBUG_VAR(blob->MomentumConserved());
  // 1.
  Particle* inpart = blob->InParticle(0);
  Vec4D mom(inpart->Momentum());
  double m02=sqr(inpart->FinalMass());
  double p02=mom.PSpat2();
  double E02=sqr(mom[0]);
  double factor=sqrt((m02+p02)/E02);
  DEBUG_VAR(factor);
  mom[0]*=factor;
  inpart->SetMomentum(mom);
  Particle_Vector daughters = blob->GetOutParticles();
  for(PVIt it=daughters.begin();it!=daughters.end();++it) {
    mom=(*it)->Momentum();
    mom[0]*=factor;
    (*it)->SetMomentum(mom);
  }
  DEBUG_VAR(blob->MomentumConserved());

  // 2.
  Poincare twiddle2rest(inpart->Momentum());
  Poincare labboost(labmom);
  labboost.Invert();

  blob->Boost(twiddle2rest);
  blob->Boost(labboost);
  DEBUG_VAR(blob->MomentumConserved());

  // 3.
  if (!m_stretcher.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob ("<<blob->Type()<<", "
	       <<blob->NInP()<<" -> "<<blob->NOutP()<<"), retrying event.\n";
    msg_Tracking()<<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  for (size_t i(0); i<blob->NOutP(); ++i)
    if (blob->OutParticle(i)->DecayBlob())
      blob->OutParticle(i)->DecayBlob()
          ->AddData("p_actual",
                    new Blob_Data<Vec4D>(blob->OutParticle(i)->Momentum()));
  DEBUG_VAR(blob->MomentumConserved());
}

Blob* FindSPBlob(Blob* startblob)
{
  if (startblob->Type()==btp::Signal_Process) {
    return startblob;
  }

  for (size_t i=0; i<startblob->NInP(); ++i) {
    if (startblob->InParticle(i)->ProductionBlob()) {
      Blob* blob = FindSPBlob(startblob->InParticle(i)->ProductionBlob());
      if (blob) return blob;
    }
  }
  return NULL;
}


typedef std::vector<std::pair<std::pair<ATOOLS::Flavour,ATOOLS::Vec4D>, Spin_Density*> > SpinDensityMap;
bool Decay_Handler_Base::DoSpecialDecayTauSC(Particle* part)
{
  if (!m_specialtauspincorr) return false;
  Blob* blob=part->ProductionBlob();
  if (blob==NULL || blob->Type()!=btp::Fragmentation) return false;
  for (size_t i=0; i<blob->NOutP(); ++i)
    if (blob->OutParticle(i)->Flav().Kfcode()!=kf_tau)
      return false;
  DEBUG_FUNC(*part);

  Blob* signal=FindSPBlob(blob);
  if (!signal) {
    PRINT_INFO("Signal blob not found.");
    return false;
  }
  Blob_Data_Base* data = (*signal)["Tau_SpinDensity"];
  SpinDensityMap* tau_spindensity = data ? data->Get<SpinDensityMap*>() : NULL;
  if (!tau_spindensity) return false;

  double bestDeltaR=1000.0; Spin_Density* sigma_tau=NULL;
  for (SpinDensityMap::iterator it=tau_spindensity->begin(); it!=tau_spindensity->end(); ++it) {
    if (it->first.first==part->Flav()) {
      double newDeltaR=part->Momentum().DR(it->first.second);
      if (newDeltaR<bestDeltaR) {
        bestDeltaR=newDeltaR;
        sigma_tau=it->second;
      }
    }
  }
  if (sigma_tau==NULL) {
    PRINT_INFO("Tau Spin_Density not found");
  }
  else {
    DEBUG_VAR(*sigma_tau);
    sigma_tau->SetParticle(part);
    Decay_Matrix* D=FillDecayTree(part->DecayBlob(), sigma_tau);
    delete D;
    return true;
  }
  return false;
}

void Decay_Handler_Base::TreatInitialBlob(ATOOLS::Blob* blob,
                                          METOOLS::Amplitude2_Tensor* amps,
                                          const Particle_Vector& origparts)
{
  DEBUG_FUNC("");
  DEBUG_VAR(*blob);
  m_decaychainend=false;
  // delete decay matrices from the previous event
  if (m_polcrosssec && !m_decaymatrices.empty()) m_decaymatrices=std::vector<METOOLS::Decay_Matrix>();
  // random shuffle, against bias in spin correlations and mixing
  Particle_Vector daughters = blob->GetOutParticles();
  std::vector<size_t> shuffled(daughters.size());
  for (size_t i=0; i<daughters.size(); ++i) shuffled[i]=i;
  for (size_t i=0; i<daughters.size(); ++i) {
    if (!daughters[i]->Flav().Stable() &&
	abs(daughters[i]->Momentum().Abs2()-
	    sqr(daughters[i]->FinalMass()))>1e-6) {
      if ((m_oserrors++)<5) {
	PRINT_INFO("Initial particle "<<daughters[i]->Flav()<<" not onshell: "
		   <<"sqrt|p^2|="<<sqrt(daughters[i]->Momentum().Abs2())
		   <<" vs. m="<<daughters[i]->FinalMass());
      }
      throw Return_Value::Retry_Event;
    }
  }

  std::shuffle(shuffled.begin(), shuffled.end(), *ran);
  // initial blobs still contain on-shell particles, stretch them off-shell
  for (size_t i=0; i<daughters.size(); ++i) {
    if (!Decays(daughters[shuffled[i]]->Flav()) ||
        daughters[shuffled[i]]->DecayBlob()) {
      continue;
    }
    CreateDecayBlob(daughters[shuffled[i]]);
  }
  SetMasses(blob, false);
  if (!m_stretcher.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob ("<<blob->Type()<<", "
	       <<blob->NInP()<<" -> "<<blob->NOutP()<<"), retrying event.\n";
    msg_Tracking()<<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  for (size_t i(0); i<blob->NOutP(); ++i)
    if (blob->OutParticle(i)->DecayBlob())
      blob->OutParticle(i)->DecayBlob()
          ->AddData("p_actual",
                    new Blob_Data<Vec4D>(blob->OutParticle(i)->Momentum()));
  DEBUG_VAR(*blob);

  for (size_t ii(0); ii<daughters.size(); ++ii) {
    size_t i=shuffled[ii];
    DEBUG_INFO("treating "<<*daughters[i]);
    m_decaychainend=false;
    if (m_spincorr) {
      if (amps && origparts[i]) {
        DEBUG_VAR(*amps);
        Decay_Matrix* D(NULL);
        if (!Decays(daughters[i]->Flav()) || !daughters[i]->DecayBlob() ||
            daughters[i]->DecayBlob()->NOutP()>0) {
          D=new Decay_Matrix(origparts[i]); // delta
        }
        else {
          DEBUG_INFO("treating: "<<origparts[i]->Flav());
          Spin_Density sigma(origparts[i],amps);
          sigma.SetParticle(daughters[i]);
          DEBUG_VAR(sigma);
          D=FillDecayTree(daughters[i]->DecayBlob(), &sigma);
          D->SetParticle(origparts[i]);

          // save all decay matrices for later transformation to desired bases of polarization definition
          // and calculation of polarized cross sections
          if (m_polcrosssec) {
            // TODO: What happens if vector bosons decay into unstable particles?
            m_decaymatrices.push_back(*D);
          }
        }
        if (amps->Contains(origparts[i]) &&
            // Contract tau spins here only if they are globally stable or decayed here,
            // i.e. not in hard decays if taus will be decayed in hadron decays
            // (this is relevant e.g. for tttautau)
            (origparts[i]->Flav().Kfcode()!=kf_tau ||
             Flavour(kf_tau).IsStable() ||
             Decays(Flavour(kf_tau)))) {
          DEBUG_INFO("contracting with D["<<D->Particle()<<"]");
          amps->Contract(D);
        }
        delete D;
      }
      else {
        Spin_Density sigma(daughters[i]);
        if (Decays(daughters[i]->Flav())) {
          if (DoSpecialDecayTauSC(daughters[i])) {
            DEBUG_INFO("did special tau spin correlation treatment");
          }
          else {
            Decay_Matrix* D=FillDecayTree(daughters[i]->DecayBlob(), &sigma);
            delete D;
          }
        }
      }
    }
    else {
      if (Decays(daughters[i]->Flav())) {
        FillDecayTree(daughters[i]->DecayBlob(), NULL);
      }
    }
    m_decaychainend=true;
  }
  if (p_softphotons && m_qedmode==2)
    AttachExtraQEDRecursively(blob);
}

Decay_Matrix* Decay_Handler_Base::FillDecayTree(Blob * blob, Spin_Density* s0)
{
  Particle* inpart = blob->InParticle(0);
  DEBUG_FUNC(inpart->RefFlav()<<" "<<inpart->Number());
  if (s0) DEBUG_VAR(*s0);
  Vec4D labmom = inpart->Momentum();

  // fill decay blob all on-shell
  Blob_Data_Base* data = (*blob)["p_onshell"];
  if (data) inpart->SetMomentum(data->Get<Vec4D>());
  else {
    msg_Error()<<METHOD<<" could not find p_onshell tag in ("<<blob->Type()<<", "
	       <<blob->NInP()<<" -> "<<blob->NOutP()<<"), retrying event.\n";
    msg_Tracking()<<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  msg_Debugging()<<*blob<<std::endl;
  Amplitude2_Tensor* amps=FillOnshellDecay(blob, s0);
  inpart->SetStatus(part_status::decayed);
  if (inpart->Info()!='M') inpart->SetInfo('D');

  Particle_Vector daughters = blob->GetOutParticles();
  std::shuffle(daughters.begin(), daughters.end(), *ran);

  if (!(blob->Type()==btp::Hadron_Decay &&
        blob->Has(blob_status::needs_showers))) {
    for (PVIt it=daughters.begin();it!=daughters.end();++it) {
      if (Decays((*it)->Flav())) {
        if (!CanDecay(inpart->Flav())) {
          msg_Error()<<METHOD<<" Particle '"<<inpart->Flav()
                     <<"' set unstable, but decay handler doesn't know how "
                     <<"to deal with it.";
          throw Return_Value::Retry_Event;
        }
        CreateDecayBlob(*it);
      }
    }
  }

  SetMasses(blob, true);
  BoostAndStretch(blob, labmom);
  DEBUG_VAR(*blob);
  if (p_softphotons) AttachExtraQED(blob);

  DEBUG_INFO("recursively treating the created daughter decay blobs:");
  if (m_spincorr) DEBUG_VAR(*amps);
  for (size_t i(0); i<daughters.size();++i) {
    // have to ignore photons from soft photon handler
    if (daughters[i]->Info()=='S') continue;
    DEBUG_VAR(daughters[i]->Flav());

    if (!Decays(daughters[i]->Flav()) ||
        (blob->Type()==btp::Hadron_Decay &&
         blob->Has(blob_status::needs_showers))) {
      DEBUG_INFO("is stable.");
      if (m_specialtauspincorr && daughters[i]->Flav().Kfcode()==kf_tau &&
          !daughters[i]->Flav().IsStable() && blob->Type()==btp::Hard_Decay &&
          rpa->gen.SoftSC()) {
        DEBUG_INFO("  keeping tau spin information for hadronic tau decays.");
        SpinDensityMap* tau_spindensity;
        Blob* spblob(FindSPBlob(blob));
        if (!spblob) THROW(fatal_error, "Internal Error 1");
        Blob_Data_Base * bdb((*spblob)["Tau_SpinDensity"]);
        if (!bdb) {
          tau_spindensity = new SpinDensityMap;
          spblob->AddData("Tau_SpinDensity",new Blob_Data<SpinDensityMap*>(tau_spindensity));
        }
        else {
          tau_spindensity = bdb->Get<SpinDensityMap*>();
        }
        tau_spindensity->push_back(make_pair(make_pair(daughters[i]->Flav(), daughters[i]->Momentum()),
                                             new Spin_Density(daughters[i],amps)));
        DEBUG_VAR(*(tau_spindensity->back().second));
        tau_spindensity->back().second->SetParticle(NULL);
      }
      else if (m_spincorr) {
        Decay_Matrix* D=new Decay_Matrix(daughters[i]);
        amps->Contract(D);
        delete D;
      }
      if (m_spincorr) DEBUG_VAR(*amps);
    }
    else {
      if (!CanDecay(inpart->Flav())) {
	msg_Error()<<METHOD<<" could not decay particle flavour "<<inpart->Flav()
		   <<" in ("<<blob->Type()<<", "
		   <<blob->NInP()<<" -> "<<blob->NOutP()<<"), retrying event.\n";
	msg_Tracking()<<*blob<<endl;
      }
      Spin_Density* si=m_spincorr?new Spin_Density(daughters[i],s0,amps):NULL;
      DEBUG_INFO("decaying with spin density "<<si);
      Decay_Matrix* Di=FillDecayTree(daughters[i]->DecayBlob(), si);
      if (si) delete si;
      if (Di) {
        amps->Contract(Di);
        delete Di;
      }
    }
  }
  DEBUG_INFO("finished daughters of "<<inpart->RefFlav()<<" "
             <<inpart->Number());
  if (m_spincorr) DEBUG_VAR(*amps);
  return m_spincorr?new Decay_Matrix(inpart,amps):NULL;
}

Amplitude2_Tensor* Decay_Handler_Base::FillOnshellDecay(Blob *blob,
                                                        Spin_Density* sigma)
{
  DEBUG_FUNC("");
  Decay_Channel* dc(NULL);
  Blob_Data_Base* data = (*blob)["dc"];
  if (data) {
    dc=data->Get<Decay_Channel*>();
  }
  else {
    Decay_Table* table=p_decaymap->FindDecay(blob->InParticle(0)->Flav());
    if (table==NULL) {
      msg_Error()<<METHOD<<"Error: Did not find "<<blob->InParticle(0)->Flav()
                 <<" in decay tables."<<endl;
      throw Return_Value::Retry_Event;
    }
    dc=table->Select();
    blob->AddData("dc",new Blob_Data<Decay_Channel*>(dc));
  }
  if (!dc) THROW(fatal_error,"No decay channel found for "
                             +blob->InParticle(0)->Flav().IDName()+".");
  msg_Debugging()<<*dc<<std::endl;

  Particle* inpart=blob->InParticle(0);
  inpart->SetStatus(part_status::decayed);
  Flavour flav; Particle* particle=NULL;
  for (size_t i=1; i<dc->Flavs().size(); ++i) {
    flav=dc->Flavs()[i];
    if (inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti()) flav = flav.Bar();
    particle = new Particle(0, flav);
    particle->SetFinalMass(Mass(flav));
    particle->SetStatus(part_status::active);
    particle->SetNumber();
    particle->SetInfo('D');
    blob->AddToOutParticles( particle );
  }

  size_t n=1+blob->NOutP();
  vector<Vec4D> moms(n);
  moms[0]=inpart->Momentum();

  Amplitude2_Tensor* amps(NULL);
  if (sigma) {
    std::vector<Particle*> parts;
    parts.insert(parts.end(), blob->InParticle(0));
    Particle_Vector outparts=blob->GetOutParticles();
    parts.insert(parts.end(), outparts.begin(), outparts.end());
    dc->GenerateKinematics(moms,
                           inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti(),
                           sigma,parts);
    amps=dc->Amps();
  }
  else {
    dc->GenerateKinematics(moms,
                           inpart->Flav().IsAnti()!=dc->GetDecaying().IsAnti());
  }
  for (size_t i=1; i<n; i++) blob->GetOutParticles()[i-1]->SetMomentum(moms[i]);
  msg_Debugging()<<*blob<<std::endl;
  return amps;
}

Cluster_Amplitude* Decay_Handler_Base::ClusterConfiguration(Blob *const bl)
{
  msg_Indent();
  p_ampl = Cluster_Amplitude::New();
  p_ampl->SetMS(this);
  if (!(bl->NInP()==1 && bl->Type()==btp::Hadron_Decay &&
	bl->InParticle(0)->GetFlow(1)==0 && bl->InParticle(0)->GetFlow(2)==0) ) {
    for (int i(0);i<bl->NInP();++i) {
      Particle *p(bl->InParticle(i));
      ColorID col(p->GetFlow(2),p->GetFlow(1));
      p_ampl->CreateLeg(-p->Momentum(),p->Flav().Bar(),col,1<<i);
    }
    p_ampl->SetNIn(bl->NInP());
  }
  for (int i(0);i<bl->NOutP();++i) {
    Particle *p(bl->OutParticle(i));
    if (p->GetFlow(1)==0 && p->GetFlow(2)==0) continue;
    ColorID col(p->GetFlow(1),p->GetFlow(2));
    p_ampl->CreateLeg(p->Momentum(),p->Flav(),col,1<<(i+p_ampl->NIn()));
  }
  while (m_cluster && p_ampl->Legs().size()>p_ampl->NIn()+2) {
    msg_Debugging()<<*p_ampl<<"\n";
    Cluster_Amplitude *ampl(p_ampl);
    p_ampl = p_ampl->InitNext();
    p_ampl->SetMS(this);
    for (size_t i(0);i<ampl->NIn();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
    p_ampl->SetNIn(ampl->NIn());
    Cluster_Leg *lij(NULL);
    for (size_t i(ampl->NIn());i<ampl->Legs().size()-1;++i) {
      Cluster_Leg *li(ampl->Leg(i));
      for (size_t j(i+1);j<ampl->Legs().size();++j) {
        Cluster_Leg *lj(ampl->Leg(j));
        ColorID nc;
        if (li->Col().m_i==0 && li->Col().m_j==0) {
          nc=lj->Col();
        }
        else if (lj->Col().m_i==0 && lj->Col().m_j==0) {
          nc=li->Col();
        }
        else if (li->Col().m_i && li->Col().m_i==lj->Col().m_j) {
          nc.m_i=lj->Col().m_i;
          nc.m_j=li->Col().m_j;
        }
        else if (li->Col().m_j && li->Col().m_j==lj->Col().m_i) {
          nc.m_i=li->Col().m_i;
          nc.m_j=lj->Col().m_j;
        }
        if (nc.m_i>=0 && nc.m_j>=0) {
          Flavour fl(kf_photon);
          if (nc.m_i && nc.m_j) fl=Flavour(kf_gluon);
          else if (nc.m_i) fl=Flavour(kf_d);
          else if (nc.m_j) fl=Flavour(kf_d).Bar();
          p_ampl->CreateLeg(li->Mom()+lj->Mom(),fl,nc,li->Id()+lj->Id());
	  size_t k(0);
	  for (k=ampl->NIn();k<ampl->Legs().size();++k)
	    if (k!=i && k!=j) break;
          p_ampl->Legs().back()->SetK(ampl->Leg(k)->Id());
          p_ampl->Prev()->SetIdNew(lj->Id());
          lij=p_ampl->Legs().back();
          break;
        }
      }
      if (lij) break;
    }
    if (lij==NULL) THROW(fatal_error,"Internal error");
    for (size_t i(ampl->NIn());i<ampl->Legs().size();++i) {
      Cluster_Leg *cl(ampl->Leg(i));
      if (cl->Id()&lij->Id()) continue;
      p_ampl->CreateLeg(cl->Mom(),cl->Flav(),cl->Col(),cl->Id());
    }
  }
  double mu2=p_ampl->Leg(0)->Mom().Abs2();
  p_ampl->SetMuF2(mu2);
  p_ampl->SetKT2(mu2);
  p_ampl->SetMuQ2(mu2);
  msg_Debugging()<<*p_ampl<<"\n";
  while (p_ampl->Prev()) {
    p_ampl=p_ampl->Prev();
    p_ampl->SetMuF2(mu2);
    p_ampl->SetKT2(mu2);
    p_ampl->SetMuQ2(mu2);
  }
  msg_Debugging()<<"}\n";
  return p_ampl;
}

void Decay_Handler_Base::CleanUp()
{
  if (p_decaymap) p_decaymap->ResetCounters();
}

bool Decay_Handler_Base::Decays(const ATOOLS::Flavour& flav)
{
  if (flav.IsStable()) return false;
  return true;
}

bool Decay_Handler_Base::CanDecay(const ATOOLS::Flavour& flav)
{
  if (p_decaymap) return p_decaymap->Knows(flav);
  else return false;
}

bool Decay_Handler_Base::AttachExtraQED(Blob* blob, size_t mode)
{
  DEBUG_FUNC("qedmode="<<m_qedmode
             <<", shower="<<blob->Has(blob_status::needs_showers)
             <<", qed="<<blob->Has(blob_status::needs_extraQED)
             <<", mode="<<mode
             <<", process="<<blob->ShortProcessName());
  if (!blob->Has(blob_status::needs_extraQED)) return false;
  if (blob->NInP()!=1) return AttachExtraQEDToProductionBlob(blob);
  if (mode==0 && m_qedmode!=1) return false;
  if (mode==1 && m_qedmode!=2) return false;
  for (size_t i(0);i<blob->NOutP();++i)
    if (blob->OutParticle(i)->Flav().Strong()) return false;
  msg_Debugging()<<*blob<<std::endl;
  msg_Debugging()<<"Momentum conserved: "<<blob->CheckMomentumConservation()
                 <<std::endl;
  if (!p_softphotons->AddRadiation(blob)) {
    msg_Error()<<METHOD<<"(): Soft photon handler failed, retrying event."
               <<std::endl;
    throw Return_Value::Retry_Event;
  }
  msg_Debugging()<<*blob<<std::endl;
  msg_Debugging()<<"Momentum conserved: "<<blob->CheckMomentumConservation()
                 <<std::endl;
  blob->UnsetStatus(blob_status::needs_extraQED);
  msg_Debugging()<<"Added anything? "<<p_softphotons->AddedAnything()
                 <<std::endl;
  return p_softphotons->AddedAnything();
}

bool Decay_Handler_Base::AttachExtraQEDToProductionBlob(Blob* blob)
{
  DEBUG_FUNC("qedmode="<<m_qedmode<<", decay "<<blob->ShortProcessName());
  return false;
}

bool Decay_Handler_Base::AttachExtraQEDRecursively(Blob* blob, bool aa)
{
  DEBUG_FUNC("qedmode="<<m_qedmode<<", decay "<<blob->ShortProcessName()
             <<", already boosted="<<aa);
  if (m_qedmode!=2) return false;
  aa+=AttachExtraQED(blob,1);
  msg_Debugging()<<"added anything: "<<aa<<std::endl;
  for (size_t i(0);i<blob->NOutP();++i) {
    if (blob->OutParticle(i)->DecayBlob()) {
      Blob * decblob(blob->OutParticle(i)->DecayBlob());
      msg_Debugging()<<blob->OutParticle(i)->Flav()<<" has "
                     <<(blob->OutParticle(i)->DecayBlob()?"a ":"no ")
                     <<"decay blob"<<std::endl;
      if (decblob) {
        if (aa) UpdateDecayBlob(decblob);
        AttachExtraQEDRecursively(decblob,aa);
      }
    }
  }
  return aa;
}

void Decay_Handler_Base::UpdateDecayBlob(Blob* blob)
{
  DEBUG_FUNC(blob->ShortProcessName());
  const Vec4D& P((*blob)["p_actual"]->Get<Vec4D>());
  const Vec4D& Pt(blob->InParticle(0)->Momentum());
  const Vec4D e(P-Pt);
  msg_Debugging()<<"P-Pt="<<e<<" ["<<e.Mass()<<"]"<<std::endl;
  const Lorentz_Ten2D lambda(MetricTensor()-2.*BuildTensor(e,e)/e.Abs2());
  msg_Debugging()<<"\\Lambda="<<std::endl<<lambda<<std::endl;
  for (size_t i(0);i<blob->NOutP();++i) {
    Vec4D mom(blob->OutParticle(i)->Momentum());
    msg_Debugging()<<blob->OutParticle(i)->Flav().IDName()<<" "
                   <<mom<<" ["<<mom.Mass()<<"] -> ";
    mom=Contraction(lambda,2,mom);
    blob->OutParticle(i)->SetMomentum(mom);
    msg_Debugging()<<mom<<" ["<<mom.Mass()<<"]"<<std::endl;
  }
  CheckOnshellness(blob);
  if (msg_LevelIsDebugging()) {
    for (size_t i(0);i<blob->NOutP();++i) {
      Vec4D mom(blob->OutParticle(i)->Momentum());
      msg_Debugging()<<blob->OutParticle(i)->Flav().IDName()<<" "
                     <<mom<<" ["<<mom.Mass()<<"]"<<std::endl;
    }
  }
  msg_Debugging()<<"Momentum conservation in decay blob of "
                 <<blob->InParticle(0)->Flav()<<": "
                 <<blob->CheckMomentumConservation()<<std::endl;
}

bool Decay_Handler_Base::CheckOnshellness(Blob* blob)
{
  std::vector<double> masses;
  bool allonshell(true);
  double accu(sqrt(Accu()));
  for (size_t i(0);i<blob->NOutP();++i) {
    masses.push_back(blob->OutParticle(i)->FinalMass());
    if (allonshell &&
        !IsEqual(blob->OutParticle(i)->Momentum().Abs2(),
                 sqr(blob->OutParticle(i)->FinalMass()),accu)) allonshell=false;
  }
  msg_Debugging()<<"masses="<<masses<<std::endl;
  if (allonshell) return true;
  msg_Debugging()<<"need to put on-shell"<<std::endl;
  m_stretcher.StretchMomenta(blob->GetOutParticles(),masses);
  return false;
}
