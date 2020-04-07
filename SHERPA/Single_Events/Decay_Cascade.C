#include "SHERPA/Single_Events/Decay_Cascade.H"

#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Phys/Blob.H"
#include "ATOOLS/Phys/Blob_List.H"
#include "ATOOLS/Phys/Momenta_Stretcher.H"
#include "ATOOLS/Phys/Cluster_Amplitude.H"
#include "ATOOLS/Math/Random.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Org/Return_Value.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Run_Parameter.H"
#include "ATOOLS/Org/Scoped_Settings.H"

#include "PHASIC++/Decays/Decay_Handler.H"
#include "PHASIC++/Decays/Decay_Clustering.H"

#include "METOOLS/SpinCorrelations/Spin_Density.H"
#include "METOOLS/SpinCorrelations/Decay_Matrix.H"
#include "METOOLS/SpinCorrelations/Amplitude2_Tensor.H"

#include "SHERPA/SoftPhysics/Soft_Photon_Handler.H"

using namespace SHERPA;
using namespace PHASIC;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;


Decay_Cascade::Decay_Cascade(const std::vector<PHASIC::Decay_Handler*>& decayhandlers, Soft_Photon_Handler* sph) :
  m_decayhandlers(decayhandlers), p_softphotons(sph)
{
  Settings& s = Settings::GetMainSettings();
  m_mass_smearing = s["DECAY_MASS_SMEARING"].SetDefault(2).Get<int>();
  m_qedmode = s["DECAY_QED_CORRECTIONS"].SetDefault(1).Get<int>();
  m_spincorr = s["DECAY_SPIN_CORRELATIONS"].SetDefault(1).Get<int>();
  m_name      = std::string("Decays");
  m_type      = eph::Perturbative;

  // TODO
  if (m_qedmode && m_spincorr) m_qedmode=2;
  /*
  if (m_qedmode>0 && !m_usemass) {
    THROW(fatal_error,std::string("QED corrections to hard decays only ")
                      +std::string("available in massive mode."));
  }
  */
}

Decay_Cascade::~Decay_Cascade() 
{
}

Return_Value::code Decay_Cascade::Treat(Blob_List * bloblist)
{
  if(bloblist->empty()) return Return_Value::Nothing;

  bool didit(false);
  p_bloblist=bloblist;
  for (size_t blit(0);blit<bloblist->size();++blit) {
    Blob* blob=(*bloblist)[blit];
    if (blob->Has(blob_status::needs_harddecays) || blob->Has(blob_status::needs_hadrondecays)) {
      DEBUG_FUNC("Treating blob "<<blob->Id());
      didit = true;
      try {
        METOOLS::Amplitude2_Tensor* amps(NULL);
        if (m_spincorr) {
          Blob* signal=bloblist->FindFirst(btp::Signal_Process);
          if (signal) {
            Blob_Data_Base* data = (*signal)["ATensor"];
            if (data) amps=data->Get<METOOLS::Amplitude2_Tensor*>();
          }
        }
        TreatInitialBlob(blob, amps);
      } catch (Return_Value::code ret) {
        return ret;
      }
      blob->UnsetStatus(blob_status::needs_harddecays);
      blob->UnsetStatus(blob_status::needs_hadrondecays);
    }
  }
  return (didit ? Return_Value::Success : Return_Value::Nothing);
}


void Decay_Cascade::TreatInitialBlob(ATOOLS::Blob* blob,
                                     METOOLS::Amplitude2_Tensor* amps)
{
  DEBUG_FUNC("");
  DEBUG_VAR(*blob);

  // TODO more elegant handling (send back to Event_Handler after removing blobs?)
  if (blob->NInP()==1 && blob->InParticle(0)->Flav().DecayHandler()) {
    // option to veto/re-do an existing decay blob, e.g. to reject exclusive
    // hadron decay channels from fragmentation
    if (blob->InParticle(0)->Flav().DecayHandler()->VetoDecayAndPrepForNew(blob)) {
      Spin_Density* sigma=new Spin_Density(blob->InParticle(0));
      FillDecayTree(blob, sigma);
      delete sigma;
      if (p_softphotons && m_qedmode==2) AttachExtraQEDRecursively(blob);
      for (auto dh: m_decayhandlers) dh->AfterTreatInitialBlob(blob, p_bloblist);
      return;
    }
  }

  // random shuffle, against bias in spin correlations and mixing
  Particle_Vector daughters = blob->GetOutParticles();
  std::vector<size_t> shuffled(daughters.size());
  for (size_t i=0; i<daughters.size(); ++i) shuffled[i]=i;
  for (size_t i=0; i<daughters.size(); ++i) {
    if (!daughters[i]->Flav().Stable() &&
	abs(daughters[i]->Momentum().Abs2()-
	    sqr(daughters[i]->FinalMass()))>1e-6) {
      PRINT_INFO("Initial particle "<<daughters[i]->Flav()<<" not onshell: "
                 <<"p^2="<<daughters[i]->Momentum().Mass()
                 <<" vs. m^2="<<daughters[i]->FinalMass());
      //      throw Return_Value::Retry_Event;
    }
  }
  random_shuffle(shuffled.begin(), shuffled.end(), *ran);

  // initial blobs still contain on-shell particles, stretch them off-shell
  for (size_t i=0; i<daughters.size(); ++i) {
    if (!daughters[shuffled[i]]->Flav().IsStable())
      daughters[shuffled[i]]->Flav().DecayHandler()->CreateDecayBlob(p_bloblist, daughters[shuffled[i]]);
  }
  if (SetMasses(blob, false)) {
    Momenta_Stretcher stretch;
    if (!stretch.StretchBlob(blob)) {
      msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl
                 <<*blob<<endl;
      throw Return_Value::Retry_Event;
    }
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
    if (amps && amps->Contains(daughters[i]->OriginalPart())) {
      Decay_Matrix* D(NULL);
      if (daughters[i]->Flav().IsStable()) {
        D=new Decay_Matrix(daughters[i]->OriginalPart()); // delta
      }
      else {
        if (!daughters[i]->DecayBlob() || daughters[i]->DecayBlob()->NOutP()>0) THROW(fatal_error,"Internal Error");
        DEBUG_INFO("treating: "<<daughters[i]->OriginalPart()->Flav());
        Spin_Density sigma(daughters[i]->OriginalPart(),amps);
        sigma.SetParticle(daughters[i]);
        DEBUG_VAR(sigma);
        D=FillDecayTree(daughters[i]->DecayBlob(), &sigma);
        D->SetParticle(daughters[i]->OriginalPart());
      }
      DEBUG_INFO("contracting with D["<<D->Particle()->Flav()<<"]");
      amps->Contract(D);
      delete D;
    }
    else {
      DEBUG_INFO("Did not find "<<daughters[i]->OriginalPart()->Flav()<<" in amps");
      Spin_Density sigma(daughters[i]);
      if (!daughters[i]->Flav().IsStable()) {
        FillDecayTree(daughters[i]->DecayBlob(), &sigma);
      }
    }
  }
  if (p_softphotons && m_qedmode==2) AttachExtraQEDRecursively(blob);

  for (auto dh: m_decayhandlers) dh->AfterTreatInitialBlob(blob, p_bloblist);

  UpdateShowerBlob(blob);
}





Decay_Matrix* Decay_Cascade::FillDecayTree(Blob * blob, Spin_Density* s0)
{
  blob->InParticle(0)->Flav().DecayHandler()->BeforeFillDecayTree(blob);
  
  Particle* inpart = blob->InParticle(0);
  DEBUG_FUNC(inpart->RefFlav()<<" "<<inpart->Number());
  if (s0) DEBUG_VAR(*s0);
  Vec4D labmom = inpart->Momentum();
  
  // fill decay blob all on-shell
  Blob_Data_Base* data = (*blob)["p_onshell"];
  if (data) inpart->SetMomentum(data->Get<Vec4D>());
  else {
    msg_Error()<<METHOD<<" could not find p_onshell, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  msg_Debugging()<<*blob<<std::endl;
  Amplitude2_Tensor* amps=inpart->Flav().DecayHandler()->FillOnshellDecay(blob, s0);
  inpart->SetStatus(part_status::decayed);
  if (inpart->Info()!='M') inpart->SetInfo('D');

  // Special case: interrupt decay cascade at partonic hadron decay
  if (blob->Type()==btp::Hadron_Decay) {
    bool partonic=false;
    for (auto p: blob->GetOutParticles()) if (p->Flav().IsQCD()) partonic=true;
    if (partonic) {
      DEBUG_INFO("Interrupt decay cascade at partonic hadron decay.");

      // TODO this just reconstructs what is done in master, but could be made more elegant?
      Particle_Vector daughters = blob->GetOutParticles();
      random_shuffle(daughters.begin(), daughters.end(), *ran);
      SetMasses(blob, true);
      BoostAndStretch(blob, labmom);
      if (p_softphotons) AttachExtraQED(blob);

      for (auto p: blob->GetOutParticles()) {
        if (amps) amps->ContractDiagonal(p);
      }
      return amps?new Decay_Matrix(inpart,amps):NULL;
    }
  }

  Particle_Vector daughters = blob->GetOutParticles();
  random_shuffle(daughters.begin(), daughters.end(), *ran);
  for (size_t i(0); i<daughters.size();++i) {
    if (daughters[i]->Flav().IsStable()) continue;
    daughters[i]->Flav().DecayHandler()->CreateDecayBlob(p_bloblist, daughters[i]);
  }

  SetMasses(blob, true);
  BoostAndStretch(blob, labmom);
  DEBUG_VAR(*blob);
  if (p_softphotons) AttachExtraQED(blob);

  DEBUG_INFO("recursively treating the created daughter decay blobs:");
  for (size_t i(0); i<daughters.size();++i) {
    // have to ignore photons from soft photon handler
    if (daughters[i]->Info()=='S') continue;
    DEBUG_VAR(daughters[i]->Flav());

    if (daughters[i]->Flav().IsStable()) {
      DEBUG_INFO("is stable");
      if (amps) amps->ContractDiagonal(daughters[i]);
      continue;
    }

    if (amps) {
      Spin_Density si(daughters[i],s0,amps);
      DEBUG_INFO("decaying with spin density "<<si);
      Decay_Matrix* Di=FillDecayTree(daughters[i]->DecayBlob(), &si);
      if (Di) {
        amps->Contract(Di);
        delete Di;
      }
      else {
        amps->ContractDiagonal(daughters[i]);
      }
    }
    else {
      DEBUG_INFO("decaying without spin correlations");
      FillDecayTree(daughters[i]->DecayBlob(), NULL);
      amps->ContractDiagonal(daughters[i]);
    }
  }
  DEBUG_INFO("finished daughters of "<<inpart->RefFlav()<<" "
             <<inpart->Number());

  return amps?new Decay_Matrix(inpart,amps):NULL;
}






void Decay_Cascade::UpdateShowerBlob(Blob* initialblob)
{
  DEBUG_FUNC(initialblob->Id());
  Blob_Data_Base * bdb((*initialblob)["ShowerBlob"]);
  if (!bdb) return;
  Blob* showerblob=bdb->Get<Blob*>();
  for (auto outpart: initialblob->GetOutParticles()) {
    if (outpart->Flav().IsStable() && outpart->DecayBlob() != showerblob) {
      THROW(fatal_error, "Inconsistent decay input blob.");
    }
    if (!outpart->Flav().Stable()) {
      AddDecayFinalState(showerblob, outpart);
    }
  }
  Decay_Clustering clustering;
  clustering.DefineInitialShowerConditions(initialblob, showerblob);
  DEBUG_VAR(*showerblob);
}



void Decay_Cascade::AddDecayFinalState(Blob* showerblob, Particle* part)
{
  if (part->DecayBlob()) {
    showerblob->RemoveInParticle(part, false);
    for (auto daughter: part->DecayBlob()->GetOutParticles()) {
      AddDecayFinalState(showerblob, daughter);
    }
  }
  else {
    showerblob->AddToInParticles(part);
  }
}



bool DecayWidthSortFunc(Particle* p1, Particle* p2) {
  double width1(0.0), width2(0.0);
  if (p1->Flav().DecayHandler()) width1=p1->Flav().DecayHandler()->DecayWidth(p1->Flav());
  if (p2->Flav().DecayHandler()) width2=p2->Flav().DecayHandler()->DecayWidth(p2->Flav());
  return width1 < width2;
}


bool Decay_Cascade::SetMasses(ATOOLS::Blob* blob, bool usefinalmass)
{
  DEBUG_FUNC(blob->GetOutParticles().size());
  Particle_Vector daughters = blob->GetOutParticles();
  if (m_mass_smearing==0 || daughters.size()==1) return false;
  sort(daughters.begin(), daughters.end(), DecayWidthSortFunc);

  double max_mass;
  if (usefinalmass) {
    max_mass=blob->InParticle(0)->FinalMass();
    for(auto daughter: daughters) {
      if (daughter->Flav().Stable()) max_mass -= daughter->FinalMass();
    }
  }
  else {
    Vec4D total(0.0,0.0,0.0,0.0);
    for(auto daughter: daughters) {
      if (!daughter->Flav().Stable()) total+=daughter->Momentum();
    }
    max_mass=total.Mass();
  }
  
  bool success=true; size_t cnt=0; bool didit=false;
  do {
    success=true;
    double max = max_mass;
    for(auto daughter: daughters) {
      if(!daughter->Flav().IsStable() && daughter->DecayBlob()) {
        if (!daughter->Flav().DecayHandler()->DiceMass(daughter,max)) {
          success=false; ++cnt;
          if (cnt>22) {
            msg_Error()<<METHOD<<" failed to set masses, retrying event."<<endl;
            throw Return_Value::Retry_Event;
          }
          break;
        }
        didit = true;
        max-=daughter->FinalMass();
      }
    }
  } while(success==false);
  return didit;
}

void Decay_Cascade::BoostAndStretch(Blob* blob, const Vec4D& labmom)
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
  Momenta_Stretcher stretch;
  if (!stretch.StretchBlob(blob)) {
    msg_Error()<<METHOD<<" failed to stretch blob, retrying event."<<endl
               <<*blob<<endl;
    throw Return_Value::Retry_Event;
  }
  for (size_t i(0); i<blob->NOutP(); ++i)
    if (blob->OutParticle(i)->DecayBlob())
      blob->OutParticle(i)->DecayBlob()
          ->AddData("p_actual",
                    new Blob_Data<Vec4D>(blob->OutParticle(i)->Momentum()));
  DEBUG_VAR(blob->MomentumConserved());
}



bool Decay_Cascade::AttachExtraQED(Blob* blob, size_t mode)
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

bool Decay_Cascade::AttachExtraQEDToProductionBlob(Blob* blob)
{
  DEBUG_FUNC("qedmode="<<m_qedmode<<", decay "<<blob->ShortProcessName());
  return false;
}

bool Decay_Cascade::AttachExtraQEDRecursively(Blob* blob, bool aa)
{
  DEBUG_FUNC("qedmode="<<m_qedmode<<", decay "<<blob->ShortProcessName()
             <<", already boosted="<<aa);
  if (m_qedmode!=2) return false;
  aa+=AttachExtraQED(blob,1);
  msg_Debugging()<<"added anything: "<<aa<<std::endl;
  for (auto part: blob->GetOutParticles()) {
    if (!part->Flav().IsStable() && part->DecayBlob()) {
      Blob * decblob(part->DecayBlob());
      if (aa) UpdateDecayBlob(decblob);
      AttachExtraQEDRecursively(decblob,aa);
    }
  }
  return aa;
}

void Decay_Cascade::UpdateDecayBlob(Blob* blob)
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

bool Decay_Cascade::CheckOnshellness(Blob* blob)
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
  Momenta_Stretcher momstretch;
  momstretch.StretchMomenta(blob->GetOutParticles(),masses);
  return false;
}


void Decay_Cascade::CleanUp(const size_t & mode)
{
  for (size_t i=0; i<m_decayhandlers.size(); ++i) {
    m_decayhandlers[i]->CleanUp();
  }
}

void Decay_Cascade::Finish(const std::string &)
{
}
